"""Reuse propagated orbits while evaluating independent design candidates."""
from __future__ import annotations

from dataclasses import asdict, replace
import hashlib
import json
import math
from pathlib import Path

import numpy as np
import pandas as pd

from .coverage import sso_inclination_deg
from .mission import EARTH_RADIUS_M, J2, J3, SelectedMission, analyze_accesses, propagate_selected_mission, target_grid
from .mission_clouds import apply_markov_clouds
from .orekit_runtime import resolve_data_path


class MissionCoverageEvaluator:
    """Evaluate hardware candidates without loading an adopted-design CSV.

    Orbit propagation depends on altitude, satellite phasing and timing, not
    detector choice, telescope count or swath. Those candidates reuse the same
    cached ephemeris and Sun positions. Only the most recent orbit is held in
    memory; compressed disk caches stay outside the repository by default.
    """

    def __init__(self, *, data_path: str | Path | None = None,
                 cache_dir: str | Path = "/private/tmp/tfg-orekit/coverage-cache",
                 base_mission: SelectedMission | None = None,
                 cloud_fraction: float = 1/6, persistence: float = 0.5,
                 mc_samples: int = 200, seed: int = 12345,
                 requirement_days: float = 7.0, backend: str = "orekit",
                 criterion: str = "p95"):
        self.base = base_mission or SelectedMission(duration_days=30)
        self.base.validate()
        if backend not in {"orekit", "scipy-reference"}:
            raise ValueError("Unknown propagation backend.")
        if criterion not in {"matlab-mean", "p95"}:
            raise ValueError("Unknown selection criterion")
        self.criterion = criterion
        self.backend, self.data_path = backend, data_path
        self.cache_dir = Path(cache_dir)
        self.cloud_parameters = dict(cloud_fraction=cloud_fraction, persistence=persistence,
                                     mc_samples=mc_samples, seed=seed, requirement_days=requirement_days)
        self._last_orbit = None
        self._last_data = None
        self._results: dict[str, dict] = {}

    def _orbit(self, mission: SelectedMission):
        inputs = {name: getattr(mission, name) for name in (
            "epoch", "altitude_km", "inclination_deg", "eccentricity", "argument_of_perigee_deg",
            "satellites", "ltan_h", "duration_days", "step_s")}
        inputs.update({"backend": self.backend, "cache_version": "j2-j3-mean-v1"})
        if self.backend == "orekit":
            data = resolve_data_path(self.data_path)
            inputs["physical_data"] = str(data)
            inputs["physical_data_mtime_ns"] = data.stat().st_mtime_ns
        key = hashlib.sha256(json.dumps(inputs, sort_keys=True).encode()).hexdigest()[:24]
        if key == self._last_orbit:
            return self._last_data
        self.cache_dir.mkdir(parents=True, exist_ok=True)
        path = self.cache_dir / f"orbit_{key}.npz"
        if path.exists():
            with np.load(path, allow_pickle=False) as saved:
                frame = pd.DataFrame(saved["ephemeris"], columns=saved["columns"].tolist())
                frame["satellite"] = frame.satellite.astype(int)
                frame["sunlit"] = frame.sunlit.astype(bool)
                sun = saved["sun"]
                provenance = json.loads(str(saved["provenance"]))
        else:
            if self.backend == "orekit":
                frame, sun, provenance = propagate_selected_mission(mission, self.data_path)
            else:
                from .mission_reference import propagate_reference
                frame, sun, provenance = propagate_reference(mission)
            np.savez_compressed(path, ephemeris=frame.to_numpy(dtype=float),
                                columns=np.asarray(frame.columns, dtype=str), sun=sun,
                                provenance=json.dumps(provenance))
        self._last_orbit, self._last_data = key, (frame, sun, provenance)
        return self._last_data

    def evaluate_probe(self, altitude_km: float, satellite_count: int, acquired_swath_km: float) -> dict:
        """A failed subset is a lower-bound rejection; a pass still needs the full grid."""
        return self.evaluate(altitude_km, satellite_count, acquired_swath_km, probe=True)

    def evaluate(self, altitude_km: float, satellite_count: int, acquired_swath_km: float, *, probe: bool = False) -> dict:
        """Return geometric and Markov-cloud metrics for total acquired width.

        The caller must enforce detector/FoV/mass limits and supply the rounded
        total acquired swath. Telescope count must not multiply this width again.
        For identical cloud realizations, widening the footprint cannot worsen
        a no-access-gap statistic; candidate searches may use that monotonicity.
        """
        inclination = sso_inclination_deg(altitude_km)
        if not math.isfinite(inclination):
            return {"feasible": False, "cloud_max_gap_p95_days": float("inf"), "reason": "SSO inclination infeasible"}
        eccentricity = -J3/(2*J2)*EARTH_RADIUS_M/(EARTH_RADIUS_M + altitude_km*1000)*math.sin(math.radians(inclination))
        mission = replace(self.base, altitude_km=float(altitude_km), inclination_deg=inclination,
                          eccentricity=eccentricity, satellites=int(satellite_count),
                          acquired_swath_km=float(acquired_swath_km))
        mission.validate()
        result_key = hashlib.sha256(json.dumps({"mission": asdict(mission), "clouds": self.cloud_parameters,
                                                "backend": self.backend, "criterion": self.criterion, "probe": probe}, sort_keys=True).encode()).hexdigest()
        if result_key in self._results:
            return self._results[result_key]
        frame, sun, provenance = self._orbit(mission)
        subset = None
        if probe:
            grid = target_grid(mission)
            latitudes = sorted(grid.latitude_deg.unique())
            longitudes = sorted(grid.longitude_deg.unique())
            subset = grid[grid.latitude_deg.isin([latitudes[0], latitudes[len(latitudes)//2], latitudes[-1]]) &
                          grid.longitude_deg.isin([longitudes[0], longitudes[len(longitudes)//2], longitudes[-1]])]
        grid_complete = not probe
        if self.criterion == 'matlab-mean' and not probe:
            grid=target_grid(mission)
            center_lon=(mission.longitude_min_deg+mission.longitude_max_deg)/2
            center_lat=(mission.latitude_min_deg+mission.latitude_max_deg)/2
            ordered=grid.assign(_lon=abs(grid.longitude_deg-center_lon),_lat=abs(grid.latitude_deg-center_lat)).sort_values(['_lon','_lat']).drop(columns=['_lon','_lat'])
            acquired_parts,metric_parts=[],[]
            for start in range(0,len(ordered),16):
                acquired,_,metrics=analyze_accesses(frame,sun,mission,targets=ordered.iloc[start:start+16],include_contacts=False)
                acquired_parts.append(acquired);metric_parts.append(metrics)
                fail=(metrics.access_count<2).any() or (metrics.observed_max_gap_days/(1-self.cloud_parameters['cloud_fraction'])>self.cloud_parameters['requirement_days']).any()
                if fail:
                    grid_complete=start+16>=len(ordered)
                    break
            acquisitions=pd.concat(acquired_parts,ignore_index=True)
            targets=pd.concat(metric_parts,ignore_index=True)
            contacts=pd.DataFrame()
        else:
            acquisitions, contacts, targets = analyze_accesses(frame, sun, mission, targets=subset,
                                                               include_contacts=self.criterion!='matlab-mean')
        clouds = None
        p95 = None
        geometric = float(targets.max_gap_days.max())
        reobserved = bool((targets.access_count >= 2).all())
        observed = float(targets.observed_max_gap_days.max()) if reobserved else 1e6
        mean_corrected = observed / (1 - self.cloud_parameters['cloud_fraction'])
        if self.criterion == 'p95':
            clouds, _ = apply_markov_clouds(acquisitions, targets, mission.duration_days*86400,
                                           epoch=mission.epoch, **self.cloud_parameters)
            p95 = float(clouds["regional_max_gap_days"]["p95"])
            selection_days = p95
        else:
            selection_days = mean_corrected
        result = {"feasible": selection_days <= self.cloud_parameters["requirement_days"],
                  "selection_criterion": self.criterion, "selection_revisit_days": selection_days,
                  "geometric_max_gap_days": geometric,
                  "geometric_observed_max_gap_days": observed,
                  "mean_corrected_revisit_days": mean_corrected,
                  "all_targets_reobserved": reobserved, "probe_only": probe, "grid_complete": grid_complete,
                  "cloud_max_gap_p95_days": p95,
                  "target_count": len(targets), "acquisition_count": len(acquisitions),
                  "unobserved_targets": int(targets.never_observed.sum()),
                  "contact_count": len(contacts) if self.criterion != "matlab-mean" else None, "clouds": clouds,
                  "propagation": provenance, "mission": asdict(mission)}
        self._results[result_key] = result
        return result
