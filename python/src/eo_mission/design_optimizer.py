"""Recompute the thesis trade space and select by maintenance-inclusive mass.

Inputs come from master.m, never selected_design.csv. Optical/spacecraft laws
and first-feasible grid policies follow MATLAB; timed Orekit accesses and a
specified Markov percentile replace its analytical mean-cloud revisit model.
"""
from __future__ import annotations

from dataclasses import asdict
import hashlib
import json
import math
from pathlib import Path

import numpy as np
import pandas as pd

from .mass import compute_dry_mass, compute_total_mass
from .mission import SelectedMission
from .mission_clouds import cloud_hardware_floor
from .optics import compute_mtf, compute_snr
from .thesis_inputs import ThesisInputs, load_thesis_inputs


def optical_grids(inputs: ThesisInputs) -> tuple[dict, dict]:
    """Evaluate all three spectral bands with the matching MATLAB catalog rows."""
    mtf, snr = {}, {}
    for band, wavelength in enumerate(inputs.wavelengths_m):
        detectors = inputs.co2 if band < 2 else inputs.o2
        for detector in range(3):
            for telescope in range(len(inputs.telescope_names)):
                common = dict(pixel_size=detectors['pitch_m'][detector],
                              gsd=inputs.gsd_m, r_obs=inputs.obstruction[telescope],
                              altitudes=inputs.altitudes_km, diameters=inputs.diameters_mm)
                mtf[band, detector, telescope] = compute_mtf(
                    lambda_=wavelength, mtf_detector=detectors['mtf'][detector],
                    mtf_alignment=inputs.alignment[telescope], mtf_req=inputs.mtf_min, **common).to_numpy()
                snr[band, detector, telescope] = compute_snr(
                    lambda_c=wavelength, eta=detectors['eta'][detector],
                    tau=inputs.transmission[telescope], snr_req=inputs.snr_min, **common).to_numpy()
    return mtf, snr


def enumerate_candidates(inputs: ThesisInputs) -> pd.DataFrame:
    """Enumerate hardware using MATLAB's binding-band/minimum-pupil policy.

    MATLAB selects MTF at lambda2 and SNR at lambda3 using corresponding detector
    rows, then ranks by dry plus maintenance propellant. Deorbit is a later
    budget and is not silently added to the optimization objective.
    """
    mtf, snr = optical_grids(inputs)
    rows, mass_cache = [], {}
    diameters = np.asarray(inputs.diameters_mm)
    for config_index, (satellites, telescopes) in enumerate(inputs.constellations):
        dry = compute_dry_mass(diameters, telescopes)
        for detector in range(3):
            for telescope, name in enumerate(inputs.telescope_names):
                mask = (mtf[1, detector, telescope] >= inputs.mtf_min) & (snr[2, detector, telescope] >= inputs.snr_min)
                for h_index, altitude in enumerate(inputs.altitudes_km):
                    feasible_pupils = np.flatnonzero(mask[h_index])
                    if not len(feasible_pupils): continue
                    d_index = int(feasible_pupils[0])
                    diameter = float(diameters[d_index])
                    if telescope == 0 and diameter > 200: continue
                    max_fov = 2 * altitude * math.tan(math.radians(telescopes * inputs.fov_deg[telescope] / 2))
                    detector_width = inputs.gsd_m * inputs.co2['pixels'][detector] / 1000
                    swaths = []
                    for width in inputs.swaths_km:
                        if width > max_fov or width > inputs.max_detectors * detector_width: continue
                        detector_count = max(telescopes, math.ceil(width / detector_width))
                        acquired = telescopes * math.ceil(width * 1000 / (inputs.gsd_m * telescopes)) * inputs.gsd_m / 1000
                        actual_width = min(detector_count * detector_width, max_fov)
                        if acquired <= actual_width:
                            swaths.append((float(width), float(acquired), detector_count, actual_width))
                    if not swaths: continue
                    properties = dry.iloc[d_index]
                    key = (altitude, diameter, telescopes)
                    if key not in mass_cache:
                        result = compute_total_mass(
                            [altitude], [properties['Masa_seca']], [properties['satelliteArea_m2']],
                            inputs.maintenance_lower, inputs.maintenance_upper).iloc[0]
                        mass_cache[key] = result
                    mass = mass_cache[key]
                    if not np.isfinite(mass['Masa_total']): continue
                    n = math.sqrt(3.986004418e14 / (6378137 + altitude * 1000)**3)
                    inclination = math.degrees(math.acos(-(2 * math.pi / (365.2421897*86400)) /
                        (1.5 * n * 1.082629989052e-3 * (6378137 / (6378137+altitude*1000))**2)))
                    rows.append({
                        'configuration_index': config_index, 'Num_Satelites': satellites,
                        'Num_Telescopios': telescopes, 'ID_Detector': detector+1,
                        'Detector_CO2': inputs.co2['names'][detector], 'Detector_O2': inputs.o2['names'][detector],
                        'Tipo_Telescopio': name, 'Altura_km': altitude, 'Diametro_Pupila_mm': diameter,
                        'Inclinacion_SSO_deg': inclination, 'LTAN_h': inputs.ltan_h, 'GSD': inputs.gsd_m,
                        'MTF_a_Nyquist_lambda2': float(mtf[1,detector,telescope][h_index,d_index]),
                        'SNR_lambda3': float(snr[2,detector,telescope][h_index,d_index]),
                        'Masa_Seca_Satelite_kg': float(properties['Masa_seca']),
                        'Area_Satelite_m2': float(properties['satelliteArea_m2']),
                        'Masa_Combustible_kg': float(mass['Masa_combustible']),
                        'Masa_Total_Satelite_kg': float(mass['Masa_total']),
                        'Masa_Total_Constelacion_kg': float(mass['Masa_total']*satellites),
                        'Num_Impulsos_Mantenimiento': int(mass['Num_impulsos']),
                        'DeltaV_Total_ms': float(mass['DeltaV_total']),
                        'admissible_swaths': swaths,
                    })
    if not rows: return pd.DataFrame()
    # Stable ties retain MATLAB's configuration/detector/telescope/altitude order.
    return pd.DataFrame(rows).sort_values('Masa_Total_Constelacion_kg', kind='stable').reset_index(drop=True)


def select_design(candidates: pd.DataFrame, evaluator, requirement_days: float) -> tuple[dict | None, list[dict]]:
    """Evaluate in mass order; the first feasible hardware row minimizes that objective."""
    evaluated, checked_cache = [], {}
    def checked(row, swath, probe=False):
        method = evaluator.evaluate_probe if probe and hasattr(evaluator, 'evaluate_probe') else evaluator.evaluate
        key=(row['Altura_km'],row['Num_Satelites'],swath[1],probe and hasattr(evaluator,'evaluate_probe'))
        if key in checked_cache: return checked_cache[key]
        result = method(row['Altura_km'], row['Num_Satelites'], swath[1])
        days = result.get('selection_revisit_days', result.get('cloud_max_gap_p95_days'))
        evaluated.append({'altitude_km': row['Altura_km'], 'satellites': row['Num_Satelites'],
                          'acquired_swath_km': swath[1], 'selection_days': days,
                          'geometric_days': result['geometric_max_gap_days'],
                          'probe_only': result.get('probe_only', False), 'feasible': bool(result['feasible'])})
        checked_cache[key]=(result,days)
        return result, days
    for index, row in enumerate(candidates.to_dict('records')):
        swaths = row['admissible_swaths']
        # With identical orbits/cloud chains, narrower footprints cannot improve
        # an access-gap maximum. Test the widest first to prune entire hardware rows.
        result, _ = checked(row, swaths[-1], probe=True)
        if not result['feasible']:
            if index % 50 == 0: print(f"DESIGN_PROGRESS hardware={index+1} rejected_by_probe", flush=True)
            continue
        result, _ = checked(row, swaths[-1])
        if not result['feasible']: continue
        lo, hi = 0, len(swaths)-1
        while lo < hi:
            mid = (lo+hi)//2
            trial, _ = checked(row, swaths[mid])
            if trial['feasible']: hi = mid
            else: lo = mid+1
        width, acquired, detector_count, actual_width = swaths[lo]
        result, days = checked(row, swaths[lo])
        if result['feasible'] and days <= requirement_days:
            solution = {k:v for k,v in row.items() if k != 'admissible_swaths'}
            solution.update(Swath_Requerido_km=width, Swath_Adquirido_km=acquired,
                            Detectores_Necesarios=detector_count, Swath_Real_Sistema_km=actual_width,
                            Dias_Cobertura=days, coverage=result)
            return solution, evaluated
    return None, evaluated


def _run_design(output_dir: str | Path, *, master_path: str | Path | None = None,
               cloud_fraction: float | None = None, persistence: float = .5,
               mc_samples: int = 200, seed: int = 12345, duration_days: float = 30,
               data_path: str | Path | None = None, criterion: str = "matlab-mean") -> Path:
    """Compute a fresh design or certify infeasibility for the specified sampled objective."""
    inputs = load_thesis_inputs(master_path)
    p = inputs.cloud_fraction if cloud_fraction is None else cloud_fraction
    if criterion not in {'matlab-mean','p95'}:
        raise ValueError('Unknown design criterion')
    if not 0 <= p < 1:
        raise ValueError('Design cloud fraction must be in [0,1)')
    if duration_days <= inputs.revisit_days:
        raise ValueError('Use a design horizon longer than the revisit requirement to avoid censored-horizon passes')
    base = SelectedMission(duration_days=duration_days, step_s=30, ltan_h=inputs.ltan_h,
                           overlap_fraction=inputs.overlap, gsd_m=inputs.gsd_m)
    floor = cloud_hardware_floor(duration_days*86400, cloud_fraction=p, persistence=persistence,
        mc_samples=mc_samples, seed=seed, epoch=base.epoch, requirement_days=inputs.revisit_days)
    candidates = enumerate_candidates(inputs)
    best, evaluated = None, []
    prune = criterion == 'p95' and floor['empirical_p95_floor_exceeds_requirement']
    if not prune and not candidates.empty:
        from .coverage_evaluator import MissionCoverageEvaluator
        evaluator = MissionCoverageEvaluator(data_path=data_path, base_mission=base, cloud_fraction=p,
            persistence=persistence, mc_samples=mc_samples, seed=seed, requirement_days=inputs.revisit_days, criterion=criterion)
        best, evaluated = select_design(candidates, evaluator, inputs.revisit_days)
    summary = {
        'status': 'feasible' if best else 'infeasible_for_configured_empirical_objective',
        'best_solution': best, 'hardware_candidates': len(candidates), 'coverage_evaluations': len(evaluated),
        'source': {'inputs': inputs.source_path,
                   'master_sha256': hashlib.sha256(Path(inputs.source_path).read_bytes()).hexdigest(),
                   'selected_design_csv_used': False},
        'inputs': asdict(inputs), 'selection_criterion': criterion, 'cloud_hardware_floor': floor,
        'objective': 'Minimum constellation dry + eight-year maintenance propellant mass; minimum optical pupil and first feasible swath follow MATLAB.',
        'coverage_change': 'Orekit timed point accesses replace the analytical orbit/revisit solver. matlab-mean uses worst observed geometric gap / (1-cloud_fraction), retaining MATLAB mean-availability correction. p95 uses regional worst-target Markov gap. Markov extremes are diagnostic in matlab-mean mode.',
        'reason': ('All hardware is pruned by the regional cloud-spell floor for this exact Monte Carlo sample.'
                   if prune else 'Mass-ordered hardware evaluation completed.'),
        'limits': 'The stochastic reliability objective differs from a mean cloud correction. This sampled-model result is not a climate/operational impossibility theorem. No orbit propagation is needed when the hardware-independent cloud floor already rejects every design.',
        'evaluated_coverage': evaluated,
    }
    out=Path(output_dir); out.mkdir(parents=True,exist_ok=True)
    saved=candidates.copy()
    if not saved.empty: saved['admissible_swaths']=saved.admissible_swaths.map(json.dumps)
    saved.to_csv(out/'hardware_candidates.csv',index=False)
    best_rows=[] if best is None else [{k:v for k,v in best.items() if k!='coverage'}]
    pd.DataFrame(best_rows, columns=[k for k in (best or {}) if k!='coverage'] or ['Es_Viable']).to_csv(out/'best_design.csv',index=False)
    target=out/'design_summary.json'
    target.write_text(json.dumps(summary,indent=2,allow_nan=False)+'\n',encoding='utf-8')
    return finalize_design(target,data_path=data_path) if best is not None else target


def finalize_design(summary_path: str | Path, *, data_path: str | Path | None = None) -> Path:
    """Verify the computed winner at 10 s and export its own mission/budgets.

    This configuration is constructed from the search result. It never imports
    the thesis's adopted-design CSV or presumes its altitude or hardware.
    """
    import yaml
    from .mission import run_selected_mission
    from .mission_budgets import contact_and_memory, end_of_life
    from .mission_visualization import build_mission_explorer
    path=Path(summary_path)
    result=json.loads(path.read_text())
    best=result['best_solution']
    if best is None: return path
    inputs=load_thesis_inputs(result['source']['inputs'])
    properties=compute_dry_mass([best['Diametro_Pupila_mm']],best['Num_Telescopios']).iloc[0]
    focal=inputs.co2['pitch_m'][best['ID_Detector']-1]*best['Altura_km']*1000/inputs.gsd_m
    best.update(Longitud_Focal_m=focal,Relacion_Focal_f_num=focal/(best['Diametro_Pupila_mm']/1000),
                Desfase_Anomalia_Verdadera_deg=360/best['Num_Satelites'],
                Excedente_Swath_km=best['Swath_Real_Sistema_km']-best['Swath_Requerido_km'],
                Excedente_Swath_porc=100*(best['Swath_Real_Sistema_km']-best['Swath_Requerido_km'])/best['Swath_Requerido_km'])
    for output,source in {'Masa_Instrumento_kg':'instrumentMass_kg','Longitud_Instrumento_m':'instrumentLength_m',
        'Area_Instrumento_m2':'instrumentArea_m2','Volumen_Instrumento_m3':'instrumentVolume_m3',
        'Potencia_Electrica_Instrumento_W':'instrumentPower_W','Masa_Seca_TM_kg':'thematicMapperDryMass_kg',
        'Masa_Seca_SEOSAT_kg':'seosatDryMass_kg','Longitud_Satelite_m':'satelliteLength_m',
        'Volumen_Satelite_m3':'satelliteVolume_m3','Estandar_U':'cubesatUnits'}.items():
        best[output]=float(properties[source])
    out=path.parent/'mission'
    out.mkdir(parents=True,exist_ok=True)
    mission=dict(best['coverage']['mission'])
    mission['step_s']=10
    clouds=result['cloud_hardware_floor']['parameters']
    configuration={'mission':mission,'clouds':{
        'cloud_fraction':clouds['cloud_fraction'], 'persistence':clouds['persistence'],
        'mc_samples':clouds['mc_samples'],'seed':clouds['seed'],
        'requirement_days':clouds['requirement_days']}}
    config_path=out/'computed_design.yaml'
    config_path.write_text(yaml.safe_dump(configuration,sort_keys=False))
    run_selected_mission(config_path,out,data_path=data_path,html=False)
    states=pd.read_csv(out/'ephemeris.csv')
    acquisitions=pd.read_csv(out/'acquisitions.csv')
    contacts=pd.read_csv(out/'contacts.csv')
    targets=pd.read_csv(out/'target_metrics.csv')
    mission_summary=json.loads((out/'mission_summary.json').read_text())
    if result['selection_criterion']=='matlab-mean':
        verified=float(targets.observed_max_gap_days.max())/(1-clouds['cloud_fraction'])
        verified_pass=bool((targets.access_count>=2).all() and verified<=clouds['requirement_days'])
    else:
        verified=mission_summary['clouds']['regional_max_gap_days']['p95']
        verified_pass=verified<=clouds['requirement_days']
    if not verified_pass:
        raise RuntimeError('The search candidate failed 10-second verification; refine the search before accepting a final design.')
    best['Dias_Cobertura']=verified
    best['coverage']['verified_10s_selection_revisit_days']=verified
    # The MATLAB contact and storage budget is a seven-day raw-data study.
    end_s=7*86400
    weekly=states[states.time_s<=end_s]
    weekly_contacts=contacts[contacts.start_s<end_s].copy()
    weekly_contacts['end_s']=weekly_contacts.end_s.clip(upper=end_s)
    weekly_contacts['duration_s']=weekly_contacts.end_s-weekly_contacts.start_s
    budgets=contact_and_memory(weekly,weekly_contacts,swath_km=best['Swath_Adquirido_km'],
        gsd_m=best['GSD'],telescopes=best['Num_Telescopios'])
    lifecycle=end_of_life(best['Altura_km'],best['Masa_Seca_Satelite_kg'],best['Area_Satelite_m2'],
                          best['Masa_Combustible_kg'],best['Num_Satelites'])
    result['mission_budgets']=budgets
    result['end_of_life']=lifecycle
    result['mission_checks_pass']=bool(budgets['link_pass'] and budgets['memory_pass'] and lifecycle['feasible'])
    result['verification']={'step_s':10,'selection_revisit_days':verified,'pass':verified_pass,
                            'source_selected_design_csv_used':False}
    mission_summary['provenance']['result_role']='independently_optimized'
    mission_summary['definitions']['thesis_baseline']='Hardware and orbit were calculated by the Python design search from MATLAB input grids. No adopted-design result was loaded.'
    mission_summary['design']={k:v for k,v in best.items() if k!='coverage'}
    mission_summary['mission_budgets']=budgets
    mission_summary['end_of_life']=lifecycle
    (out/'mission_summary.json').write_text(json.dumps(mission_summary,indent=2,allow_nan=False)+'\n')
    html=build_mission_explorer(states,acquisitions,contacts,mission_summary,out/'mission_explorer.html')
    result['mission_explorer']=str(html.resolve())
    path.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    pd.DataFrame([{k:v for k,v in best.items() if k!='coverage'}]).to_csv(path.parent/'best_design.csv',index=False)
    return path


def run_design(output_dir: str | Path, **kwargs) -> Path:
    """Mark each fresh run explicitly so previous results cannot look current."""
    from datetime import datetime, timezone
    out=Path(output_dir); out.mkdir(parents=True,exist_ok=True)
    target=out/'design_summary.json'
    started=datetime.now(timezone.utc).isoformat()
    target.write_text(json.dumps({'status':'running','started_utc':started,
                                 'selection_criterion':kwargs.get('criterion','matlab-mean')},indent=2)+'\n')
    try:
        result=_run_design(out,**kwargs)
        summary=json.loads(result.read_text())
        summary['started_utc']=started
        summary['completed_utc']=datetime.now(timezone.utc).isoformat()
        result.write_text(json.dumps(summary,indent=2,allow_nan=False)+'\n')
        return result
    except Exception as error:
        summary=json.loads(target.read_text())
        summary.update(status='failed',error=str(error),started_utc=started)
        target.write_text(json.dumps(summary,indent=2,allow_nan=False)+'\n')
        raise
