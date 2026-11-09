"""Check candidate independence and orbit-cache reuse without a JVM."""
from dataclasses import replace

import numpy as np
import pandas as pd
import pytest

from eo_mission.coverage import sso_inclination_deg
from eo_mission.coverage_evaluator import MissionCoverageEvaluator
from eo_mission.mission import EARTH_RADIUS_M, SelectedMission


def test_candidate_inputs_replace_orbit_and_swath_without_an_adopted_csv(tmp_path, monkeypatch):
    from eo_mission import mission_reference
    calls = []
    def propagate(mission):
        calls.append(mission)
        xyz = EARTH_RADIUS_M + mission.altitude_km*1000
        frame = pd.DataFrame({"time_s": [0., 100.], "satellite": [1, 1],
                              "x_m": [xyz, xyz], "y_m": [0., 0.], "z_m": [0., 0.],
                              "latitude_deg": [0., 0.], "longitude_deg": [0., 0.],
                              "altitude_km": [mission.altitude_km]*2, "sunlit": [True, True]})
        return frame, np.tile([1.5e11, 0., 0.], (2, 1)), {"backend": "test fixture"}
    monkeypatch.setattr(mission_reference, "propagate_reference", propagate)
    base = replace(SelectedMission(), satellites=1, duration_days=100/86400,
                   latitude_min_deg=-0.001, latitude_max_deg=0.001,
                   longitude_min_deg=-0.001, longitude_max_deg=0.001,
                   station_latitude_deg=0, station_longitude_deg=0)
    kwargs = dict(base_mission=base, cache_dir=tmp_path, backend="scipy-reference",
                  cloud_fraction=0, mc_samples=3)
    evaluator = MissionCoverageEvaluator(**kwargs)
    first = evaluator.evaluate(700, 1, 100)
    assert calls[0].inclination_deg == pytest.approx(sso_inclination_deg(700))
    assert calls[0].altitude_km == 700
    assert first["mission"]["acquired_swath_km"] == 100
    assert first["target_count"] == 4
    assert first["feasible"]
    wider = evaluator.evaluate(700, 1, 200)
    assert len(calls) == 1
    assert wider["cloud_max_gap_p95_days"] <= first["cloud_max_gap_p95_days"]
    fresh = MissionCoverageEvaluator(**kwargs)
    fresh.evaluate(700, 1, 150)
    assert len(calls) == 1
    fresh.evaluate(710, 1, 150)
    assert len(calls) == 2


def test_mean_cloud_criterion_uses_max_geometric_gap_and_not_p95(tmp_path, monkeypatch):
    import eo_mission.coverage_evaluator as module
    def analyzed(*args, **kwargs):
        metrics=pd.DataFrame(dict(target=['A','B'],max_gap_days=[6.5,6.0],
            observed_max_gap_days=[5.5,5.0],access_count=[3,4],never_observed=[False,False]))
        return pd.DataFrame(),pd.DataFrame(),metrics
    monkeypatch.setattr(module,'analyze_accesses',analyzed)
    def forbidden(*args, **kwargs): raise AssertionError('P95 clouds changed mean selection')
    monkeypatch.setattr(module,'apply_markov_clouds',forbidden)
    evaluator=MissionCoverageEvaluator(cache_dir=tmp_path,cloud_fraction=1/6,criterion='matlab-mean')
    monkeypatch.setattr(evaluator,'_orbit',lambda mission:(pd.DataFrame(),np.array([]),{}))
    result=evaluator.evaluate(700,2,240)
    assert result['selection_revisit_days']==pytest.approx(5.5/(1-1/6))
    assert result['feasible'] and result['cloud_max_gap_p95_days'] is None
    assert result['geometric_max_gap_days']==6.5


def test_p95_is_an_explicit_different_selection_objective(tmp_path, monkeypatch):
    import eo_mission.coverage_evaluator as module
    metrics=pd.DataFrame(dict(target=['A'],max_gap_days=[5.0],observed_max_gap_days=[5.0],access_count=[4],never_observed=[False]))
    monkeypatch.setattr(module,'analyze_accesses',lambda *a,**kw:(pd.DataFrame(),pd.DataFrame(),metrics))
    monkeypatch.setattr(module,'apply_markov_clouds',lambda *a,**kw:({'regional_max_gap_days':{'p95':11.}},pd.DataFrame()))
    evaluator=MissionCoverageEvaluator(cache_dir=tmp_path,cloud_fraction=1/6,criterion='p95')
    monkeypatch.setattr(evaluator,'_orbit',lambda mission:(pd.DataFrame(),np.array([]),{}))
    result=evaluator.evaluate(700,2,240)
    assert result['mean_corrected_revisit_days']==6
    assert result['selection_revisit_days']==11 and not result['feasible']


@pytest.mark.parametrize('gap,expected_complete',[(5.,True),(8.,False)])
def test_mean_full_grid_can_reject_early_but_cannot_accept_a_partial_grid(tmp_path,monkeypatch,gap,expected_complete):
    import eo_mission.coverage_evaluator as module
    calls=[]
    def analyze(frame,sun,mission,*,targets=None,include_contacts=True):
        calls.append(list(targets.target))
        metrics=targets.copy()
        metrics['max_gap_days']=gap+.5
        metrics['observed_max_gap_days']=gap
        metrics['access_count']=3
        metrics['never_observed']=False
        return pd.DataFrame(columns=['target','start_s','end_s']),pd.DataFrame(),metrics
    monkeypatch.setattr(module,'analyze_accesses',analyze)
    evaluator=MissionCoverageEvaluator(cache_dir=tmp_path,cloud_fraction=1/6,criterion='matlab-mean')
    monkeypatch.setattr(evaluator,'_orbit',lambda mission:(pd.DataFrame(),np.array([]),{}))
    result=evaluator.evaluate(700,2,240)
    assert result['grid_complete']==expected_complete
    assert result['feasible']==expected_complete
    if expected_complete:
        assert result['target_count']==273
        assert len({target for batch in calls for target in batch})==273
    else:
        assert len(calls)==1 and result['target_count']==16
