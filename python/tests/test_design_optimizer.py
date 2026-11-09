"""Parity of the independent design search with the authoritative MATLAB inputs."""
from dataclasses import replace
import json
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from eo_mission.design_optimizer import enumerate_candidates, optical_grids, select_design
from eo_mission.thesis_inputs import _scalar, load_thesis_inputs

ROOT = Path(__file__).resolve().parents[2]


@pytest.fixture(scope='module')
def computed_optics():
    inputs = load_thesis_inputs()
    return inputs, optical_grids(inputs)


@pytest.mark.parametrize('quantity', ['SNR', 'MTF'])
@pytest.mark.parametrize('band,detector,telescope', [(b,d,t) for b in range(3) for d in range(3) for t in range(4)])
def test_all_optical_grids_match_independent_matlab_exports(computed_optics, quantity, band, detector, telescope):
    inputs, (mtf, snr) = computed_optics
    export_detector = detector + 1 + (3 if band == 2 else 0)
    file = ROOT/'Compiled data/Pipeline'/quantity/f'{quantity}_Lambda{band+1}_Detector{export_detector}_Telescopio{telescope+1}_resultados.csv'
    reference = pd.read_csv(file, index_col=0)
    np.testing.assert_array_equal(reference.index.str.replace('Alt_','').astype(float), inputs.altitudes_km)
    np.testing.assert_array_equal(reference.columns.str.replace('Diam_','').astype(float), inputs.diameters_mm)
    actual = (mtf if quantity == 'MTF' else snr)[band,detector,telescope]
    np.testing.assert_allclose(actual, reference.to_numpy(), rtol=2e-12, atol=1e-12, equal_nan=True)


def test_inputs_are_read_from_master_not_the_selected_result(tmp_path):
    master = ROOT/'MATLAB_Code/master.m'
    altered = tmp_path/'master.m'
    altered.write_text(master.read_text().replace('config.groundSampleDistance_m = 80;', 'config.groundSampleDistance_m = 90;'))
    inputs = load_thesis_inputs(altered)
    assert inputs.gsd_m == 90
    assert len(inputs.altitudes_km) == 81 and len(inputs.diameters_mm) == 181
    assert inputs.cloud_fraction == 1/6
    assert inputs.constellations == ((1,1),(2,1),(2,2),(1,2),(3,1))
    assert inputs.co2['names'][1] == 'H2RG'
    assert inputs.co2['pixels'][1] == 2048
    assert inputs.fov_deg == (10,3,3,8)


def test_numeric_input_reader_never_executes_source():
    with pytest.raises(ValueError):
        _scalar("__import__('os').system('false')")


def test_candidate_calculation_does_not_read_any_saved_solution(monkeypatch):
    def forbidden(*args, **kwargs): raise AssertionError('A saved solution was read')
    monkeypatch.setattr(pd, 'read_csv', forbidden)
    inputs = replace(load_thesis_inputs(), altitudes_km=(690,), diameters_mm=(37,38,39), constellations=((2,2),))
    candidates = enumerate_candidates(inputs)
    h2rg = candidates[(candidates.ID_Detector == 2)&(candidates.Tipo_Telescopio == 'Refractivo')].iloc[0]
    assert h2rg.Diametro_Pupila_mm == 38
    assert h2rg.Masa_Seca_Satelite_kg == pytest.approx(3.3417048)
    assert h2rg.Masa_Total_Constelacion_kg == pytest.approx(6.78401841709255)
    assert h2rg.MTF_a_Nyquist_lambda2 == pytest.approx(0.252525753867215)
    assert h2rg.SNR_lambda3 == pytest.approx(737.998946206156)


def test_mass_ordered_selection_rejects_failed_coverage_and_uses_first_feasible_swath():
    rows = pd.DataFrame([
        {'Altura_km':500,'Num_Satelites':1,'Masa_Total_Constelacion_kg':1,'admissible_swaths':[(20,20,1,25)]},
        {'Altura_km':600,'Num_Satelites':2,'Masa_Total_Constelacion_kg':2,'admissible_swaths':[(20,20,1,25),(30,30,2,35)]},
        {'Altura_km':700,'Num_Satelites':3,'Masa_Total_Constelacion_kg':3,'admissible_swaths':[(40,40,2,45)]},
    ])
    class Coverage:
        def evaluate(self,h,n,w):
            gap = 6 if w>=30 else 9
            return dict(feasible=gap<=7,cloud_max_gap_p95_days=gap,geometric_max_gap_days=gap-1)
    best, calls = select_design(rows,Coverage(),7)
    assert best['Altura_km']==600 and best['Swath_Requerido_km']==30
    assert len(calls)==3


def test_failed_run_does_not_leave_a_stale_feasible_result(tmp_path,monkeypatch):
    import eo_mission.design_optimizer as optimizer
    target=tmp_path/'design_summary.json'
    target.write_text('{"status":"feasible","best_solution":{"Altura_km":999}}')
    def fail(*args,**kwargs): raise RuntimeError('Test propagation failed')
    monkeypatch.setattr(optimizer,'_run_design',fail)
    with pytest.raises(RuntimeError,match='Test propagation failed'):
        optimizer.run_design(tmp_path)
    state=json.loads(target.read_text())
    assert state['status']=='failed' and state['error']=='Test propagation failed'
    assert 'best_solution' not in state
