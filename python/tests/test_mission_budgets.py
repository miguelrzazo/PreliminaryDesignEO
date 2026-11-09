"""Engineering budgets compared with MATLAB and simple conservation cases."""
from pathlib import Path
import numpy as np
import pandas as pd
import pytest

from eo_mission.mass import compute_dry_mass, compute_total_mass
from eo_mission.mission_budgets import clipped_ground_segments, contact_and_memory, end_of_life


def test_reentry_budget_matches_matlab_with_explicit_solver_tolerance():
    dry=compute_dry_mass([38],2).iloc[0]
    total=compute_total_mass([690],dry.Masa_seca,dry.satelliteArea_m2).iloc[0]
    result=end_of_life(690,dry.Masa_seca,dry.satelliteArea_m2,total.Masa_combustible,2)
    reference=pd.read_csv(Path(__file__).parents[2]/'Compiled data/Mission/EndOfLife/end_of_life.csv').iloc[0]
    assert result['feasible']
    assert result['simulated_maintenance_manoeuvres']==reference.SimulatedMaintenanceManoeuvres
    assert result['end_altitude_km']==pytest.approx(reference.EndAltitude_km,abs=.002)
    # MATLAB ode45 crosses atmosphere-layer discontinuities; piecewise quadrature
    # yields the same decay model with small numerical integration differences.
    assert result['natural_decay_days']==pytest.approx(reference.NaturalDecay_days,abs=.1)
    assert result['deorbit_perigee_km']==pytest.approx(reference.DeorbitPerigee_km,abs=.01)
    assert result['deorbit_delta_v_m_s']==pytest.approx(reference.DeorbitDeltaV_ms,abs=.002)
    assert result['total_satellite_mass_kg']==pytest.approx(reference.TotalSatelliteMass_kg,abs=3e-6)
    assert result['active_reentry_days']==pytest.approx(365,abs=1e-7)


def test_rectangular_region_clips_a_crossing_and_rejects_dateline_wrap():
    distance=clipped_ground_segments([24,26],[-100,-100])[0]
    assert distance==pytest.approx(6371*np.radians(1))
    assert clipped_ground_segments([40,40],[179,-179])[0]==0
    assert clipped_ground_segments([0,0],[0,1])[0]==0
    assert clipped_ground_segments([40,40],[-100,-100])[0]==0


def test_raw_data_and_no_contact_backlog_conserve_generated_bits():
    states=pd.DataFrame(dict(satellite=[1,1,1],time_s=[0.,10.,20.],latitude_deg=[40.,40.,40.],longitude_deg=[-100.,-99.99,-99.98]))
    contacts=pd.DataFrame(columns=['satellite','start_s','end_s','duration_s'])
    result=contact_and_memory(states,contacts,swath_km=240,gsd_m=80,telescopes=2)
    row=result['per_satellite'][0]
    assert result['active_pixels_per_telescope']==1500
    assert row['data_gb']==pytest.approx(row['trace_km']*1000/80*3000*36/8e9)
    assert row['peak_memory_gb']==pytest.approx(row['data_gb'])
    assert row['minimum_rate_mbps'] is None and not result['link_pass']


def test_download_uses_time_step_and_keeps_unbounded_backlog():
    states=pd.DataFrame(dict(satellite=[1,1,1],time_s=[0.,10.,20.],latitude_deg=[40.,40.,40.],longitude_deg=[-100.,-99.99,-99.98]))
    contacts=pd.DataFrame(dict(satellite=[1],start_s=[0.],end_s=[20.],duration_s=[20.]))
    result=contact_and_memory(states,contacts,swath_km=240,gsd_m=80,telescopes=2,rate_mbps=.01,memory_gb=0)
    row=result['per_satellite'][0]
    assert row['peak_memory_gb']==pytest.approx(row['data_gb']-20*.01*1e6/8e9)
    assert not result['memory_pass']
