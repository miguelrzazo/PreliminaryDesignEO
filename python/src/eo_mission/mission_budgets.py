"""Thesis end-of-life, raw-data and memory budgets for a computed design."""
from __future__ import annotations
import math
import numpy as np
import pandas as pd
from scipy.optimize import brentq

from .mass import _decay_cycle_time, _hohmann_dv, _R_EARTH, _MU_EARTH, _ISP, _G0

DAY = 86400.0


def end_of_life(altitude_km: float, dry_mass_kg: float, area_m2: float,
                maintenance_fuel_kg: float, satellites: int) -> dict:
    """Integrate the same autonomous circular-equivalent decay law as MATLAB.

    Piecewise quadrature and root inversion solve the scalar altitude ODE,
    retaining its mass changes at manoeuvres and original 365-day deorbit target.
    """
    high, low = 1.02*altitude_km*1000, .98*altitude_km*1000
    area = 1.2*area_m2
    mass = dry_mass_kg+maintenance_fuel_kg
    remaining = 8*365.25*DAY
    count = 0
    while True:
        cycle = _decay_cycle_time(low, high, mass, area)
        if remaining <= cycle:
            end_height = brentq(lambda h: _decay_cycle_time(h, high, mass, area)-remaining,
                               low, np.nextafter(high, low), xtol=1e-5)
            break
        remaining -= cycle
        count += 1
        mass /= math.exp(_hohmann_dv(_R_EARTH+low, _R_EARTH+high)/(_ISP*_G0))
    natural = _decay_cycle_time(100000., end_height, mass, area)/DAY
    def decay_days(perigee):
        return _decay_cycle_time(100000., (end_height+perigee)/2, dry_mass_kg, area)/DAY
    if decay_days(100000.) > 365:
        return dict(feasible=False, reason='The circular-equivalent model cannot meet 365 days within its perigee bounds',
                    end_altitude_km=end_height/1000, natural_decay_days=natural)
    perigee = end_height if decay_days(end_height) <= 365 else brentq(lambda hp: decay_days(hp)-365,100000.,end_height,xtol=1e-5)
    apogee_radius = _R_EARTH+end_height
    transfer_axis = _R_EARTH+(end_height+perigee)/2
    dv = math.sqrt(_MU_EARTH/apogee_radius)-math.sqrt(_MU_EARTH*(2/apogee_radius-1/transfer_axis))
    fuel = dry_mass_kg*math.expm1(dv/(_ISP*_G0))*1.1
    final = dry_mass_kg+maintenance_fuel_kg+fuel
    return dict(feasible=True, end_altitude_km=end_height/1000, simulated_maintenance_manoeuvres=count,
                natural_decay_days=natural, deorbit_perigee_km=perigee/1000, deorbit_delta_v_m_s=dv,
                maintenance_fuel_kg=maintenance_fuel_kg, deorbit_fuel_kg=fuel,
                total_satellite_mass_kg=final, total_constellation_mass_kg=final*satellites,
                active_reentry_days=decay_days(perigee), model='MATLAB circular-equivalent drag law; piecewise quadrature/inversion')


def clipped_ground_segments(latitude,longitude,bounds=(25,49,-125,-66)) -> np.ndarray:
    """MATLAB rectangular segment clipping and 6371-km haversine distances."""
    lat,lon=np.asarray(latitude),np.asarray(longitude)
    dy,dx=np.diff(lat),np.diff(lon)
    start,end=np.zeros(len(dx)),np.ones(len(dx))
    valid=(np.abs(dx)<=180) & np.isfinite(lat[:-1]) & np.isfinite(lat[1:]) & np.isfinite(lon[:-1]) & np.isfinite(lon[1:])
    for x,d,lower,upper in ((lon[:-1],dx,bounds[2],bounds[3]),(lat[:-1],dy,bounds[0],bounds[1])):
        moving=d!=0
        valid &= moving | ((x>=lower)&(x<=upper))
        left=np.divide(lower-x,d,out=np.full_like(d,-np.inf,dtype=float),where=moving)
        right=np.divide(upper-x,d,out=np.full_like(d,np.inf,dtype=float),where=moving)
        start=np.maximum(start,np.minimum(left,right))
        end=np.minimum(end,np.maximum(left,right))
    valid &= start<=end
    lat_a,lat_b=lat[:-1]+start*dy,lat[:-1]+end*dy
    lon_a,lon_b=lon[:-1]+start*dx,lon[:-1]+end*dx
    a=np.sin(np.radians(lat_b-lat_a)/2)**2+np.cos(np.radians(lat_a))*np.cos(np.radians(lat_b))*np.sin(np.radians(lon_b-lon_a)/2)**2
    lengths=2*6371*np.arctan2(np.sqrt(np.maximum(0,a)),np.sqrt(np.maximum(0,1-a)))
    return np.where(valid,lengths,0.)


def contact_and_memory(ephemeris: pd.DataFrame,contacts: pd.DataFrame,*,swath_km:float,
                       gsd_m:float,telescopes:int,rate_mbps:float=10,memory_gb:float=2) -> dict:
    """Keep MATLAB raw acquisition and sequential unbounded-backlog equations."""
    active=math.ceil(swath_km*1000/(gsd_m*telescopes))
    rows=[]
    for satellite,frame in ephemeris.groupby('satellite',sort=True):
        frame=frame.sort_values('time_s')
        times=frame.time_s.to_numpy()
        lengths=clipped_ground_segments(frame.latitude_deg.to_numpy(),frame.longitude_deg.to_numpy())
        generated=lengths*1000/gsd_m*(telescopes*active)*(3*12)/8e9
        events=contacts[contacts.satellite==satellite]
        contact_s=float(events.duration_s.sum())
        mask=np.zeros(len(times),dtype=bool)
        for event in events.itertuples(): mask |= (times>=event.start_s)&(times<=event.end_s)
        required=peak=0.
        for i,amount in enumerate(generated,1):
            required=max(0.,required+amount-mask[i]*rate_mbps*1e6*(times[i]-times[i-1])/8e9)
            peak=max(peak,required)
        data=float(generated.sum())
        rows.append(dict(satellite=int(satellite),trace_km=float(lengths.sum()),data_gb=data,
                         contact_s=contact_s,minimum_rate_mbps=data*8000/contact_s if contact_s else None,
                         peak_memory_gb=peak))
    rates=[r['minimum_rate_mbps'] for r in rows]
    return dict(per_satellite=rows,active_pixels_per_telescope=active,bits_per_pixel=36,
                acquired_swath_km=telescopes*active*gsd_m/1000,
                required_rate_mbps=max(rates) if all(x is not None for x in rates) else None,
                reference_rate_mbps=rate_mbps,memory_capacity_gb=memory_gb,
                link_pass=all(x is not None and x<=rate_mbps for x in rates),
                memory_pass=all(r['peak_memory_gb']<=memory_gb for r in rows),
                model='MATLAB raw rectangular-region acquisition; three 12-bit bands; cloud losses do not reduce conservative raw data sizing')
