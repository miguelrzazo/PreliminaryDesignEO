"""A standalone interactive view of propagated mission results."""
from __future__ import annotations

import json
from html import escape
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd


_SATELLITE_COLORS = ("#67b7ff", "#ffc56e", "#b9a0ff", "#68d5b5")
_BACKGROUND = "#101c2c"
_TEXT = "#e8eef5"
_MUTED = "#a9b8c9"
_EPHEMERIS_COLUMNS = {
    "time_s", "satellite", "x_m", "y_m", "z_m", "latitude_deg",
    "longitude_deg", "altitude_km", "sunlit",
}


def _sample_rows(frame: pd.DataFrame, limit: int) -> pd.DataFrame:
    if len(frame) <= limit:
        return frame
    return frame.iloc[np.unique(np.linspace(0, len(frame) - 1, limit).astype(int))]


def _validate_ephemeris(ephemeris: pd.DataFrame) -> pd.DataFrame:
    missing = _EPHEMERIS_COLUMNS.difference(ephemeris.columns)
    if missing:
        raise ValueError("Missing ephemeris columns: " + ", ".join(sorted(missing)))
    if ephemeris.empty:
        raise ValueError("The mission explorer requires propagated ephemeris samples.")
    numeric = sorted(_EPHEMERIS_COLUMNS - {"sunlit"})
    if not np.isfinite(ephemeris[numeric].to_numpy(dtype=float)).all():
        raise ValueError("Ephemeris coordinates and times must be finite.")
    if (ephemeris["altitude_km"] <= 0).any():
        raise ValueError("Ephemeris samples must lie above the Earth.")
    if (ephemeris["latitude_deg"].abs() > 90).any():
        raise ValueError("Ephemeris latitude is outside [-90, 90] degrees.")
    if (ephemeris["longitude_deg"].abs() > 180).any():
        raise ValueError("Ephemeris longitude is outside [-180, 180] degrees.")
    result = ephemeris.sort_values(["satellite", "time_s"]).copy()
    if result.duplicated(["satellite", "time_s"]).any():
        raise ValueError("Each satellite must have one state per sample time.")
    return result


def _earth_position(latitude_deg: Any, longitude_deg: Any) -> tuple[np.ndarray, ...]:
    """WGS84 surface position in kilometres at geodetic latitude/longitude."""
    latitude = np.radians(latitude_deg)
    longitude = np.radians(longitude_deg)
    equatorial_radius = 6378.137
    flattening = 1 / 298.257223563
    eccentricity_squared = flattening * (2 - flattening)
    normal_radius = equatorial_radius / np.sqrt(1 - eccentricity_squared * np.sin(latitude) ** 2)
    return (
        normal_radius * np.cos(latitude) * np.cos(longitude),
        normal_radius * np.cos(latitude) * np.sin(longitude),
        normal_radius * (1 - eccentricity_squared) * np.sin(latitude),
    )


def _target_frame(summary: dict, accesses: pd.DataFrame) -> pd.DataFrame:
    targets = pd.DataFrame(summary.get("targets", summary.get("target_metrics", [])))
    if targets.empty and not accesses.empty:
        columns = ["target", "latitude_deg", "longitude_deg"]
        if set(columns).issubset(accesses.columns):
            targets = accesses[columns].drop_duplicates("target").copy()
            targets["access_count"] = targets["target"].map(accesses.groupby("target").size())
    if not targets.empty:
        required = {"target", "latitude_deg", "longitude_deg"}
        if not required.issubset(targets.columns):
            raise ValueError("Mission targets require target, latitude_deg and longitude_deg.")
        if not np.isfinite(targets[["latitude_deg", "longitude_deg"]].to_numpy(dtype=float)).all():
            raise ValueError("Mission target coordinates must be finite.")
    return targets



def _number(value: Any, digits: int = 0) -> str:
    try:
        return f"{float(value):,.{digits}f}" if np.isfinite(float(value)) else "—"
    except (TypeError, ValueError):
        return "—"


def _validate_events(events: pd.DataFrame) -> None:
    if events.empty:
        return
    required = {"satellite", "target", "start_s", "end_s", "duration_s"}
    if not required.issubset(events.columns):
        raise ValueError("Mission event tables require satellite, target, start_s, end_s and duration_s.")
    numeric = events[["start_s", "end_s", "duration_s"]].to_numpy(dtype=float)
    if not np.isfinite(numeric).all() or (numeric < 0).any():
        raise ValueError("Event times and durations must be finite and nonnegative.")
    if (events["end_s"] < events["start_s"]).any():
        raise ValueError("An event cannot end before it starts.")


def _json_default(value: Any) -> Any:
    if isinstance(value, np.generic):
        return value.item()
    if isinstance(value, (np.ndarray, pd.Series)):
        return value.tolist()
    if isinstance(value, Path):
        return str(value)
    raise TypeError(f"Cannot serialize {type(value).__name__}")


def build_mission_explorer(
    ephemeris: pd.DataFrame,
    accesses: pd.DataFrame,
    contacts: pd.DataFrame,
    summary: dict,
    output_path: str | Path,
) -> Path:
    """Write an offline interactive HTML view of actual mission samples/events.

    Input positions are Earth-fixed metres. Coordinates for targets and the
    station are geodetic degrees. State and event decimation affects display
    only; metric cards retain the counts from the full analysis. The document
    embeds its renderer and data and requires no plotting package or network.
    """
    states = _validate_ephemeris(ephemeris)
    _validate_events(accesses)
    _validate_events(contacts)
    targets = _target_frame(summary, accesses)
    mission = summary.get("mission", {})
    propagation = summary.get("propagation", {})
    metrics = summary.get("metrics", {})
    state_columns = ["time_s", "x_m", "y_m", "z_m", "latitude_deg", "longitude_deg", "altitude_km", "sunlit"]
    groups = [
        {"satellite": int(satellite), "states": _sample_rows(group, 5000)[state_columns].to_numpy().tolist()}
        for satellite, group in states.groupby("satellite", sort=True)
    ]
    duration_s = float(states["time_s"].max())
    if duration_s <= 0:
        raise ValueError("Mission visualization requires a positive time span.")
    sampled_accesses = _sample_rows(accesses.sort_values("start_s"), 10000) if not accesses.empty else accesses
    sampled_contacts = _sample_rows(contacts.sort_values("start_s"), 5000) if not contacts.empty else contacts
    payload = {
        "satellites": groups, "summary": summary, "duration_s": duration_s,
        "coastlines": json.loads(Path(__file__).with_name("coastlines.json").read_text())["lines"],
        "targets": json.loads(targets.to_json(orient="records")),
        "accesses": json.loads(sampled_accesses.to_json(orient="records")),
        "contacts": json.loads(sampled_contacts.to_json(orient="records")),
    }
    payload = json.loads(json.dumps(payload, default=_json_default), parse_constant=lambda _: None)
    data = json.dumps(payload, separators=(",", ":"), allow_nan=False).replace("<", chr(92) + "u003c").replace(">", chr(92) + "u003e").replace("&", chr(92) + "u0026")
    total_events = len(accesses) + len(contacts)
    displayed_events = len(sampled_accesses) + len(sampled_contacts)
    event_note = (f"Showing {displayed_events:,} of {total_events:,} event intervals."
                  if displayed_events < total_events else f"{total_events:,} event intervals across the propagated horizon.")
    cards = "".join(
        f'<div class="metric"><span>{escape(label)}</span><strong>{escape(value)}</strong><small>{escape(unit)}</small></div>'
        for label, value, unit in (
            ("Nominal altitude", _number(mission.get("altitude_km")), "km above reference radius"),
            ("Acquired swath", _number(mission.get("acquired_swath_km")), "km · nadir footprint"),
            ("Ground targets", _number(metrics.get("target_count", len(targets))), "points in the acquisition grid"),
            ("Clear-sky acquisitions", _number(metrics.get("access_count", len(accesses))), "sampled daylight windows"),
            ("Ground contacts", _number(metrics.get("contact_count", len(contacts))), "station elevation windows"),
        )
    )
    badges = [
        f"Epoch {mission.get('epoch', 'Unspecified')}",
        str(propagation.get("model", "See propagation configuration")),
        f"{_number(duration_s / 86400, 1)} days propagated",
        f"{_number(propagation.get('step_s'), 1)} s sampling",
        str(mission.get("station_name", "Configured ground station")),
    ]
    if summary.get("provenance", {}).get("result_role"):
        if summary['provenance']['result_role'] == 'independently_optimized':
            badges.append("Computed design · independent search")
            design=summary.get('design',{})
            badges.append(f"{_number(design.get('Diametro_Pupila_mm'), 0)} mm pupil · {_number(design.get('Num_Telescopios'), 0)} telescopes/satellite")
        else:
            badges.append("Adopted design · validation case")
    metadata = "".join(f"<span>{escape(item)}</span>" for item in badges)
    provenance = summary.get("provenance", {})
    source_note = provenance.get("description", "") if isinstance(provenance, dict) else str(provenance)
    clouds = summary.get("clouds")
    cloud_panel = ""
    if clouds:
        parameters = clouds["parameters"]
        distribution = clouds["regional_max_gap_days"]
        cloud_panel = (
            '<section class="panel"><h2>Persistent cloud scenarios</h2><p class="caption">'
            f'Mean worst-target gap: {_number(distribution["mean"], 2)} days · '
            f'95th percentile: {_number(distribution["p95"], 2)} days. '
            f'Cloud fraction {_number(parameters["cloud_fraction"]*100, 1)}%, persistence '
            f'{_number(parameters["persistence"], 2)}, {_number(parameters["mc_samples"])} realizations. '
            'Daily clouds are shared across the whole study region and applied to the actual acquisition intervals. '
            'Gaps include censored horizon boundaries. This is a finite-horizon Monte Carlo scenario, '
            'not a weather-calibrated or long-term revisit guarantee.</p></section>'
        )
    replacements = {
        "__SATELLITES__": _number(mission.get("satellites", len(groups))),
        "__METADATA__": metadata, "__METRICS__": cards, "__EVENT_NOTE__": event_note, "__CLOUDS__": cloud_panel,
        "__PROVENANCE__": f"<p>{escape(str(source_note))}</p>" if source_note else "",
        "__DATA__": data,
    }
    html = _HTML
    for token, replacement in replacements.items():
        html = html.replace(token, replacement)
    output = Path(output_path)
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(html, encoding="utf-8")
    return output


_HTML = r"""<!doctype html>
<html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width, initial-scale=1">
<title>Earth observation mission explorer</title><style>
:root{color-scheme:dark;--text:#e8eef5;--muted:#a9b8c9;--edge:#293c52;--panel:#101c2c;--mint:#80e1c2}*{box-sizing:border-box}body{margin:0;background:#091321;color:var(--text);font:15px/1.6 -apple-system,BlinkMacSystemFont,"Segoe UI",sans-serif}main{max-width:1440px;margin:auto;padding:42px 32px 28px}.eyebrow{font-size:12px;letter-spacing:.16em;text-transform:uppercase;color:var(--mint);font-weight:650}h1{font-size:clamp(30px,4vw,48px);line-height:1.12;letter-spacing:-.04em;font-weight:620;margin:10px 0 14px}.intro{color:var(--muted);max-width:830px;margin-bottom:22px}.metadata{display:flex;flex-wrap:wrap;gap:8px;margin:0 0 26px}.metadata span{font-size:12px;color:#c0cfdd;border:1px solid var(--edge);padding:5px 11px;border-radius:5px}.metrics{display:grid;grid-template-columns:repeat(5,minmax(0,1fr));gap:14px;margin-bottom:22px}.metric{border-top:1px solid #55728d;padding:13px 4px 10px}.metric span,.metric small{display:block;color:var(--muted);font-size:12px}.metric strong{display:block;font-size:29px;font-weight:550;letter-spacing:-.035em;margin:4px 0}.panel{background:var(--panel);border:1px solid var(--edge);border-radius:10px;overflow:hidden;margin-bottom:18px;min-width:0}.split{display:grid;grid-template-columns:1.35fr 1fr;gap:18px;align-items:stretch}.panel h2{font-size:18px;line-height:1.4;font-weight:550;letter-spacing:-.02em;margin:20px 23px 6px}.caption{font-size:12px;color:var(--muted);padding:0 23px 18px;margin:0}.legend{display:flex;gap:17px;flex-wrap:wrap;padding:4px 23px 7px;font-size:12px;color:var(--muted)}.legend span::before{content:"";display:inline-block;width:9px;height:9px;border-radius:50%;background:var(--c);margin-right:6px}.legend .diamond::before{border-radius:0;transform:rotate(45deg)}canvas{display:block;width:100%;outline:none}#globe{height:465px;cursor:grab;touch-action:none}#globe:active{cursor:grabbing}#targets{height:410px}#timeline{height:240px}.controls{display:flex;gap:14px;align-items:center;margin:0 23px 14px}.controls button{background:#203a50;color:#e6f4ff;border:1px solid #48657e;padding:7px 14px;min-width:74px;border-radius:5px;font:inherit;font-size:12px;cursor:pointer}.controls button:focus-visible,input:focus-visible{outline:2px solid #80e1c2;outline-offset:3px}.controls input{flex:1;min-width:40px;accent-color:#80e1c2}.clock{font-variant-numeric:tabular-nums;font-size:12px;min-width:95px;color:#d6e4ee}.live{font-size:11px;color:var(--muted);padding:0 23px 12px;min-height:25px}.notes{border-top:1px solid var(--edge);margin-top:20px;padding-top:18px;display:grid;grid-template-columns:1fr 1fr;gap:32px;color:var(--muted);font-size:12px}.notes strong{color:var(--text);font-weight:550}.notes p{margin:7px 0}footer{font-size:11px;color:#8091a4;padding-top:18px}.tooltip{position:fixed;z-index:5;pointer-events:none;display:none;max-width:290px;background:#1b2c40;border:1px solid #526a84;border-radius:6px;padding:9px 12px;font-size:12px;box-shadow:0 5px 25px #0006;white-space:pre-line}noscript{display:block;color:#ffc56e;padding:18px}@media(max-width:1050px){.split{grid-template-columns:1fr}.metrics{grid-template-columns:repeat(3,1fr)}#targets{height:360px}}@media(max-width:600px){main{padding:26px 12px}.metrics{grid-template-columns:repeat(2,1fr);gap:8px}.notes{grid-template-columns:1fr;gap:12px}#globe{height:370px}.controls{gap:8px}.clock{min-width:80px}}
</style></head><body><main>
<header><div class="eyebrow">Mission analysis · __SATELLITES__ satellites</div><h1>Earth observation mission</h1><p class="intro">Explore the propagated constellation, its daylight acquisition opportunities and its ground contacts. Every position and event comes from the exported mission analysis.</p><div class="metadata">__METADATA__</div></header>
<section class="metrics" aria-label="Mission metrics">__METRICS__</section>
__CLOUDS__
<div class="split"><section class="panel" aria-label="Interactive orbital globe"><h2>Orbital motion</h2><div id="satellite-legend" class="legend"></div><canvas id="globe" tabindex="0" aria-label="Earth-fixed globe. Drag to rotate; arrow keys rotate the view."></canvas><div class="controls"><button id="play" aria-label="Play or pause mission animation">Play</button><input id="time" type="range" min="0" value="0" aria-label="Mission time"><span id="clock" class="clock"></span></div><div id="live" class="live" aria-live="polite"></div><p class="caption">Drag to rotate; arrow keys also rotate the globe. Trails show up to the preceding 100 minutes. The surface is the WGS84 ellipsoid with Natural Earth coastline outlines and a coordinate graticule; terrain is not modelled. Display samples are decimated; analysis tables retain full resolution.</p></section>
<section class="panel" aria-label="Ground acquisition target grid"><h2>Acquisition grid</h2><div class="legend"><span style="--c:#82cfc7">Observed targets · brighter means more acquisitions</span><span style="--c:#b4bdc8">Crosses · no acquisition in horizon</span></div><canvas id="targets" aria-label="Acquisition target coordinates and counts"></canvas><p class="caption">Each marker is a configured ground target. Hover for counts and sampled gap information. Coordinate samples do not establish continuous-area coverage.</p></section></div>
<section class="panel" aria-label="Acquisition and contact windows"><h2>Observation and downlink opportunities</h2><canvas id="timeline" aria-label="Timeline of acquisitions and ground contacts"></canvas><p class="caption">__EVENT_NOTE__ Hover for event start and duration; click the timeline to move the globe. Acquisition and contact opportunities can overlap. Scheduling and radio-link feasibility are separate analyses.</p></section>
<section class="notes"><div><strong>How to read this analysis</strong><p>Acquisitions are clear-sky geometric opportunities with the configured target daylight test. Cloud realizations are not painted onto the globe. A contact is a station elevation opportunity, not evidence of available downlink service.</p><p>Interval boundaries depend on the exported sampling and crossing interpolation. Gaps truncated by the analysis boundary are censored; they do not establish a guaranteed revisit bound beyond this horizon.</p></div><div><strong>Scope and reproducibility</strong><p>The globe uses Earth-fixed coordinates. The Earth stays fixed while the satellites move. Animation decimation does not alter acquisition counts or gap statistics.</p>__PROVENANCE__<p>Units: kilometres for globe geometry and altitude, degrees for geographic coordinates, hours from the stated epoch for event times.</p></div></section><footer>Offline document · Embedded data and renderer · Natural Earth public-domain coastlines · No remote scripts, map tiles or fonts.</footer><noscript>Enable JavaScript to explore the mission. The source CSV and JSON files contain the underlying results.</noscript>
</main><div id="tooltip" class="tooltip" role="tooltip"></div><script type="application/json" id="mission-data">__DATA__</script><script>
"use strict";
const DATA=JSON.parse(document.getElementById("mission-data").textContent);
const COLORS=["#67b7ff","#ffc56e","#b9a0ff","#68d5b5"], MINT="#80e1c2", MUTED="#a9b8c9";
const globe=document.getElementById("globe"), targetCanvas=document.getElementById("targets"), timeline=document.getElementById("timeline"), tooltip=document.getElementById("tooltip"), slider=document.getElementById("time");
const groups=DATA.satellites, mission=DATA.summary.mission||{}, q=1-1/298.257223563;
const times=groups[0].states.map(s=>s[0]); slider.max=String(times.length-1); slider.step="1";
let frame=0, playing=false, lastAnimation=0, centerLon=-105*Math.PI/180, centerLat=30*Math.PI/180;
let markerHits=[],targetHits=[],eventHits=[],drag=null;
function n(value,digits=1){return value===null||value===undefined?"—":Number(value).toLocaleString("en-US",{maximumFractionDigits:digits,minimumFractionDigits:digits});}
function elapsed(seconds){const days=Math.floor(seconds/86400), hours=(seconds-days*86400)/3600;return days?`${days}d ${n(hours,1)}h`:`${n(hours,2)} h`;}
function surface(lat,lon){let p=lat*Math.PI/180,l=lon*Math.PI/180,e2=1-q*q,N=1/Math.sqrt(1-e2*Math.sin(p)**2);return [N*Math.cos(p)*Math.cos(l),N*Math.cos(p)*Math.sin(l),N*(1-e2)*Math.sin(p)];}
function project(point,r,cx,cy){const [x,y,z]=point,cl=Math.cos(centerLon),sl=Math.sin(centerLon),cp=Math.cos(centerLat),sp=Math.sin(centerLat);return [cx+r*(-sl*x+cl*y),cy-r*(-sp*cl*x-sp*sl*y+cp*z),cp*cl*x+cp*sl*y+sp*z];}
function setup(canvas){const rect=canvas.getBoundingClientRect(),dpr=Math.min(window.devicePixelRatio||1,2);canvas.width=Math.round(rect.width*dpr);canvas.height=Math.round(rect.height*dpr);const c=canvas.getContext("2d");c.scale(dpr,dpr);c.clearRect(0,0,rect.width,rect.height);return [c,rect.width,rect.height];}
function line(c,points,color,width=1){c.beginPath();let pen=false;for(const p of points){if(!p){pen=false;continue;}if(pen)c.lineTo(p[0],p[1]);else c.moveTo(p[0],p[1]);pen=true;}c.strokeStyle=color;c.lineWidth=width;c.stroke();}
function pointShape(c,x,y,color,diamond=false,size=5){c.beginPath();if(diamond){c.moveTo(x,y-size);c.lineTo(x+size,y);c.lineTo(x,y+size);c.lineTo(x-size,y);c.closePath();}else c.arc(x,y,size,0,2*Math.PI);c.fillStyle=color;c.fill();}
function stateAt(states,time){let lo=0,hi=states.length-1;while(lo<hi){let mid=Math.ceil((lo+hi)/2);if(states[mid][0]<=time)lo=mid;else hi=mid-1;}return [states[lo],lo];}
function drawGlobe(){const [c,w,h]=setup(globe),r=Math.min(w*.39,h*.44),cx=w/2,cy=h*.5,polar=Math.sqrt(Math.sin(centerLat)**2+q*q*Math.cos(centerLat)**2);
const glow=c.createRadialGradient(cx,cy,r*.6,cx,cy,r*1.13);glow.addColorStop(0,"#143349");glow.addColorStop(.88,"#183b53");glow.addColorStop(1,"#101c2c");c.fillStyle=glow;c.fillRect(0,0,w,h);
c.beginPath();c.ellipse(cx,cy,r,r*polar,0,0,Math.PI*2);const sea=c.createRadialGradient(cx-r*.32,cy-r*.3,r*.1,cx,cy,r*1.1);sea.addColorStop(0,"#224a62");sea.addColorStop(1,"#0a2235");c.fillStyle=sea;c.fill();c.strokeStyle="#55758c";c.lineWidth=1;c.stroke();
for(let lat=-60;lat<=60;lat+=30){const p=[];for(let lon=-180;lon<=180;lon+=3){let v=project(surface(lat,lon),r,cx,cy);p.push(v[2]>=0?v:null);}line(c,p,lat===0?"#4a6b80":"#31536b",lat===0?1.1:.7);}
for(let lon=-180;lon<180;lon+=30){const p=[];for(let lat=-90;lat<=90;lat+=3){let v=project(surface(lat,lon),r,cx,cy);p.push(v[2]>=0?v:null);}line(c,p,"#31536b",.7);}
for(const segment of DATA.coastlines||[]){const points=segment.map(([lon,lat])=>{const p=project(surface(lat,lon),r,cx,cy);return p[2]>=0?p:null;});line(c,points,"#abc9bfaa",.85);}
for(const target of DATA.targets){const p=project(surface(target.latitude_deg,target.longitude_deg),r,cx,cy);if(p[2]>0){c.fillStyle="#c5d5df99";c.fillRect(p[0]-1,p[1]-1,2,2);}}
if(Number.isFinite(mission.station_latitude_deg)&&Number.isFinite(mission.station_longitude_deg)){const p=project(surface(mission.station_latitude_deg,mission.station_longitude_deg),r,cx,cy);if(p[2]>0){pointShape(c,p[0],p[1],MINT,true,5);c.fillStyle=MINT;c.font="11px system-ui";c.fillText(mission.station_name||"Ground station",p[0]+10,p[1]-7);}}
markerHits=[];const t=times[frame],live=[];
for(let gi=0;gi<groups.length;gi++){const group=groups[gi],[s,index]=stateAt(group.states,t),color=COLORS[gi%COLORS.length],trail=[];
for(let j=Math.max(0,index-250);j<=index;j++){const row=group.states[j];if(row[0]<s[0]-6000)continue;const p=project([row[1]/6378137,row[2]/6378137,row[3]/6378137],r,cx,cy);const outside=((p[0]-cx)/r)**2+((p[1]-cy)/(r*polar))**2>1.003;trail.push(p[2]>=0||outside?p:null);}line(c,trail,color,1.8);
const p=project([s[1]/6378137,s[2]/6378137,s[3]/6378137],r,cx,cy),outside=((p[0]-cx)/r)**2+((p[1]-cy)/(r*polar))**2>1.003;
if(p[2]>=0||outside){c.shadowBlur=13;c.shadowColor=color;pointShape(c,p[0],p[1],color,gi%2===1,5.5);c.shadowBlur=0;c.fillStyle=color;c.font="12px system-ui";c.fillText(`S${group.satellite}`,p[0]+10,p[1]+4);markerHits.push({x:p[0],y:p[1],text:`Satellite ${group.satellite}\nT + ${elapsed(s[0])}\n${n(s[4],2)}° N · ${n(s[5],2)}° E\n${n(s[6])} km · ${s[7]?"Sunlit":"Eclipse"}`});}
live.push(`S${group.satellite}: ${n(s[6])} km · ${s[7]?"sunlit":"eclipse"}`);}
c.fillStyle=MUTED;c.font="11px system-ui";c.fillText("Earth-fixed · WGS84",20,h-15);document.getElementById("clock").textContent=elapsed(t);document.getElementById("live").textContent=live.join("     /     ");}
function drawTargets(){const [c,w,h]=setup(targetCanvas),pad={left:58,right:24,top:25,bottom:60},x0=pad.left,x1=w-pad.right,y0=pad.top,y1=h-pad.bottom,lonMin=-127,lonMax=-64,latMin=23,latMax=51;
const px=lon=>x0+(lon-lonMin)/(lonMax-lonMin)*(x1-x0),py=lat=>y1-(lat-latMin)/(latMax-latMin)*(y1-y0);c.font="11px system-ui";c.textAlign="center";
for(let lon=-120;lon<=-70;lon+=10){const x=px(lon);line(c,[[x,y0],[x,y1]],"#26394d");c.fillStyle=MUTED;c.fillText(String(lon),x,y1+21);}c.textAlign="right";
for(let lat=25;lat<=50;lat+=5){const y=py(lat);line(c,[[x0,y],[x1,y]],"#26394d");c.fillStyle=MUTED;c.fillText(String(lat),x0-10,y+4);}c.textAlign="center";c.fillText("Longitude (° E)",(x0+x1)/2,h-13);c.save();c.translate(16,(y0+y1)/2);c.rotate(-Math.PI/2);c.fillText("Latitude (° N)",0,0);c.restore();targetHits=[];
const maxCount=Math.max(1,...DATA.targets.map(t=>t.access_count||0));
for(const target of DATA.targets){const x=px(target.longitude_deg),y=py(target.latitude_deg),count=target.access_count||0,frac=count/maxCount;if(count){pointShape(c,x,y,`rgb(${Math.round(67+80*frac)},${Math.round(116+105*frac)},${Math.round(158+54*frac)})`,false,4.5);}else{line(c,[[x-4,y-4],[x+4,y+4]],"#b4bdc8",1.5);line(c,[[x-4,y+4],[x+4,y-4]],"#b4bdc8",1.5);}targetHits.push({x,y,text:`${target.target}\n${n(target.latitude_deg,2)}° N · ${n(target.longitude_deg,2)}° E\n${count} acquisition windows\nLargest sampled gap: ${n(target.max_gap_days,3)} days`});}
if(!DATA.targets.length){c.fillStyle=MUTED;c.fillText("No target coordinates exported",w/2,h/2);}}
function drawTimeline(){const [c,w,h]=setup(timeline),left=Math.min(195,w*.36),right=25,top=30,bottom=48,maxTime=Math.max(1,DATA.duration_s),lanes=[];for(const g of groups){lanes.push({satellite:g.satellite,kind:"acquisition",label:`S${g.satellite} acquisition`,color:COLORS[lanes.length/2%COLORS.length]});lanes.push({satellite:g.satellite,kind:"contact",label:`S${g.satellite} ground contact`,color:MINT});}
const px=s=>left+s/maxTime*(w-left-right),row=i=>top+(i+.5)*(h-top-bottom)/Math.max(1,lanes.length);c.font="11px system-ui";c.textAlign="center";
for(let i=0;i<=6;i++){const time=i*maxTime/6,x=px(time);line(c,[[x,top],[x,h-bottom]],"#26394d");c.fillStyle=MUTED;c.fillText(n(time/3600,0),x,h-bottom+21);}c.fillText("Elapsed hours from mission epoch",(left+w-right)/2,h-9);eventHits=[];
for(let i=0;i<lanes.length;i++){const lane=lanes[i],y=row(i);c.textAlign="right";c.fillStyle="#c4d2df";c.fillText(lane.label,left-14,y+4);line(c,[[left,y],[w-right,y]],"#26394d");const events=lane.kind==="contact"?DATA.contacts:DATA.accesses;c.lineWidth=lane.kind==="contact"?8:5;
for(const event of events){if(event.satellite!==lane.satellite)continue;const a=px(event.start_s),b=px(event.end_s);line(c,[[a,y],[Math.max(a+.8,b),y]],lane.color,lane.kind==="contact"?8:5);eventHits.push({x0:a,x1:Math.max(a+2,b),y,event,lane});}}
line(c,[[px(times[frame]),top-7],[px(times[frame]),h-bottom+2]],"#e5edf5bb",1.2);timeline.dataset.left=String(left);timeline.dataset.right=String(right);}
function draw(){drawGlobe();drawTimeline();}
function hover(event,hits,interval=false){const rect=event.currentTarget.getBoundingClientRect(),x=event.clientX-rect.left,y=event.clientY-rect.top;let hit=interval?hits.find(p=>Math.abs(p.y-y)<9&&x>=p.x0-2&&x<=p.x1+2):hits.find(p=>Math.hypot(p.x-x,p.y-y)<11);if(!hit){tooltip.style.display="none";return;}tooltip.textContent=interval?`${hit.lane.label}\n${hit.event.target}\nStart T + ${elapsed(hit.event.start_s)}\nDuration ${n(hit.event.duration_s)} s`:hit.text;tooltip.style.display="block";tooltip.style.left=Math.min(window.innerWidth-300,event.clientX+14)+"px";tooltip.style.top=Math.max(8,Math.min(window.innerHeight-125,event.clientY+12))+"px";}
function setFrame(index){frame=Math.max(0,Math.min(times.length-1,index));slider.value=String(frame);draw();}
slider.addEventListener("input",()=>setFrame(Number(slider.value)));document.getElementById("play").addEventListener("click",()=>{playing=!playing;document.getElementById("play").textContent=playing?"Pause":"Play";});
function tick(now){if(playing&&now-lastAnimation>75){if(frame>=times.length-1)setFrame(0);else setFrame(Math.min(times.length-1,frame+Math.max(1,Math.floor(times.length/500))));lastAnimation=now;}requestAnimationFrame(tick);}requestAnimationFrame(tick);
globe.addEventListener("pointerdown",event=>{drag=[event.clientX,event.clientY,centerLon,centerLat];globe.setPointerCapture(event.pointerId);});globe.addEventListener("pointermove",event=>{if(drag){centerLon=drag[2]-(event.clientX-drag[0])*.006;centerLat=Math.max(-1.45,Math.min(1.45,drag[3]+(event.clientY-drag[1])*.006));drawGlobe();tooltip.style.display="none";}else hover(event,markerHits);});globe.addEventListener("pointerup",()=>{drag=null;});globe.addEventListener("pointercancel",()=>{drag=null;});
globe.addEventListener("keydown",event=>{if(!["ArrowLeft","ArrowRight","ArrowUp","ArrowDown"].includes(event.key))return;event.preventDefault();centerLon+=event.key==="ArrowLeft"?-.1:event.key==="ArrowRight"?.1:0;centerLat=Math.max(-1.45,Math.min(1.45,centerLat+(event.key==="ArrowUp"?.1:event.key==="ArrowDown"?-.1:0)));drawGlobe();});
targetCanvas.addEventListener("pointermove",event=>hover(event,targetHits));timeline.addEventListener("pointermove",event=>hover(event,eventHits,true));for(const canvas of [globe,targetCanvas,timeline])canvas.addEventListener("pointerleave",()=>tooltip.style.display="none");
timeline.addEventListener("click",event=>{const rect=timeline.getBoundingClientRect(),left=Number(timeline.dataset.left),right=Number(timeline.dataset.right),t=Math.max(0,Math.min(DATA.duration_s,(event.clientX-rect.left-left)/(rect.width-left-right)*DATA.duration_s));let index=0;while(index+1<times.length&&times[index+1]<=t)index++;setFrame(index);});
const legend=document.getElementById("satellite-legend");groups.forEach((group,index)=>{const item=document.createElement("span");item.textContent=`Satellite ${group.satellite}`;item.style.setProperty("--c",COLORS[index%COLORS.length]);if(index%2)item.className="diamond";legend.append(item);});const stationItem=document.createElement("span");stationItem.textContent=mission.station_name||"Ground station";stationItem.style.setProperty("--c",MINT);stationItem.className="diamond";legend.append(stationItem);
new ResizeObserver(()=>{drawTargets();draw();}).observe(document.querySelector("main"));drawTargets();draw();
</script></body></html>
"""
