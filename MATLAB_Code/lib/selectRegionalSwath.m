function [width_km,worstRevisit_days] = selectRegionalSwath(altitude_km,satelliteCount,telescopeCount,widths_km,params,cache)
%SELECTREGIONALSWATH Apply the sampled CONUS requirement before mass ranking.
width_km = NaN; worstRevisit_days = NaN;
for candidate_km = widths_km(:).'
    acquired_km = telescopeCount*ceil(candidate_km*1000/(params.GSD*telescopeCount))*params.GSD/1000;
    key = sprintf('%.17g_', [altitude_km,satelliteCount,acquired_km, ...
        params.Cov_Requirement,params.solapamiento,params.cobertura_nubes]);
    if isKey(cache,key)
        worst = cache(key);
    else
        worst = regionalRevisit(altitude_km,satelliteCount,acquired_km,params);
        cache(key) = worst;
    end
    if worst <= params.Cov_Requirement
        width_km = candidate_km; worstRevisit_days = worst; return;
    end
end
end

function worst = regionalRevisit(altitude_km,satelliteCount,acquired_km,params)
a = 6378.137e3+altitude_km*1000;
n = sqrt(3.986004418e14/a^3);
sunRate = 2*pi/(365.2421897*86400);
inc = acos(-sunRate/(1.5*n*1.082629989052e-3*(6378.137e3/a)^2));
psi = atand((1-params.solapamiento)*acquired_km/(2*altitude_km));
worst = 0;
for latitude = 25:2:49
    revisit = params.regionalRevisitSolver([a,0,inc,0,0],latitude,1,0, ...
        satelliteCount,1,0,1,0,[1,params.Cov_Requirement],psi)/(1-params.cobertura_nubes);
    assert(isfinite(revisit) && revisit > 0,'eo:RegionalRevisit','Invalid regional revisit.');
    worst = max(worst,revisit);
    if worst > params.Cov_Requirement, return; end
end
end
