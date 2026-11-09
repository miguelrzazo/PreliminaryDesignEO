function results = runFrozenOrbitStudy(config)
%RUNFROZENORBITSTUDY Compare continuous latitude arcs in the J2/J3 mean model.
% Ascending/descending excursions and their joint envelope are separate metrics.
a = config.earthRadius_m+config.selectedAltitude_km*1000;
mu = config.earthGravitationalParameter_m3_s2;
J2 = config.earthJ2; J3 = config.earthJ3;
n = sqrt(mu/a^3);
inc = acos(-2*a^(7/2)*config.sunNodalRate_rad_s / ...
    (3*J2*config.earthRadius_m^2*sqrt(mu)));
eFrozen = -J3/(2*J2)*(config.earthRadius_m/a)*sin(inc);
time_days = (0:config.missionDuration_days).';
p = a*(1-eFrozen^2);
alpha = -3*n*J2*(config.earthRadius_m/p)^2*(5/4*sin(inc)^2-1)*time_days*86400;
latitudes_deg = linspace(config.latitudeBounds_deg(1),config.latitudeBounds_deg(2),241);
caseNames = {'frozen','circular'};
summaryRows = zeros(2,9);
results.config = config;
for caseIndex = 1:2
    if caseIndex == 1
        ex = zeros(size(alpha)); ey = eFrozen*ones(size(alpha));
    else
        ex = eFrozen*sin(alpha); ey = eFrozen*(1-cos(alpha));
    end
    eccentricity = hypot(ex,ey); omega = atan2(ey,ex);
    [minimum_km,maximum_km,perArc_km] = computeAcquisitionHeightEnvelope( ...
        eccentricity,omega,latitudes_deg,rad2deg(inc),a,config.earthRadius_m);
    envelope_km = maximum_km-minimum_km;
    deviation_km = max(abs([minimum_km,maximum_km]-config.selectedAltitude_km),[],2);
    assert(all(isfinite([minimum_km;maximum_km;perArc_km(:)])), ...
        'eo:FrozenResult','Non-finite height result.');
    assert(all(perArc_km <= envelope_km+1e-10,'all'), ...
        'eo:FrozenEnvelope','An arc cannot exceed the joint envelope.');
    evolution = table(time_days,eccentricity,rad2deg(omega),minimum_km,maximum_km, ...
        perArc_km(:,1),perArc_km(:,2),envelope_km,deviation_km, ...
        'VariableNames',{'MissionDay','Eccentricity','ArgumentOfPerigee_deg', ...
        'MinHeight_km','MaxHeight_km','AscendingArcVariation_km', ...
        'DescendingArcVariation_km','CombinedEnvelope_km','MaxDeviationFromNominal_km'});
    results.(caseNames{caseIndex}) = evolution;
    summaryRows(caseIndex,:) = [min(eccentricity),max(eccentricity),min(minimum_km), ...
        max(maximum_km),max(perArc_km,[],1),max(envelope_km), ...
        max(deviation_km),config.selectedAltitude_km];
end
results.summary = array2table(summaryRows,'VariableNames', ...
    {'MinEccentricity','MaxEccentricity','MinHeight_km','MaxHeight_km', ...
    'MaxAscendingArcVariation_km','MaxDescendingArcVariation_km', ...
    'MaxCombinedEnvelope_km','MaxDeviationFromNominal_km','NominalAltitude_km'});
results.summary.Case = caseNames.';
assert(abs(summaryRows(1,5)-summaryRows(1,6))<1e-9, ...
    'eo:FrozenSymmetry','Frozen northbound and southbound arcs must be symmetric.');
[~,~,~,profile_km] = computeAcquisitionHeightEnvelope(eFrozen,pi/2, ...
    latitudes_deg,rad2deg(inc),a,config.earthRadius_m);
results.profile = table(latitudes_deg(:),profile_km(1,1:241).', ...
    fliplr(profile_km(1,242:end)).','VariableNames', ...
    {'Latitude_deg','AscendingHeight_km','DescendingHeight_km'});
if ~config.writeOutputs, return; end
if ~exist(config.outputDirectory,'dir'), mkdir(config.outputDirectory); end
writetable(results.profile,fullfile(config.outputDirectory,'acquisition_arcs.csv'));
writetable(results.frozen,fullfile(config.outputDirectory,'frozen_evolution.csv'));
writetable(results.circular,fullfile(config.outputDirectory,'circular_evolution.csv'));
writetable(results.summary,fullfile(config.outputDirectory,'frozen_comparison.csv'));
end
