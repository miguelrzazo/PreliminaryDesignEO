function results = runRgtStudy(config)
%RUNRGTSTUDY Preserve the original J2/J3 SSO rates and RGT bisection criterion.
% Independent sampled rates and repeat-day residuals are array operations.
% Root refinement remains sequential and retains the original tolerances.
validateattributes(config.altitudes_km,{'numeric'},{'vector','nonempty','real','finite','positive'});
validateattributes(config.maxRepeatDays,{'numeric'},{'scalar','integer','positive','finite'});
validateattributes(config.maxIterations,{'numeric'},{'scalar','integer','positive','finite'});
validateattributes(config.rootTolerance,{'numeric'},{'scalar','real','positive','finite'});
validateattributes(config.rootBounds_m,{'numeric'},{'vector','numel',2,'real','positive','finite'});
assert(config.rootBounds_m(1)<config.rootBounds_m(2),'eo:RgtBounds','Root bounds must be increasing.');
altitudes_km=config.altitudes_km(:).';
[inclinations_deg,nodalRates_rad_s,precessionRates_rad_s]=rgtRates(altitudes_km*1e3,config);
nodalPeriods_s=2*pi./nodalRates_rad_s;
ratios=nodalRates_rad_s./(config.earthRotation_rad_s-precessionRates_rad_s);
repeatDays=1:config.maxRepeatDays;
% Enumerate only integer orbit-count/day-count candidates in the original order.
candidates=cell(config.maxRepeatDays,1);
for dayCount=repeatDays
    orbitCounts=ceil(dayCount*min(ratios)):floor(dayCount*max(ratios));
    candidates{dayCount}=[orbitCounts(:),repmat(dayCount,numel(orbitCounts),1)];
end
candidates=vertcat(candidates{:});
solutions=zeros(size(candidates,1),6); solutionCount=0;
for caseIndex=1:size(candidates,1)
    orbitCount=candidates(caseIndex,1); dayCount=candidates(caseIndex,2);
    trialAltitude_m=solveRepeatAltitude(orbitCount,dayCount,config);
    if isnan(trialAltitude_m), continue; end
    [inclination_deg,nodalRate_rad_s]=rgtRates(trialAltitude_m,config);
    solutionCount=solutionCount+1;
    solutions(solutionCount,:)=[trialAltitude_m/1e3,inclination_deg,orbitCount,dayCount, ...
        2*pi/nodalRate_rad_s,orbitCount/dayCount];
end
solutions=sortrows(solutions(1:solutionCount,:),1);
% One altitude column against one repeat-day row replaces the nested residual sweep.
orbitCounts=round(ratios(:).*repeatDays);
errors_rad_s=abs(nodalRates_rad_s(:) - ...
    (config.earthRotation_rad_s-precessionRates_rad_s(:)).*orbitCounts./repeatDays);
residualDrift_km_day=min(errors_rad_s,[],2).' * config.earthRadius_m * ...
    (2*pi/config.earthRotation_rad_s)/1e3;
results.config=config; results.altitudes_km=altitudes_km;
results.inclinations_deg=inclinations_deg; results.nodalPeriods_s=nodalPeriods_s;
results.solutions=solutions; results.residualDrift_km_day=residualDrift_km_day;
if isempty(solutions)
    results.nearestSolution=[];
else
    [~,nearestIndex]=min(abs(solutions(:,1)-config.referenceAltitude_km));
    results.nearestSolution=solutions(nearestIndex,:);
end
if ~config.writeOutputs, return; end
if ~exist(config.outputDirectory,'dir'), mkdir(config.outputDirectory); end
if ~isempty(results.nearestSolution)
    nearest = array2table(results.nearestSolution,'VariableNames', ...
        {'Altitude_km','Inclination_deg','Orbits','Days','NodalPeriod_s','Ratio'});
    writetable(nearest,fullfile(config.outputDirectory,'nearest_rgt.csv'));
end

figure('Position',[100 100 700 420]);
semilogy(altitudes_km,residualDrift_km_day,'b-','LineWidth',1.2); hold on;
for solutionIndex=1:size(solutions,1)
    [~,altitudeIndex]=min(abs(altitudes_km-solutions(solutionIndex,1)));
    plot(solutions(solutionIndex,1),residualDrift_km_day(altitudeIndex),'g.','MarkerSize',18);
end
xline(config.referenceAltitude_km,'r--','LineWidth',1.5); hold off;
xlabel('Altura [km]'); ylabel('Deriva longitudinal por día sideral [km/día sideral]');
title('Condición RGT para órbita SSO con J_2'); grid on; set(gca,'FontSize',11);
exportgraphics(gcf,fullfile(config.outputDirectory,'RGT_residual.pdf'),'ContentType','vector');
end

function [inclination_deg,nodalRate_rad_s,precessionRate_rad_s,frozenEccentricity] = rgtRates(altitudes_m,config)
semimajorAxes_m=config.earthRadius_m+altitudes_m;
meanMotions_rad_s=sqrt(config.earthGravitationalParameter_m3_s2./semimajorAxes_m.^3);
inclinationCosines=-(2*semimajorAxes_m.^(7/2)*config.sunNodalRate_rad_s) / ...
    (3*config.earthJ2*config.earthRadius_m^2*sqrt(config.earthGravitationalParameter_m3_s2));
inclination_deg=acosd(max(-1,min(1,inclinationCosines)));
frozenEccentricity=-config.earthJ3/(2*config.earthJ2) * ...
    (config.earthRadius_m./semimajorAxes_m).*sind(inclination_deg);
semilatusRecta_m=semimajorAxes_m.*(1-frozenEccentricity.^2);
precessionRate_rad_s=-(3*config.earthJ2*meanMotions_rad_s*config.earthRadius_m^2 ./ ...
    (2*semilatusRecta_m.^2)).*cosd(inclination_deg);
perigeeRates_rad_s=(3*config.earthJ2*meanMotions_rad_s*config.earthRadius_m^2 ./ ...
    (4*semilatusRecta_m.^2)).*(4-5*sind(inclination_deg).^2);
anomalyCorrections_rad_s=(3*config.earthJ2*meanMotions_rad_s*config.earthRadius_m^2 .* ...
    sqrt(1-frozenEccentricity.^2)./(4*semilatusRecta_m.^2)).*(2-3*sind(inclination_deg).^2);
nodalRate_rad_s=meanMotions_rad_s+anomalyCorrections_rad_s+perigeeRates_rad_s;
end

function residual = rgtResidual(altitude_m,orbitCount,dayCount,config)
[~,nodalRate_rad_s,precessionRate_rad_s]=rgtRates(altitude_m,config);
residual=nodalRate_rad_s-(config.earthRotation_rad_s-precessionRate_rad_s)*orbitCount/dayCount;
end

function trialAltitude_m=solveRepeatAltitude(orbitCount,dayCount,config)
    lowerAltitude_m=config.rootBounds_m(1); upperAltitude_m=config.rootBounds_m(2);
    lowerResidual=rgtResidual(lowerAltitude_m,orbitCount,dayCount,config);
    upperResidual=rgtResidual(upperAltitude_m,orbitCount,dayCount,config);
    if lowerResidual*upperResidual>0, trialAltitude_m=NaN; return; end
    for iterationIndex=1:config.maxIterations
        trialAltitude_m=(lowerAltitude_m+upperAltitude_m)/2;
        trialResidual=rgtResidual(trialAltitude_m,orbitCount,dayCount,config);
        if abs(trialResidual)<config.rootTolerance, break; end
        if lowerResidual*trialResidual<0
            upperAltitude_m=trialAltitude_m; upperResidual=trialResidual;
        else
            lowerAltitude_m=trialAltitude_m; lowerResidual=trialResidual;
        end
    end
end
