function results = runDeorbitStudy(config)
%RUNDEORBITSTUDY Integrate maintenance cycles and circular-equivalent end-of-life decay.
validateattributes(config.missionDuration_years,{'numeric'},{'scalar','real','finite','positive'});
validateattributes(config.dryMass_kg,{'numeric'},{'scalar','real','finite','positive'});
validateattributes(config.dragArea_m2,{'numeric'},{'scalar','real','finite','positive'});
validateattributes(config.maintenanceDrop_percent,{'numeric'},{'scalar','real','finite','>',0,'<',100});
assert(config.nominalAltitude_km>config.reentryAltitude_km,'eo:DeorbitBounds','Nominal altitude must exceed reentry altitude.');
results.config=config;
nominalAltitude_km=config.nominalAltitude_km; missionDuration_years=config.missionDuration_years;
maintenanceDrop_percent=config.maintenanceDrop_percent; reentryAltitude_km=config.reentryAltitude_km;
deorbitRequirement_days=config.deorbitRequirement_days; dryMass_kg=config.dryMass_kg;
dragArea_m2=1.2*config.dragArea_m2; dragCoefficient=config.dragCoefficient; specificImpulse_s=config.specificImpulse_s;
earthGravitationalParameter_m3_s2 = 3.986004418e14; earthRadius_m = 6378e3; standardGravity_m_s2 = 9.80665;
nominalAltitude_m = nominalAltitude_km * 1e3;
resetAltitude_m = 1.02*nominalAltitude_m;
postMissionHorizon_s = 100*365.25*86400;
maintenanceThreshold_m = nominalAltitude_km * (1 - maintenanceDrop_percent / 100) * 1e3;
missionDuration_s = missionDuration_years * 365.25 * 24 * 3600;
timeSegments_s = {0}; altitudeSegments_m = {resetAltitude_m};
elapsedTime_s = 0; currentAltitude_m = resetAltitude_m; currentMass_kg = dryMass_kg;
totalDeltaV_m_s = 0; impulseCount = 0;
previousPath=path; pathCleanup=onCleanup(@() path(previousPath));
addpath(fullfile(fileparts(mfilename('fullpath')),'lib'));
if isfield(config,'maintenanceFuel_kg')
    maintenanceFuel_kg = config.maintenanceFuel_kg;
else
    mass = computeTotalMass(nominalAltitude_km,dryMass_kg,config.dragArea_m2,0.98,1.02);
    maintenanceFuel_kg = mass.fuelMass_kg;
end
initialMass_kg = dryMass_kg + maintenanceFuel_kg;
currentMass_kg = initialMass_kg;
fprintf('Calculo preliminar:\n  - Masa total inicial: %.2f kg\n\n', initialMass_kg);
fprintf('Iniciando simulacion de la mision de %d a~nos...\n', missionDuration_years);
maintenanceTargetAxis_m = earthRadius_m + maintenanceThreshold_m;
while elapsedTime_s < missionDuration_s
    initialSemimajorAxis_m = earthRadius_m + currentAltitude_m;
    options = odeset('Events', @(t,semimajorAxis_m) decayEvent(t, semimajorAxis_m, maintenanceTargetAxis_m), 'RelTol', 1e-7);
    [integrationTime_s, integratedAxes_m] = ode45(@(t,semimajorAxis_m) decayODE(t, semimajorAxis_m, earthRadius_m, earthGravitationalParameter_m3_s2, dragCoefficient, dragArea_m2, currentMass_kg), ...
                           [0, missionDuration_s - elapsedTime_s], initialSemimajorAxis_m, options);
    timeSegments_s{end+1,1} = elapsedTime_s + integrationTime_s(2:end);
    altitudeSegments_m{end+1,1} = integratedAxes_m(2:end) - earthRadius_m;
    elapsedTime_s = elapsedTime_s + integrationTime_s(end); currentAltitude_m = integratedAxes_m(end) - earthRadius_m;
    if elapsedTime_s >= missionDuration_s; break; end
    impulseCount = impulseCount + 1;
    fprintf('Impulso #%d en t = %.2f a~nos. Altitud: %.1f km\n', impulseCount, elapsedTime_s / (365.25*24*3600), currentAltitude_m/1000);
    initialRadius_m = earthRadius_m + currentAltitude_m; targetRadius_m = earthRadius_m + resetAltitude_m;
    impulseDeltaV_m_s = hohmannDeltaV(initialRadius_m, targetRadius_m, earthGravitationalParameter_m3_s2);
    currentMass_kg = currentMass_kg / exp(impulseDeltaV_m_s / (specificImpulse_s * standardGravity_m_s2));
    currentAltitude_m = resetAltitude_m;
    timeSegments_s{end+1,1} = elapsedTime_s; altitudeSegments_m{end+1,1} = currentAltitude_m;
end
missionTime_s = vertcat(timeSegments_s{:});
missionAltitudes_m = vertcat(altitudeSegments_m{:});
fprintf('Simulacion de mision completada.\n\n');
massBeforeDeorbit_kg = currentMass_kg;
reentryTargetAxis_m = earthRadius_m + reentryAltitude_km * 1e3;
reentryOptions = odeset('Events', @(t,semimajorAxis_m) decayEvent(t, semimajorAxis_m, reentryTargetAxis_m), 'RelTol', 1e-7);
fprintf('Simulando decaimiento natural (sin impulso)...\n');
naturalInitialAxis_m = earthRadius_m + missionAltitudes_m(end);
[naturalTime_s, naturalAxes_m] = ode45(@(t,semimajorAxis_m) decayODE(t, semimajorAxis_m, earthRadius_m, earthGravitationalParameter_m3_s2, dragCoefficient, dragArea_m2, massBeforeDeorbit_kg), ...
                               [0, postMissionHorizon_s], naturalInitialAxis_m, reentryOptions);
assert(naturalAxes_m(end)-earthRadius_m <= reentryAltitude_km*1000+1, ...
    'eo:NaturalDecayHorizon','Natural decay exceeds the integration horizon.');
naturalReentryTime_days = naturalTime_s(end) / (24 * 3600);
fprintf('  - Reentrada natural estimada en %.2f dias (%.2f a~nos).\n\n', naturalReentryTime_days, naturalReentryTime_days/365.25);
fprintf('Buscando impulso minimo de deorbitacion para %.0f dias...\n', deorbitRequirement_days);
deorbitApogeeRadius_m = earthRadius_m + missionAltitudes_m(end);
lowerPerigee_m = reentryAltitude_km * 1e3;
upperPerigee_m = missionAltitudes_m(end);
for iterationIndex = 1:50
    trialPerigee_m = 0.5 * (lowerPerigee_m + upperPerigee_m);
    trialPerigeeRadius_m = earthRadius_m + trialPerigee_m;
    trialTransferAxis_m = 0.5 * (deorbitApogeeRadius_m + trialPerigeeRadius_m);
    [trialTime_s, ~] = ode45(@(t,semimajorAxis_m) decayODE(t, semimajorAxis_m, earthRadius_m, earthGravitationalParameter_m3_s2, dragCoefficient, dragArea_m2, dryMass_kg), ...
                       [0, postMissionHorizon_s], trialTransferAxis_m, reentryOptions);
    trialReentryTime_days = trialTime_s(end) / (24 * 3600);
    if trialReentryTime_days <= deorbitRequirement_days
        lowerPerigee_m = trialPerigee_m;
    else
        upperPerigee_m = trialPerigee_m;
    end
end
deorbitPerigee_m = lowerPerigee_m;
deorbitPerigeeRadius_m = earthRadius_m + deorbitPerigee_m;
finalCircularSpeed_m_s = sqrt(earthGravitationalParameter_m3_s2 / deorbitApogeeRadius_m);
deorbitTransferAxis_m = 0.5 * (deorbitApogeeRadius_m + deorbitPerigeeRadius_m);
transferApogeeSpeed_m_s = sqrt(earthGravitationalParameter_m3_s2 * (2/deorbitApogeeRadius_m - 1/deorbitTransferAxis_m));
deorbitDeltaV_m_s = finalCircularSpeed_m_s - transferApogeeSpeed_m_s;
deorbitFuel_kg = dryMass_kg * (exp(deorbitDeltaV_m_s / (specificImpulse_s * standardGravity_m_s2)) - 1) * 1.1;
fprintf('Maniobra de Deorbitacion:\n');
fprintf('  - Perigeo objetivo: %.1f km\n', deorbitPerigee_m / 1e3);
fprintf('  - Delta-V minimo: %.2f m/s\n', deorbitDeltaV_m_s);
fprintf('  - Combustible reservado (10%% margen): %.3f kg\n', deorbitFuel_kg);
fprintf('Simulando reentrada con impulso...\n');
[activeTime_s, activeAxes_m] = ode45(@(t,semimajorAxis_m) decayODE(t, semimajorAxis_m, earthRadius_m, earthGravitationalParameter_m3_s2, dragCoefficient, dragArea_m2, dryMass_kg), ...
                               [0, postMissionHorizon_s], deorbitTransferAxis_m, reentryOptions);
activeReentryTime_days = activeTime_s(end) / (24 * 3600);
fprintf('  - Reentrada estimada en %.2f dias.\n\n', activeReentryTime_days);
totalTime_s = [missionTime_s; elapsedTime_s + activeTime_s(2:end)];
totalAltitudes_m = [missionAltitudes_m; activeAxes_m(2:end) - earthRadius_m];
time_years = totalTime_s / (365.25 * 24 * 3600);
altitudes_km = totalAltitudes_m / 1000;
results.time_s=totalTime_s; results.altitudes_m=totalAltitudes_m;
results.maintenanceFuel_kg=maintenanceFuel_kg; results.deorbitFuel_kg=deorbitFuel_kg;
results.deorbitDeltaV_m_s=deorbitDeltaV_m_s; results.deorbitPerigee_m=deorbitPerigee_m;
results.naturalReentryTime_days=naturalReentryTime_days; results.activeReentryTime_days=activeReentryTime_days;
results.impulseCount=impulseCount;
results.summary = table(nominalAltitude_km,missionAltitudes_m(end)/1000,impulseCount, ...
    naturalReentryTime_days,deorbitPerigee_m/1000,deorbitDeltaV_m_s,maintenanceFuel_kg, ...
    deorbitFuel_kg,dryMass_kg+maintenanceFuel_kg+deorbitFuel_kg, ...
    config.satelliteCount*(dryMass_kg+maintenanceFuel_kg+deorbitFuel_kg), ...
    'VariableNames',{'NominalAltitude_km','EndAltitude_km','SimulatedMaintenanceManoeuvres', ...
    'NaturalDecay_days','DeorbitPerigee_km','DeorbitDeltaV_ms','MaintenanceFuel_kg', ...
    'DeorbitFuel_kg','TotalSatelliteMass_kg','TotalConstellationMass_kg'});
assert(activeReentryTime_days <= deorbitRequirement_days+0.1, ...
    'eo:DeorbitRequirement','Deorbit time exceeds the model limit.');
if ~config.writeOutputs, return; end
if ~exist(config.outputDirectory,'dir'), mkdir(config.outputDirectory); end
writetable(results.summary,fullfile(config.outputDirectory,'end_of_life.csv'));
figure('Name', 'Evolucion de la Altitud con Reentrada');
hold on; grid on;
plot(time_years, altitudes_km, 'b-', 'LineWidth', 2, ...
     'DisplayName', sprintf('Trayectoria ($M_{seca}$=%dkg, A=%.1fm$^2$)', dryMass_kg, dragArea_m2));
yline(nominalAltitude_km, '--g', 'LineWidth', 1.5, 'DisplayName', 'Altitud Inicial');
yline(maintenanceThreshold_m/1e3, '--r', 'LineWidth', 1.5, 'DisplayName', 'Altitud Umbral');
xline(missionDuration_years, '--m', 'LineWidth', 1.5, 'DisplayName', 'Fin de Mision');
plot(missionDuration_years, altitudes_km(find(time_years >= missionDuration_years, 1)), ...
     'r*', 'MarkerSize', 10, 'LineWidth', 2, ...
     'DisplayName', sprintf('Impulso Deorbit ($\\Delta V$=%.1f m/s)', deorbitDeltaV_m_s));
title(sprintf('Evolucion de la Altitud durante %d a~nos y Reentrada', missionDuration_years), 'FontSize', 16, 'Interpreter', 'latex');
xlabel('Tiempo (a~nos)', 'FontSize', 14, 'Interpreter', 'latex');
ylabel('Altitud (km)', 'FontSize', 14, 'Interpreter', 'latex');
xlim([0, time_years(end) * 1.05]); ylim([0, nominalAltitude_km * 1.05]);
legend('Location', 'best', 'Interpreter', 'latex', 'FontSize', 12);
hold off;
deorbitFigure = figure('Name', 'Decaimiento Post-Mision');
hold on; grid on;
semilogx(naturalTime_s / (24*3600), (naturalAxes_m - earthRadius_m) / 1000, ...
    'k--', 'LineWidth', 2, 'DisplayName', 'Decaimiento Natural (sin impulso)');
semilogx(activeTime_s / (24*3600), (activeAxes_m - earthRadius_m) / 1000, ...
    'r-', 'LineWidth', 2, 'DisplayName', sprintf('Con impulso deorbitado (%.1f m/s)', deorbitDeltaV_m_s));
xline(deorbitRequirement_days, 'b:', 'LineWidth', 1.5, 'DisplayName', 'Limite 365 dias');
title('Decaimiento Post-Mision', 'FontSize', 16, 'Interpreter', 'latex');
xlabel('Tiempo tras fin de mision (dias)', 'FontSize', 14, 'Interpreter', 'latex');
ylabel('Altitud (km)', 'FontSize', 14, 'Interpreter', 'latex');
legend('Location', 'best', 'Interpreter', 'latex', 'FontSize', 12);
hold off;
exportgraphics(deorbitFigure,fullfile(config.outputDirectory,'deorbit.jpg'),'Resolution',300);
end

function axisRate_m_s = decayODE(~, semimajorAxis_m, earthRadius_m, earthGravitationalParameter_m3_s2, dragCoefficient, dragArea_m2, satelliteMass_kg)
    altitude_m = semimajorAxis_m - earthRadius_m; if altitude_m < 0; altitude_m = 0; end
    density_kg_m3 = getAtmosphericDensity(altitude_m);
    axisRate_m_s = -dragCoefficient * dragArea_m2 * density_kg_m3 * sqrt(earthGravitationalParameter_m3_s2 * semimajorAxis_m) / satelliteMass_kg;
end
function [value, isterminal, direction] = decayEvent(~, semimajorAxis_m, targetAxis_m)
    value = semimajorAxis_m - targetAxis_m; isterminal = 1; direction = -1;
end
function deltaV_m_s = hohmannDeltaV(initialRadius_m, targetRadius_m, earthGravitationalParameter_m3_s2)
    initialSpeed_m_s = sqrt(earthGravitationalParameter_m3_s2 / initialRadius_m); targetSpeed_m_s = sqrt(earthGravitationalParameter_m3_s2 / targetRadius_m); transferAxis_m = (initialRadius_m + targetRadius_m) / 2;
    firstTransferSpeed_m_s = sqrt(earthGravitationalParameter_m3_s2 * (2/initialRadius_m - 1/transferAxis_m));
    secondTransferSpeed_m_s = sqrt(earthGravitationalParameter_m3_s2 * (2/targetRadius_m - 1/transferAxis_m));
    deltaV_m_s = abs(firstTransferSpeed_m_s - initialSpeed_m_s) + abs(targetSpeed_m_s - secondTransferSpeed_m_s);
end
function density_kg_m3 = getAtmosphericDensity(altitude_m)
    % Preserve this script's original exponential-atmosphere table and lookup.
    altitude_km = max(altitude_m / 1e3, 0);
    baseAltitudes_km = [0, 25, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120, ...
                 130, 140, 150, 180, 200, 250, 300, 350, 400, 450, ...
                 500, 600, 700, 800, 900, 1000];
    baseDensities_kg_m3 = [1.225, 3.899e-2, 1.774e-2, 3.972e-3, 1.057e-3, ...
                3.206e-4, 8.770e-5, 1.905e-5, 3.396e-6, 5.297e-7, ...
                9.661e-8, 2.438e-8, 8.484e-9, 3.845e-9, 2.070e-9, ...
                5.464e-10, 2.789e-10, 7.248e-11, 2.418e-11, ...
                9.518e-12, 3.725e-12, 1.585e-12, 6.967e-13, ...
                1.454e-13, 3.614e-14, 1.170e-14, 5.245e-15, 3.019e-15];
    scaleHeights_km = [7.249, 6.349, 6.682, 7.554, 8.382, 7.714, 6.549, ...
            5.799, 5.382, 5.877, 7.263, 9.473, 12.636, 16.149, ...
            22.523, 29.740, 37.105, 45.546, 53.628, 53.298, ...
            58.515, 60.828, 63.822, 71.835, 88.667, 124.64, ...
            181.05, 268.00];
    layerIndex = find(altitude_km >= baseAltitudes_km, 1, 'last');
    if isempty(layerIndex)
        layerIndex = 1;
    elseif layerIndex == numel(baseAltitudes_km)
        layerIndex = numel(baseAltitudes_km);
    end
    density_kg_m3 = baseDensities_kg_m3(layerIndex) * exp(-(altitude_km - baseAltitudes_km(layerIndex)) / scaleHeights_km(layerIndex));
end
