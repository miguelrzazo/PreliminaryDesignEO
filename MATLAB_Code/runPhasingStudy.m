function results = runPhasingStudy(config)
%RUNPHASINGSTUDY Preserve the original mean-motion drift and polar animation.
% This remains the source's mean-anomaly approximation, not a Kepler solver.
plotDefaultsCleanup = configureThesisPlotDefaults();
mu_m3_s2 = config.earthGravitationalParameter_m3_s2;
earthRadius_m = config.earthRadius_m;
finalRadius_m = earthRadius_m + config.finalAltitude_m;
phaseApogeeRadius_m = earthRadius_m + config.phaseApogeeAltitude_m;
phaseSemimajorAxis_m = (finalRadius_m+phaseApogeeRadius_m)/2;
phaseEccentricity = (phaseApogeeRadius_m-finalRadius_m)/(phaseApogeeRadius_m+finalRadius_m);
finalPeriod_s = 2*pi*sqrt(finalRadius_m^3/mu_m3_s2);
phasePeriod_s = 2*pi*sqrt(phaseSemimajorAxis_m^3/mu_m3_s2);
finalMeanMotion_rad_s = 2*pi/finalPeriod_s;
phaseMeanMotion_rad_s = 2*pi/phasePeriod_s;
driftRate_rad_s = finalMeanMotion_rad_s-phaseMeanMotion_rad_s;
assert(driftRate_rad_s>0,'eo:PhasingDirection','The original positive-drift model requires a higher phasing apogee.');
validateattributes(config.timeStep_s,{'numeric'},{'scalar','positive','finite'});
phasingDuration_s = pi/abs(driftRate_rad_s);
% The source loop samples the first step at/after pi, then increments its clock.
time_s = (0:ceil(phasingDuration_s/config.timeStep_s))*config.timeStep_s;
finalAngles_rad = finalMeanMotion_rad_s*time_s;
phaseAngles_rad = phaseMeanMotion_rad_s*time_s;
phaseRadii_m = (phaseSemimajorAxis_m*(1-phaseEccentricity^2)) ./ ...
    (1+phaseEccentricity*cos(phaseAngles_rad));
results.config=config; results.time_s=time_s;
results.phasingDuration_s=phasingDuration_s; results.finalPeriod_s=finalPeriod_s;
results.phasePeriod_s=phasePeriod_s; results.phaseOrbitCount=phasingDuration_s/phasePeriod_s;
results.finalPositions_m=[finalRadius_m*cos(finalAngles_rad);finalRadius_m*sin(finalAngles_rad)];
results.phasePositions_m=[phaseRadii_m.*cos(phaseAngles_rad);phaseRadii_m.*sin(phaseAngles_rad)];
results.accumulatedPhase_rad=driftRate_rad_s*time_s;
results.displayedDeploymentTime_s=time_s(end)+config.timeStep_s;
if ~config.writeOutputs, return; end
if ~exist(config.outputDirectory,'dir'), mkdir(config.outputDirectory); end
figure('Name','Simulación de Desfase Orbital','NumberTitle','off'); hold on; axis equal; grid on;
[x,y,z]=sphere;
surf(x*earthRadius_m/1e3,y*earthRadius_m/1e3,z*earthRadius_m/1e3, ...
    'FaceColor','#0077be','EdgeColor','none','DisplayName','Tierra');
angles_rad=linspace(0,2*pi,300);
plot(finalRadius_m/1e3*cos(angles_rad),finalRadius_m/1e3*sin(angles_rad),'g--', ...
    'LineWidth',1.5,'DisplayName','Orbita final');
phaseCurve_m=phaseSemimajorAxis_m*(1-phaseEccentricity^2)./(1+phaseEccentricity*cos(angles_rad));
plot(phaseCurve_m/1e3.*cos(angles_rad),phaseCurve_m/1e3.*sin(angles_rad),'b', ...
    'LineWidth',1.5,'DisplayName','Orbita de fase');
firstSatellite=plot(NaN,NaN,'go','MarkerFaceColor','g','DisplayName','Satelite 1');
phaseVehicle=plot(NaN,NaN,'ro','MarkerFaceColor','r','DisplayName','Etapa/Satelite 2');
secondSatellite=plot(NaN,NaN,'mo','MarkerFaceColor','m','DisplayName','Satelite 2 (Desplegado)');
titleHandle=title('Inicio de Mision','Interpreter','latex','FontSize',14);
legend('Location','northeastoutside','AutoUpdate','off','Interpreter','latex');
xlim([-1.2 1.2]*phaseApogeeRadius_m/1e3); ylim([-1.2 1.2]*phaseApogeeRadius_m/1e3);
xlabel('X (km)'); ylabel('Y (km)');
% Rendering is sequential; the trajectories themselves were computed as arrays.
for sampleIndex=1:numel(time_s)
    set(firstSatellite,'XData',results.finalPositions_m(1,sampleIndex)/1e3, ...
        'YData',results.finalPositions_m(2,sampleIndex)/1e3);
    set(phaseVehicle,'XData',results.phasePositions_m(1,sampleIndex)/1e3, ...
        'YData',results.phasePositions_m(2,sampleIndex)/1e3);
    titleHandle.String=sprintf('Fase de Deriva | Tiempo: %.1f h | Desfase: %.1f deg', ...
        time_s(sampleIndex)/3600,rad2deg(results.accumulatedPhase_rad(sampleIndex)));
    drawnow limitrate;
end
titleHandle.String=sprintf('Desfase de 180 deg alcanzado en %.1f horas.',results.displayedDeploymentTime_s/3600);
set(secondSatellite,'XData',results.phasePositions_m(1,end)/1e3, ...
    'YData',results.phasePositions_m(2,end)/1e3,'Visible','on');
set(phaseVehicle,'Visible','off'); hold off;
exportgraphics(gcf,fullfile(config.outputDirectory,'phasing.jpg'),'Resolution',300);
end
