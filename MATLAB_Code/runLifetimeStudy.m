function results = runLifetimeStudy(config)
%RUNLIFETIMESTUDY Compute the historical approximation and optionally plot it.
previousPath = path;
pathCleanup = onCleanup(@() path(previousPath));
addpath(fullfile(fileparts(mfilename('fullpath')),'lib'));
results.config = config;
results.altitudes_km = config.altitudes_km;
results.lifetime_days = computeLifetimeApproximation(config.altitudes_km, ...
    config.anchorAltitudes_km,config.anchorLifetimes_days);
if ~config.writeOutputs, return; end
figure('Position',[100 100 800 480]);
plot(results.altitudes_km,results.lifetime_days,'-','LineWidth',2);
set(gca,'YScale','log'); xlabel('Altura orbital (km)');
ylabel('Dias de vida media en orbita'); grid on;
set(gca,'GridLineStyle','-','GridAlpha',0.2,'XMinorGrid','off','YMinorGrid','off');
end
