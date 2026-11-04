%% Editable inputs — historical approximate orbital lifetimes
% Preserve the empirical anchors; these values are not propagation results.
config = struct();
config.altitudes_km = 400:50:2000;
config.anchorAltitudes_km = [400 500 700 900 2000];
config.anchorLifetimes_days = [365 3650 36500 365000 3650000];
config.writeOutputs = true;
results = runLifetimeStudy(config);
