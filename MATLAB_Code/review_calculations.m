function review_calculations
% Run the tutor-review calculations without replacing published thesis assets.
script_dir = fileparts(mfilename('fullpath'));
root = fileparts(script_dir);
out = fullfile(root, 'review-results');
if ~exist(out, 'dir'), mkdir(out); end
addpath(script_dir, fullfile(script_dir, 'lib'));
set(groot, 'defaultFigureVisible', 'off');
diary(fullfile(out, 'matlab.log'));
diary_cleanup = onCleanup(@() diary('off'));
fprintf('MATLAB %s\n', version);
disp(ver);

% Check the corrected SSO relation with the independent nodal-rate formula.
a = 6378.137e3 + 630e3;
mu = 3.986004418e14; J2 = 1.082629989052e-3;
solar_rate = 2*pi/(365.2421897*86400);
n = sqrt(mu/a^3);
inc = acos(-solar_rate/(1.5*n*J2*(6378.137e3/a)^2));
assert(abs(-1.5*n*J2*(6378.137e3/a)^2*cos(inc)-solar_rate) < 1e-18);
fprintf('CHECK SSO(630 km): %.8f deg\n', rad2deg(inc));
revisit = CoverageRevisitCalc(80, 630, 190, 7, 2048, 10, ...
    2, 2, 2, 'Refractivo', .05, 1/6, 3);
assert(isscalar(revisit) && (isnan(revisit) || revisit > 0));
fprintf('CHECK corrected reference revisit: %.9f days\n', revisit);

% Keep the 630-km comparison even if the corrected sweep picks another orbit.
reference_frozen = calcularFrozen(630, fullfile(out, 'Frozen630'));
previous_output = getenv('TFG_OUTPUT_ROOT');
previous_gsd = getenv('TFG_GSD');
previous_plots = getenv('TFG_SKIP_PLOTS');
environment_cleanup = onCleanup(@() restore_environment(previous_output, previous_gsd, previous_plots));
setenv('TFG_OUTPUT_ROOT', fullfile(out, 'Pipeline'));
setenv('TFG_GSD', '80');
setenv('TFG_SKIP_PLOTS', '1');
run_master(root);
solutions = readtable(fullfile(out, 'Pipeline', 'OptimumConfigs', ...
    'reporte_completo_configuraciones.csv'));
solutions = solutions(logical(solutions.Es_Viable), :);
assert(~isempty(solutions), 'Corrected sweep produced no viable designs.');
solutions = sortrows(solutions, 'Masa_Total_Constelacion_kg');
best = solutions(1,:);
writetable(best, fullfile(out, 'selected_design.csv'));
disp(best);
design = table2struct(best);
design.GSD = 80;
selected_frozen = calcularFrozen(design.Altura_km, fullfile(out, 'FrozenSelected'));
rgt = calcularRGT(design.Altura_km, fullfile(out, 'RGT'));
end_of_life = Deorbit(design, fullfile(out, 'EndOfLife'));

% Revisit each sampled latitude; 25 deg N is a hypothesis to check.
latitudes = (25:2:49).';
revisit_days = zeros(size(latitudes));
a = 6378.137e3 + design.Altura_km*1000;
n = sqrt(mu/a^3);
inc = acos(-solar_rate/(1.5*n*J2*(6378.137e3/a)^2));
psi = atand(.95*design.Swath_Adquirido_km/(2*design.Altura_km));
for k = 1:numel(latitudes)
    revisit_days(k) = RevisitCalc([a,0,inc,0,0], latitudes(k), ...
        1,0,design.Num_Satelites,1,0,1,0,[1,7],psi)/(1-1/6);
end
assert(all(isfinite(revisit_days)) && all(revisit_days <= 7), ...
    'Selected design fails the sampled regional revisit requirement.');
writetable(table(latitudes,revisit_days), fullfile(out, 'latitude_revisit.csv'));
fprintf('CHECK sampled worst-latitude revisit: %.9f days\n', max(revisit_days));

% Compare the original propagator with SGP4 using identical initial elements.
keplerian = analisiscontacto(design, fullfile(out, 'ContactsKeplerian'), 'two-body-keplerian');
sgp4 = analisiscontacto(design, fullfile(out, 'ContactsSGP4'), 'sgp4');
report = struct('matlab_version', version, 'source_commit', getenv('GITHUB_SHA'), ...
    'selected', table2struct(best), ...
    'reference_frozen', table2struct(reference_frozen), ...
    'selected_frozen', table2struct(selected_frozen), 'rgt', table2struct(rgt), ...
    'end_of_life', table2struct(end_of_life), ...
    'latitudes_deg', latitudes, 'revisit_days', revisit_days, ...
    'keplerian', table2struct(keplerian), 'sgp4', table2struct(sgp4));
report_json = jsonencode(report);
fid = fopen(fullfile(out, 'summary.json'), 'w');
assert(fid >= 0, 'Cannot write calculation summary.');
fprintf(fid, '%s\n', report_json); fclose(fid);
fprintf('REVIEW_SUMMARY %s\n', report_json);
fprintf('REVIEW_CALCULATIONS_COMPLETE\n');
end

function run_master(root)
run(fullfile(root, 'MATLAB_Code', 'master.m'));
end

function restore_environment(output, gsd, plots)
setenv('TFG_OUTPUT_ROOT', output);
setenv('TFG_GSD', gsd);
setenv('TFG_SKIP_PLOTS', plots);
end
