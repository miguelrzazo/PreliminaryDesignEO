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
calcularFrozen(630, fullfile(out, 'Frozen630'));
previous_output = getenv('TFG_OUTPUT_ROOT');
previous_gsd = getenv('TFG_GSD');
environment_cleanup = onCleanup(@() restore_environment(previous_output, previous_gsd));
setenv('TFG_OUTPUT_ROOT', fullfile(out, 'Pipeline'));
setenv('TFG_GSD', '80');
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
calcularFrozen(design.Altura_km, fullfile(out, 'FrozenSelected'));

% Revisit each sampled latitude; 25 deg N is a hypothesis to check.
latitudes = (25:2:49).';
revisit_days = zeros(size(latitudes));
a = 6378.137e3 + design.Altura_km*1000;
n = sqrt(mu/a^3);
inc = acos(-solar_rate/(1.5*n*J2*(6378.137e3/a)^2));
psi = atand(.95*design.Swath_Requerido_km/(2*design.Altura_km));
for k = 1:numel(latitudes)
    revisit_days(k) = RevisitCalc([a,0,inc,0,0], latitudes(k), ...
        1,0,design.Num_Satelites,1,0,1,0,[1,7],psi)/(1-1/6);
end
writetable(table(latitudes,revisit_days), fullfile(out, 'latitude_revisit.csv'));
fprintf('CHECK sampled worst-latitude revisit: %.9f days\n', max(revisit_days));

% Compare the original propagator with SGP4 using identical initial elements.
analisiscontacto(design, fullfile(out, 'ContactsKeplerian'), 'two-body-keplerian');
analisiscontacto(design, fullfile(out, 'ContactsSGP4'), 'sgp4');
fprintf('REVIEW_CALCULATIONS_COMPLETE\n');
end

function run_master(root)
run(fullfile(root, 'MATLAB_Code', 'master.m'));
end

function restore_environment(output, gsd)
setenv('TFG_OUTPUT_ROOT', output);
setenv('TFG_GSD', gsd);
end
