% Recompute and plot the existing fixed-height GSD sensitivity samples.
script_dir = fileparts(mfilename('fullpath'));
repo_root = fileparts(script_dir);
results_root = fullfile(repo_root, 'Compiled data');
addpath(fullfile(repo_root, 'MATLAB_Code', 'lib'));

summary = readtable(fullfile(results_root, 'GSD_sensitivity_summary.csv'));
summary = sortrows(summary, 'GSD_m');

for i = 1:height(summary)
    [dry_mass, ~, ~, ~, ~, ~, ~, ~, ~, ~, area_sat, ~] = ...
        calcularMasaSeca(summary.Diametro_Pupila_mm(i), summary.Num_Telescopios(i));
    [~, fuel_final] = calcularMasaTotal(summary.Altura_km(i), dry_mass, area_sat);
    summary.Masa_Seca_Satelite_kg(i) = dry_mass;
    summary.Masa_Combustible_kg(i) = fuel_final;
    summary.Masa_Total_Constelacion_kg(i) = summary.Num_Satelites(i) * (dry_mass + fuel_final);
end
writetable(summary, fullfile(results_root, 'GSD_sensitivity_summary.csv'));

gsd = summary.GSD_m;
dry_constellation = 2 * summary.Masa_Seca_Satelite_kg;
fuel_constellation = 2 * summary.Masa_Combustible_kg;
diameter_mm = summary.Diametro_Pupila_mm;

fig = figure('Visible', 'off', 'Position', [100, 100, 900, 520]);

yyaxis left;
mass_bars = bar(gsd, [dry_constellation, fuel_constellation], 0.72, 'stacked', 'LineWidth', 0.4);
mass_bars(1).FaceColor = [0.25, 0.45, 0.70];  % masa seca: azul grisáceo
mass_bars(2).FaceColor = [0.90, 0.55, 0.20];  % combustible: naranja
mass_bars(1).EdgeColor = [0.15, 0.25, 0.40];
mass_bars(2).EdgeColor = [0.55, 0.30, 0.05];
hold on;
nominal_line = xline(80, ':', 'GSD nominal', 'LabelVerticalAlignment', 'bottom', ...
    'Interpreter', 'latex', 'Color', [0.2 0.2 0.2]);
ylabel('Masa de seleccion de la constelacion [kg]', 'Interpreter', 'latex');
set(gca, 'YColor', 'k');

yyaxis right;
diameter_line = plot(gsd, diameter_mm, 'k-o', 'LineWidth', 2, 'MarkerSize', 7, ...
    'MarkerFaceColor', 'k');
ylabel('Diametro pupila [mm]', 'Interpreter', 'latex');
set(gca, 'YColor', 'k');

xlabel('GSD [m]', 'Interpreter', 'latex');
title(sprintf('Sensibilidad al GSD a altura fija de %g km', summary.Altura_km(1)), ...
    'Interpreter', 'latex', 'Color', 'k');
legend([mass_bars(1), mass_bars(2), nominal_line, diameter_line], ...
    {'Masa seca', 'Combustible de mantenimiento', 'GSD nominal', 'Diametro pupila'}, ...
    'Location', 'northwest', 'Interpreter', 'latex', 'TextColor', 'black', 'Color', 'white');
grid on;
xlim([25, 125]);

fig.Color = 'white';
set(findall(fig, 'Type', 'Axes'), 'Color', 'white', 'XColor', 'black');

latex_dir = fullfile(repo_root, 'Latex_Code', '4.Payload');
exportgraphics(fig, fullfile(latex_dir, 'GSD_sensitivity.pdf'), 'ContentType', 'vector');
exportgraphics(fig, fullfile(results_root, 'GSD_sensitivity.pdf'), 'ContentType', 'vector');
exportgraphics(fig, fullfile(results_root, 'GSD_sensitivity.jpg'), 'Resolution', 300);

close(fig);
fprintf('GSD sensitivity data and plot saved at fixed altitude %.0f km.\n', summary.Altura_km(1));
