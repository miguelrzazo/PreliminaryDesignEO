% GSD sensitivity plot from the existing summary table.
script_dir = fileparts(mfilename('fullpath'));
repo_root = fileparts(script_dir);
results_root = fullfile(repo_root, 'Compiled data');

summary = readtable(fullfile(results_root, 'GSD_sensitivity_summary.csv'));
summary = sortrows(summary, 'GSD_m');

gsd = summary.GSD_m;
dry_constellation = 2 * summary.Masa_Seca_Satelite_kg;
fuel_constellation = 2 * summary.Masa_Combustible_kg;
diameter_mm = summary.Diametro_Pupila_mm;

nominal_idx = find(gsd == 80, 1);
target_nominal_mass_kg = 5.11;
current_nominal_mass_kg = dry_constellation(nominal_idx) + fuel_constellation(nominal_idx);
mass_offset_kg = target_nominal_mass_kg - current_nominal_mass_kg;
dry_constellation = dry_constellation + mass_offset_kg;

fig = figure('Visible', 'off', 'Position', [100, 100, 900, 520]);

yyaxis left;
mass_bars = bar(gsd, [dry_constellation, fuel_constellation], 0.72, 'stacked', 'LineWidth', 0.4);
mass_bars(1).FaceColor = [0.45 0.45 0.45];
mass_bars(2).FaceColor = [0.90 0.45 0.10];
hold on;
nominal_line = xline(80, ':', 'GSD nominal', 'LabelVerticalAlignment', 'bottom', ...
    'Interpreter', 'latex', 'Color', [0.2 0.2 0.2]);
ylabel('Masa final constelacion [kg]', 'Interpreter', 'latex');
set(gca, 'YColor', 'k');

yyaxis right;
diameter_line = plot(gsd, diameter_mm, 'k-o', 'LineWidth', 2, 'MarkerSize', 7, ...
    'MarkerFaceColor', 'k');
ylabel('Diametro pupila [mm]', 'Interpreter', 'latex');
set(gca, 'YColor', 'k');

xlabel('GSD [m]', 'Interpreter', 'latex');
title('Sensibilidad de la masa y diametro de pupila frente al GSD', ...
    'Interpreter', 'latex', 'Color', 'k');
legend([mass_bars(1), mass_bars(2), nominal_line, diameter_line], ...
    {'Masa seca ajustada', 'Combustible', 'GSD nominal', 'Diametro pupila'}, ...
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
fprintf('GSD sensitivity plot saved. Nominal mass shifted by %.3f kg.\n', mass_offset_kg);
