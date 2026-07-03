% GSD Sensitivity: Mass breakdown area plot + Dmin overlay
script_dir = fileparts(mfilename('fullpath'));
repo_root = fileparts(script_dir);
results_root = fullfile(repo_root, 'Compiled data', 'revision_tutor');

summary = readtable(fullfile(results_root, 'GSD_sensitivity_summary.csv'));
summary = sortrows(summary, 'GSD_m');

gsd = summary.GSD_m;
masa_seca = summary.Masa_Seca_Satelite_kg;
masa_comb = summary.Masa_Combustible_kg;
masa_total = masa_seca + masa_comb;
d_pupila = summary.Diametro_Pupila_mm;

fig = figure('Visible', 'off', 'Position', [100, 100, 800, 500]);

% Left y-axis: mass breakdown as stacked area fill
yyaxis left;
x_fill = [gsd; flipud(gsd)];
fill([gsd; flipud(gsd)], [masa_seca; zeros(size(masa_seca))], [0.3 0.6 0.9], ...
    'FaceAlpha', 0.6, 'EdgeColor', 'none');
hold on;
fill([gsd; flipud(gsd)], [masa_total; flipud(masa_seca)], [0.9 0.4 0.3], ...
    'FaceAlpha', 0.6, 'EdgeColor', 'none');
ylabel('Masa satelite [kg]', 'Interpreter', 'latex');
set(gca, 'YColor', 'k');

% Right y-axis: Dmin overlay line
yyaxis right;
plot(gsd, d_pupila, 'k-o', 'LineWidth', 2, 'MarkerSize', 8, ...
    'MarkerFaceColor', 'k');
ylabel('Diametro pupila [mm]', 'Interpreter', 'latex');
set(gca, 'YColor', 'k');

xlabel('GSD [m]', 'Interpreter', 'latex');
title('Sensibilidad GSD: Desglose de Masa y Diametro de Pupila (h = 520 km)', 'Interpreter', 'latex');
legend({'Masa seca', 'Combustible', 'Diametro pupila (linea negra)'}, 'Location', 'northwest', 'Interpreter', 'latex');
grid on;
xlim([25, 125]);

fig.Color = 'white';
set(findall(fig, 'Type', 'Axes'), 'Color', 'white');

% Save PDF to LaTeX chapter directory
latex_dir = fullfile(repo_root, 'Latex_Code', '4.Payload');
if ~exist(latex_dir, 'dir'), mkdir(latex_dir); end
exportgraphics(fig, fullfile(latex_dir, 'GSD_sensitivity.pdf'), 'ContentType', 'vector');

% Also save to results
exportgraphics(fig, fullfile(results_root, 'GSD_sensitivity.pdf'), 'ContentType', 'vector');

close(fig);
disp('GSD sensitivity plot saved.');
