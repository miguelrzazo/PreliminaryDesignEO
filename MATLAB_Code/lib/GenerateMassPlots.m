function GenerateMassPlots(results_table)
% GenerateMassPlots: Crea gráficos de masa para las configuraciones viables.
%
% Descripción:
% Esta función genera dos tipos de gráficos para cada configuración única
% encontrada en la tabla de resultados:
% 1. Comparativa de Masa Seca: Muestra la masa seca estimada. Se guarda
%    en la carpeta 'Masa_seca'.
% 2. Desglose de Masa vs. Altura: Muestra la masa total, seca y de 
%    combustible. Se guarda en la carpeta 'Masa_total'.
%
% Input:
% results_table - Tabla con todos los datos de las soluciones viables.

%% ========================================================================
% 1. CONFIGURACIÓN INICIAL
% ========================================================================
dry_mass_plot_dir = 'Masa_seca';
total_mass_plot_dir = 'Masa_total';
if ~exist(dry_mass_plot_dir, 'dir'), mkdir(dry_mass_plot_dir); end
if ~exist(total_mass_plot_dir, 'dir'), mkdir(total_mass_plot_dir); end

fprintf('Generando gráficas de masa para las soluciones encontradas...\n');

%% ========================================================================
% 2. GENERACIÓN DE GRÁFICOS POR CONFIGURACIÓN
% ========================================================================
unique_configs = unique(results_table(:, {'Num_Satelites', 'Num_Telescopios', 'ID_Detector', 'Tipo_Telescopio'}), 'rows');

% Iterar sobre cada configuración única
for i = 1:height(unique_configs)
    current_config = unique_configs(i, :);
    N_sat = current_config.Num_Satelites;
    N_tel = current_config.Num_Telescopios;
    det_id = current_config.ID_Detector;
    tel_name_cell = current_config.Tipo_Telescopio;
    tel_name = tel_name_cell{1};

    % Filtrar datos para la configuración actual
    config_mask = results_table.Num_Satelites == N_sat & ...
                  results_table.Num_Telescopios == N_tel & ...
                  results_table.ID_Detector == det_id & ...
                  strcmp(results_table.Tipo_Telescopio, tel_name);
    
    data_for_plot = results_table(config_mask, :);
    data_for_plot = sortrows(data_for_plot, 'Altura_km');

    if isempty(data_for_plot), continue; end

    % --- GRÁFICO 1: COMPARATIVA DE MASA SECA ---
    fig1 = figure('Visible', 'off', 'Position', [100, 100, 900, 600]);
    hold on;

    % Plot de las estimaciones base (Thematic Mapper y SEOSAT)
    plot(data_for_plot.Altura_km, data_for_plot.Masa_Seca_TM_kg, '--o', 'LineWidth', 1.5, 'DisplayName', 'Estimacion (Thematic Mapper)');
    plot(data_for_plot.Altura_km, data_for_plot.Masa_Seca_SEOSAT_kg, '--s', 'LineWidth', 1.5, 'DisplayName', 'Estimacion (SEOSAT)');
    
    % Plot de la masa propuesta para la configuración actual
    plot(data_for_plot.Altura_km, data_for_plot.Masa_Seca_Satelite_kg, '-d', 'LineWidth', 2, 'MarkerSize', 5, 'DisplayName', sprintf('Propuesta (%d Tel)', N_tel));
    
    % Si la configuración actual es de 2 telescopios, añade la línea de 1 telescopio para comparar
    if N_tel >= 2
        % Estimar la masa seca para una configuración de 1 telescopio dividiendo por 1.5
        % según la lógica de 'calcularMasaSeca.m'.
        masa_seca_estimada_1_tel = data_for_plot.Masa_Seca_Satelite_kg / 1.5;
    
        % Graficar la masa seca estimada para 1 telescopio con un estilo diferente
        plot(data_for_plot.Altura_km, masa_seca_estimada_1_tel, '--p', 'LineWidth', 1.5, 'DisplayName', 'Promedio (1 Telescopio)');
    end
    
    hold off;
    
    % Configuración del gráfico (título, etiquetas, etc.)
    title_str = sprintf('Comparativa de Masa Seca vs. Altura\nConfig: %d Sat, %d Tel, Det %d, %s', N_sat, N_tel, det_id, tel_name);
    title(title_str, 'Interpreter', 'latex');
    xlabel('Altura Orbital [km]', 'Interpreter', 'latex');
    ylabel('Masa Seca [kg]', 'Interpreter', 'latex');
    legend('Location', 'northwest', 'Interpreter', 'latex');
    grid on;
    axis tight;
    
    % Guardar el gráfico
    fig1.Color = 'white';
    set(findall(fig1, 'Type', 'Axes'), 'Color', 'white');
    plot_filename1 = sprintf('MasaSeca_Comparativa_%dsat_%dtel_Det%d_%s.jpg', N_sat, N_tel, det_id, strrep(tel_name, ' ', '_'));
    saveas(fig1, fullfile(dry_mass_plot_dir, plot_filename1));
    exportgraphics(fig1, fullfile(dry_mass_plot_dir, strrep(plot_filename1, '.jpg', '.pdf')), 'ContentType', 'vector');
    close(fig1);

    % --- GRÁFICO 2: DESGLOSE DE MASA VS ALTURA (teórico + viable) ---
    fig2 = figure('Visible', 'off', 'Position', [100, 100, 900, 600]);
    hold on;

    % Separar teóricos y viables
    if ismember('Es_Viable', data_for_plot.Properties.VariableNames)
        viable_mask = data_for_plot.Es_Viable == true;
    else
        viable_mask = true(size(data_for_plot,1), 1);  % compatibilidad
    end
    theo_mask = ~viable_mask;

    % Curva teórica (gris claro punteado)
    if any(theo_mask)
        plot(data_for_plot.Altura_km(theo_mask), data_for_plot.Masa_Total_Satelite_kg(theo_mask), ...
            'Color', [0.55 0.55 0.55], 'LineStyle', ':', 'LineWidth', 1.5, 'DisplayName', 'Masa total (teorica)');
        plot(data_for_plot.Altura_km(theo_mask), data_for_plot.Masa_Seca_Satelite_kg(theo_mask), ...
            'Color', [0.65 0.65 0.80], 'LineStyle', ':', 'LineWidth', 1.2, 'DisplayName', 'Masa seca (teorica)');
        plot(data_for_plot.Altura_km(theo_mask), data_for_plot.Masa_Combustible_kg(theo_mask), ...
            'Color', [0.60 0.75 0.60], 'LineStyle', ':', 'LineWidth', 1.2, 'DisplayName', 'Masa combustible (teorica)');
    end

    % Curva viable (color sólido)
    if any(viable_mask)
        plot(data_for_plot.Altura_km(viable_mask), data_for_plot.Masa_Total_Satelite_kg(viable_mask), ...
            'r-d', 'LineWidth', 2, 'MarkerSize', 7, 'MarkerFaceColor', 'r', 'DisplayName', 'Masa total (puntos viables)');
        plot(data_for_plot.Altura_km(viable_mask), data_for_plot.Masa_Seca_Satelite_kg(viable_mask), ...
            'b--s', 'LineWidth', 1.5, 'DisplayName', 'Masa seca (viable)');
        plot(data_for_plot.Altura_km(viable_mask), data_for_plot.Masa_Combustible_kg(viable_mask), ...
            'g:^', 'LineWidth', 1.5, 'DisplayName', 'Masa combustible (viable)');

        [min_mass, min_idx] = min(data_for_plot.Masa_Total_Satelite_kg(viable_mask));
        viable_alturas = data_for_plot.Altura_km(viable_mask);
        min_altura = viable_alturas(min_idx);
        plot(min_altura, min_mass, 'ko', 'MarkerSize', 12, 'MarkerFaceColor', 'r', 'HandleVisibility', 'off');
        text_str = sprintf('Minimo viable: %.3f kg @ %d km', min_mass, min_altura);
        text(min_altura + 25, min_mass, text_str, 'VerticalAlignment', 'middle', ...
             'FontWeight', 'bold', 'BackgroundColor', 'white', 'EdgeColor', 'k');
    end

    hold off;
    
    % Configuración del gráfico
    title_str = sprintf('Desglose de Masas vs. Altura\nConfig: %d Sat, %d Tel, Det %d, %s', N_sat, N_tel, det_id, tel_name);
    title(title_str, 'Interpreter', 'latex');
    xlabel('Altura Orbital [km]', 'Interpreter', 'latex');
    ylabel('Masa [kg]', 'Interpreter', 'latex');
    legend('Location', 'northwest', 'Interpreter', 'latex');
    grid on;
    axis tight;
    
    if any(viable_mask)
        viable_alts = data_for_plot.Altura_km(viable_mask);
        xlim([max(400, min(data_for_plot.Altura_km)), max(viable_alts) + 10]);
    end
    
    % Guardar el gráfico (solo PDF)
    plot_filename2 = sprintf('MasaDesglose_%dsat_%dtel_Det%d_%s', N_sat, N_tel, det_id, strrep(tel_name, ' ', '_'));
    fig2.Color = 'white';
    set(findall(fig2, 'Type', 'Axes'), 'Color', 'white');
    exportgraphics(fig2, fullfile(total_mass_plot_dir, [plot_filename2 '.pdf']), 'ContentType', 'vector');
    close(fig2);
end

fprintf('Gráficas de masa guardadas en las carpetas ''%s'' y ''%s''.\n', dry_mass_plot_dir, total_mass_plot_dir);

end
