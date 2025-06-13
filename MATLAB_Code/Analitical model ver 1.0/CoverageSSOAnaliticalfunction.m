function CoverageSSOAnaliticalfunction(GSD, alturas_orbitales, swaths_km, Cov_Requirement, current_Npix, current_fov_limit, N_sat, N_telescopes, detector_type, telescope_index, telescope_name, solapamiento, cobertura_nubes, max_detectores)
 
    
    % Constants based on sun-synchronous orbit characteristics
    earth_radius = 6371; % km
    continental_US_width = 4100; % km (width from Atlantic to Pacific coast)
    earth_circumference = 2 * pi * earth_radius; % km
    sidereal_day = 86164; % seconds (sidereal day)
    
    % Use input factors
    daylight_factor = 0.5; % Only daytime passes are useful
    
    % Configuration name based on current parameters
    config_name = sprintf('%dSat_%dTel_%s_Det%d', N_sat, N_telescopes, strrep(telescope_name, ' ', ''), detector_type);
    
    % Check if configuration respects max_detectores limit
    if detector_type > max_detectores
        fprintf('Detector %d exceeds max detectors limit (%d)\n', detector_type, max_detectores);
        return;
    end
    
    % Calculate swath limits for each number of detectors (1 to max_detectores)
    detector_swath_limits = [];
    for n_det = 1:max_detectores
        swath_limit = n_det * GSD * current_Npix / 1000; % Convert to km (GSD in m, result in km)
        detector_swath_limits(end+1) = swath_limit;
    end
    
    % Preallocate results
    coverage_days = zeros(length(alturas_orbitales), length(swaths_km));
    
    for h = 1:length(alturas_orbitales)
        for s = 1:length(swaths_km)
            height = alturas_orbitales(h);
            swath = swaths_km(s);
            
            % Check FOV constraint: swath must be achievable with given FOV limit
            % For nadir pointing: swath = 2 * height * tan(FOV/2)
            max_achievable_swath = 2 * height * tand(current_fov_limit/2);
            
            % Check swath limit constraint for current detector type
            max_detector_swath = detector_type * GSD * current_Npix / 1000; % km
            
            if swath > max_achievable_swath || swath > max_detector_swath
                % Constraints exceeded, mark as NaN
                coverage_days(h,s) = NaN;
                continue;
            end
            
            % Calculate orbital parameters for sun-synchronous orbit
            orbital_velocity = sqrt(398600.4418 / (earth_radius + height)); % km/s
            orbital_period = 2 * pi * (earth_radius + height) / orbital_velocity; % seconds
            orbits_per_day = sidereal_day / orbital_period;
            
            % Calculate effective swath considering overlap and configuration
            effective_swath = swath * N_telescopes * (1 - solapamiento);
            
            % Calculate ground track spacing due to Earth rotation
            ground_track_spacing = earth_circumference / orbits_per_day;
            
            % Number of passes needed to cover US width
            passes_needed = ceil(continental_US_width / effective_swath);
            
            % Time between useful passes (considering only daylight passes)
            time_between_passes = orbital_period / daylight_factor; % seconds
            
            % Total time to cover US width (for single satellite)
            theoretical_coverage_time = passes_needed * time_between_passes; % seconds
            theoretical_coverage_days = theoretical_coverage_time / sidereal_day;
            
            % Apply cloud cover factor (increases time needed)
            base_days = theoretical_coverage_days / (1 - cobertura_nubes);
            
            % Apply satellite constellation factor
            constellation_days = base_days / N_sat;
            
            % Efficiency factor based on configuration complexity
            if N_telescopes > 1 || detector_type > 1
                efficiency_factor = 0.9; % 10% efficiency loss for complex configurations
                final_days = constellation_days / efficiency_factor;
            else
                final_days = constellation_days;
            end
            
            % Apply coverage requirement threshold
            if final_days > Cov_Requirement
                coverage_days(h,s) = NaN;
            else
                coverage_days(h,s) = final_days;
            end
        end
    end
    
    % Save results to txt and csv files
    txt_filename = fullfile('coverage', sprintf('coverage_%s.txt', config_name));
    csv_filename = fullfile('coverage', sprintf('coverage_%s.csv', config_name));
    
    % Create header for txt file
    header_txt = sprintf('Coverage analysis for %s\n', config_name);
    header_txt = [header_txt sprintf('Heights (km): %s\n', mat2str(alturas_orbitales))];
    header_txt = [header_txt sprintf('Swaths (km): %s\n', mat2str(swaths_km))];
    header_txt = [header_txt sprintf('Coverage requirement: %.0f days\n', Cov_Requirement)];
    header_txt = [header_txt sprintf('Continental US width: %.0f km\n', continental_US_width)];
    header_txt = [header_txt sprintf('GSD: %.2f m\n', GSD)];
    header_txt = [header_txt sprintf('Pixels per detector: %d\n', current_Npix)];
    header_txt = [header_txt sprintf('Cloud coverage factor: %.2f\n', cobertura_nubes)];
    header_txt = [header_txt sprintf('Swath overlap factor: %.2f\n', solapamiento)];
    header_txt = [header_txt sprintf('FOV limit: %.1f degrees\n', current_fov_limit)];
    header_txt = [header_txt sprintf('Max detectors: %d\n', max_detectores)];
    header_txt = [header_txt sprintf('Number of satellites: %d\n', N_sat)];
    header_txt = [header_txt sprintf('Number of telescopes: %d\n', N_telescopes)];
    header_txt = [header_txt sprintf('Detector type: %d\n', detector_type)];
    header_txt = [header_txt sprintf('Max swath for this detector: %.2f km\n\n', detector_type * GSD * current_Npix / 1000)];
    
    % Write txt file with header
    fid = fopen(txt_filename, 'w');
    fprintf(fid, '%s', header_txt);
    fclose(fid);
    
    % Append data to txt file
    writematrix(coverage_days, txt_filename, 'Delimiter', 'tab', 'WriteMode', 'append');
    
    % Save csv file
    writematrix(coverage_days, csv_filename);
    
    % Generate single heatmap with detector limit lines
    fig = figure('Visible', 'off');
    
    % Plot the matrix as imagesc to allow overlay of lines
    imagesc(swaths_km, alturas_orbitales, coverage_days);
    colormap(parula);
    
    % Set NaN values to white
    current_colormap = colormap;
    colormap([1 1 1; current_colormap]); % Add white at beginning
    
    % Add vertical dashed lines for detector limits (1 to max_detectores)
    hold on;
    for n_det = 1:max_detectores
        swath_limit = detector_swath_limits(n_det);
        if swath_limit <= max(swaths_km) && swath_limit >= min(swaths_km)
            xline(swath_limit, '--k', 'LineWidth', 2, 'Alpha', 0.8);
            text(swath_limit, max(alturas_orbitales)*0.95, sprintf('%d Det', n_det), ...
                 'HorizontalAlignment', 'center', 'FontSize', 10, 'Color', 'black', 'FontWeight', 'bold');
        end
    end
    hold off;
    
    % Formatting
    title_text = sprintf('Coverage Days - %s\nBlanco: >%d dias', strrep(config_name, '_', ' '), Cov_Requirement);
    title(title_text, 'Interpreter', 'latex', 'FontSize', 14);
    xlabel('Swath (km)', 'Interpreter', 'latex', 'FontSize', 12);
    ylabel('Height (km)', 'Interpreter', 'latex', 'FontSize', 12);
    
    % Add colorbar
    c = colorbar;
    ylabel(c, 'Coverage Days', 'Interpreter', 'latex', 'FontSize', 11);
    
    % Set axis properties
    set(gca, 'YDir', 'normal');
    axis tight;
    
    % Save heatmap
    png_filename = fullfile('coverage', sprintf('heatmap_%s.png', config_name));
    print(fig, png_filename, '-dpng', '-r300');
    
    % Close figure
    close(fig);
    
    % Display summary information
    fprintf('Coverage analysis completed for: %s\n', strrep(config_name, '_', ' '));
    fprintf('Height range: %.0f-%.0f km\n', min(alturas_orbitales), max(alturas_orbitales));
    fprintf('Swath range: %.0f-%.0f km\n', min(swaths_km), max(swaths_km));
    fprintf('Coverage requirement threshold: %.0f days\n', Cov_Requirement);
    fprintf('Detector type: %d, Max swath: %.2f km\n', detector_type, detector_type * GSD * current_Npix / 1000);
    fprintf('GSD: %.2f m, Npix: %d, FOV limit: %.1f°\n', GSD, current_Npix, current_fov_limit);
    fprintf('Results saved in ''coverage'' folder\n');
end
