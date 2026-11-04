function [h_min, h_max] = referenceHeightEnvelope(eccentricity, omega, ...
        lat_arc, inclination_deg, a, R_E)
    h_passes = zeros(2, numel(lat_arc));
    for pass_idx = 1:2
        if pass_idx == 1
            lat_pass = lat_arc;
            u = deg2rad(asind(sind(atand((1-1/298.257223563)^2*tand(lat_pass))) / sind(inclination_deg)));
        else
            lat_pass = fliplr(lat_arc);
            u = pi - deg2rad(asind(sind(atand((1-1/298.257223563)^2*tand(lat_pass))) / sind(inclination_deg)));
        end
        nu = u - omega;
        r_orbit = a * (1 - eccentricity^2) ./ ...
            (1 + eccentricity * cos(nu));
        R_ellipsoid = ellipsoidRadius(R_E, lat_pass);
        h_passes(pass_idx, :) = (r_orbit - R_ellipsoid) / 1e3;
    end
    h_min = min(h_passes, [], 'all');
    h_max = max(h_passes, [], 'all');
end

function radius = ellipsoidRadius(equatorial_radius, latitude)
    flattening = 1 / 298.257223563;
    polar_radius = equatorial_radius * (1 - flattening);
    radius = sqrt(((equatorial_radius^2 .* cosd(latitude)).^2 + ...
        (polar_radius^2 .* sind(latitude)).^2) ./ ...
        ((equatorial_radius .* cosd(latitude)).^2 + ...
        (polar_radius .* sind(latitude)).^2));
end

