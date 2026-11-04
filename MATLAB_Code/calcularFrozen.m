function summary = calcularFrozen(h_selected, output_dir)
% Compare first-order J2/J3 mean orbits over the CONUS latitude bounds.
% Each direction uses a continuous 25--49 deg N arc. Ignoring longitude
% clipping gives a latitude-envelope bound, not an exact geographic pass.
if nargin < 1, h_selected = 630; end
if nargin < 2, output_dir = fullfile(pwd, 'review-results', 'Frozen'); end
if ~exist(output_dir, 'dir'), mkdir(output_dir); end
validateattributes(h_selected, {'numeric'}, {'scalar','finite','positive'});

mu = 3.986004418e14;
R_E = 6378.137e3;
J2 = 1.08263e-3;
J3 = -2.5327e-6;
w_sol = 1.991e-7;
a = R_E + h_selected * 1000;
n = sqrt(mu/a^3);
inc = acos(-2*a^(7/2)*w_sol/(3*J2*R_E^2*sqrt(mu)));
e_f = -J3/(2*J2)*(R_E/a)*sin(inc);
t_days = (0:8*365.25).';
lat_geodetic = linspace(25, 49, 241);
f = 1/298.257223563;
lat_geocentric = atand((1-f)^2*tand(lat_geodetic));
u_asc = asin(sind(lat_geocentric)/sin(inc));
u = [u_asc; pi-u_asc];
b = R_E*(1-f);
surface_radius = sqrt(((R_E^2*cosd(lat_geodetic)).^2 + ...
    (b^2*sind(lat_geodetic)).^2)./ ...
    ((R_E*cosd(lat_geodetic)).^2 + (b*sind(lat_geodetic)).^2));

% Use the same first-order fixed point for both initial conditions.
p = a*(1-e_f^2);
alpha_rate = -3*n*J2*(R_E/p)^2*(5/4*sin(inc)^2-1);
alpha = alpha_rate*t_days*86400;
case_names = {'frozen', 'circular'};
summary_rows = zeros(2, 9);
for c = 1:2
    if c == 1
        ex = zeros(size(alpha)); ey = e_f*ones(size(alpha));
    else
        ex = e_f*sin(alpha); ey = e_f*(1-cos(alpha));
    end
    ecc = hypot(ex, ey);
    omega = atan2(ey, ex);
    rows = zeros(numel(t_days), 9);
    for k = 1:numel(t_days)
        nu = u-omega(k);
        radius = a*(1-ecc(k)^2)./(1+ecc(k)*cos(nu));
        heights = (radius-surface_radius)/1000;
        per_pass = max(heights, [], 2)-min(heights, [], 2);
        envelope = max(heights(:))-min(heights(:));
        rows(k,:) = [t_days(k), ecc(k), rad2deg(omega(k)), ...
            min(heights(:)), max(heights(:)), per_pass.', envelope, ...
            max(abs(heights(:)-h_selected))];
        assert(all(per_pass <= envelope+1e-10), 'Pass exceeds joint envelope.');
    end
    assert(all(isfinite(rows(:))), 'Non-finite height result.');
    summary_rows(c,:) = [min(ecc), max(ecc), min(rows(:,4)), ...
        max(rows(:,5)), max(rows(:,6)), max(rows(:,7)), ...
        max(rows(:,8)), max(rows(:,9)), h_selected];
    names = {'MissionDay','Eccentricity','ArgumentOfPerigee_deg', ...
        'MinHeight_km','MaxHeight_km','AscendingArcVariation_km', ...
        'DescendingArcVariation_km','CombinedEnvelope_km', ...
        'MaxDeviationFromNominal_km'};
    writetable(array2table(rows, 'VariableNames', names), ...
        fullfile(output_dir, [case_names{c}, '_evolution.csv']));
    fprintf('FROZEN_RESULT %s: emax=%.9g, asc=%.6f km, desc=%.6f km, envelope=%.6f km\n', ...
        case_names{c}, summary_rows(c,2), summary_rows(c,5:7));
end
names = {'MinEccentricity','MaxEccentricity','MinHeight_km','MaxHeight_km', ...
    'MaxAscendingArcVariation_km','MaxDescendingArcVariation_km', ...
    'MaxCombinedEnvelope_km','MaxDeviationFromNominal_km','NominalAltitude_km'};
summary = array2table(summary_rows, 'VariableNames', names);
summary.Case = case_names.';
writetable(summary, fullfile(output_dir, 'frozen_comparison.csv'));
end
