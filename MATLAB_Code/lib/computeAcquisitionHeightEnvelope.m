function [minimumHeights_km,maximumHeights_km,perArcVariation_km,heights_km] = computeAcquisitionHeightEnvelope( ...
    eccentricities,perigeeAngles_rad,arcLatitudes_deg,inclination_deg,semimajorAxis_m,earthRadius_m)
%COMPUTEACQUISITIONHEIGHTENVELOPE Evaluate both passes for every time sample.
% Time-dependent elements are columns; the two latitude passes form one row.
% Preserve the source ellipsoid and true-anomaly construction.
latitudes_deg=[arcLatitudes_deg(:).',fliplr(arcLatitudes_deg(:).')];
geocentricLatitudes_deg=atand((1-1/298.257223563)^2*tand(arcLatitudes_deg(:).'));
ascendingAngles_rad=asin(sind(geocentricLatitudes_deg)/sind(inclination_deg));
descendingAngles_rad=pi-asin(sind(fliplr(geocentricLatitudes_deg))/sind(inclination_deg));
anomalies_rad=[ascendingAngles_rad,descendingAngles_rad]-perigeeAngles_rad(:);
radii_m=semimajorAxis_m*(1-eccentricities(:).^2) ./ (1+eccentricities(:).*cos(anomalies_rad));
polarRadius_m=earthRadius_m*(1-1/298.257223563);
ellipsoidRadii_m=sqrt(((earthRadius_m^2.*cosd(latitudes_deg)).^2 + ...
    (polarRadius_m^2.*sind(latitudes_deg)).^2) ./ ...
    ((earthRadius_m.*cosd(latitudes_deg)).^2+(polarRadius_m.*sind(latitudes_deg)).^2));
heights_km=(radii_m-ellipsoidRadii_m)/1e3;
minimumHeights_km=min(heights_km,[],2); maximumHeights_km=max(heights_km,[],2);
arcSize = numel(arcLatitudes_deg);
perArcVariation_km = [max(heights_km(:,1:arcSize),[],2)-min(heights_km(:,1:arcSize),[],2), ...
    max(heights_km(:,arcSize+1:end),[],2)-min(heights_km(:,arcSize+1:end),[],2)];
end
