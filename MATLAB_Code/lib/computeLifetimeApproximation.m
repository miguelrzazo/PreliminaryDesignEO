function lifetime_days = computeLifetimeApproximation(altitudes_km, anchorAltitudes_km, anchorLifetimes_days)
%COMPUTELIFETIMEAPPROXIMATION Preserve the original empirical log-log segments.
% This is an approximate table fit, not an orbital propagator. Below the first
% anchor the historical result is zero; above the final anchor it remains NaN.
validateattributes(altitudes_km,{'numeric'},{'real','finite'});
validateattributes(anchorAltitudes_km,{'numeric'},{'vector','real','positive','finite'});
validateattributes(anchorLifetimes_days,{'numeric'},{'vector','positive','finite','numel',numel(anchorAltitudes_km)});
assert(all(diff(anchorAltitudes_km)>0),'eo:LifetimeAnchors','Lifetime anchor altitudes must be strictly increasing.');
lifetime_days = zeros(size(altitudes_km));
mask = altitudes_km >= anchorAltitudes_km(1);
lifetime_days(mask) = exp(interp1(log(anchorAltitudes_km),log(anchorLifetimes_days),log(altitudes_km(mask))));
end
