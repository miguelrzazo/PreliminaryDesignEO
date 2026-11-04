function cases = analysisCases(configurationCount, detectorCount, telescopeCount, altitudeCount)
%ANALYSISCASES Enumerate metadata cases; altitude is the fastest-varying index.
if nargin < 4, altitudeCount = 1; end
[altitudeIndices,telescopeIndices,detectorIndices,configurationIndices] = ...
    ndgrid(1:altitudeCount,1:telescopeCount,1:detectorCount,1:configurationCount);
cases = [configurationIndices(:),detectorIndices(:),telescopeIndices(:),altitudeIndices(:)];
end
