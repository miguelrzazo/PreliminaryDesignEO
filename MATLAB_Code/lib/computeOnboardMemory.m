function memory = computeOnboardMemory(segmentDataBits, accessMask, capacity_GB, downloadPerSample_GB)
%COMPUTEONBOARDMEMORY Preserve the original backlog and overflow accumulation.
% Time remains sequential; satellites are updated together at each time step.
% Required backlog uses its unbounded previous state, as in the original model.
validateattributes(accessMask,{'logical'},{'2d','nonempty'});
validateattributes(segmentDataBits,{'cell'},{'vector'});
validateattributes(capacity_GB,{'numeric'},{'scalar','real','finite','nonnegative'});
validateattributes(downloadPerSample_GB,{'numeric'},{'scalar','real','finite','nonnegative'});
satelliteCount = size(accessMask,1); sampleCount = size(accessMask,2);
assert(numel(segmentDataBits)==satelliteCount,'eo:MemoryAxes','One data series is required per satellite.');
if sampleCount == 1
    assert(all(cellfun(@isempty,segmentDataBits)),'eo:MemoryAxes','A single sample has no data segments.');
    generated_GB = zeros(satelliteCount,0);
else
    generated_GB = cell2mat(cellfun(@(bits) bits(:).'/8e9,segmentDataBits(:),'UniformOutput',false));
end
assert(isequal(size(generated_GB),[satelliteCount,sampleCount-1]),'eo:MemoryAxes','Data segments do not match time samples.');
memory.bounded_GB = zeros(satelliteCount,sampleCount);
memory.required_GB = zeros(satelliteCount,sampleCount);
memory.overflow_GB = zeros(satelliteCount,1);
for sampleIndex = 2:sampleCount
    requested_GB = max(0,memory.required_GB(:,sampleIndex-1) + generated_GB(:,sampleIndex-1) - ...
        accessMask(:,sampleIndex)*downloadPerSample_GB);
    memory.required_GB(:,sampleIndex) = requested_GB;
    memory.overflow_GB = memory.overflow_GB + max(0,requested_GB-capacity_GB);
    memory.bounded_GB(:,sampleIndex) = min(capacity_GB,requested_GB);
end
end
