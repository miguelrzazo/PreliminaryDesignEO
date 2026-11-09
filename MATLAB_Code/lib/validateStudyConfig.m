function validateStudyConfig(config)
%VALIDATESTUDYCONFIG Reject invalid study inputs before any output is created.
validateattributes(config,{'struct'},{'scalar'});
validateattributes(config.writeOutputs,{'logical'},{'scalar'});
assert(ischar(config.outputRoot) && (isrow(config.outputRoot) || isempty(config.outputRoot)), ...
    'eo:OutputRoot','outputRoot must be a character vector.');
validateattributes(config.groundSampleDistance_m,{'numeric'},{'scalar','real','finite','positive'});
% Existing selection compares band 2 MTF with band 3 SNR. Their roles are fixed.
validateattributes(config.wavelengths_m,{'numeric'},{'vector','numel',3,'real','finite','positive'});
validateDesignAxes(config.altitudes_km,config.pupilDiameters_mm);
validateattributes(config.swaths_km,{'numeric'},{'vector','nonempty','real','finite','positive'});
axes = {config.altitudes_km,config.pupilDiameters_mm,config.swaths_km};
assert(all(cellfun(@(axis) numel(unique(axis)) == numel(axis),axes)), ...
    'eo:DuplicateAxes','Design axes must contain unique coordinates.');
validateattributes(config.snrRequirement,{'numeric'},{'scalar','real','finite','nonnegative'});
validateattributes(config.mtfRequirement,{'numeric'},{'scalar','real','finite','nonnegative'});
validateattributes(config.revisitRequirement_days,{'numeric'},{'scalar','real','finite','positive'});
assert(config.revisitRequirement_days == fix(config.revisitRequirement_days), ...
    'eo:WholeDayHorizon','The current revisit solver requires a whole-day horizon.');
validateattributes(config.maxDetectorCount,{'numeric'},{'scalar','integer','finite','positive'});
validateattributes(config.overlapFraction,{'numeric'},{'scalar','real','finite','>=',0,'<',1});
validateattributes(config.cloudFraction,{'numeric'},{'scalar','real','finite','>=',0,'<',1});
validateattributes(config.localTimeAscendingNode_h,{'numeric'},{'scalar','real','finite','>=',0,'<',24});
validateattributes(config.maintenanceLowerFactor,{'numeric'},{'scalar','real','finite','positive','<',1});
validateattributes(config.maintenanceUpperFactor,{'numeric'},{'scalar','real','finite','>=',1});
validateattributes(config.satelliteConfigurations,{'numeric'}, ...
    {'2d','size',[NaN 2],'nonempty','integer','positive','finite'});
assert(all(ismember(config.satelliteConfigurations(:,2),[1 2])), ...
    'eo:TelescopeCount','The existing model supports one or two telescopes per satellite.');
validateDetectors(config.co2Detectors);
validateDetectors(config.o2Detectors);
assert(iscellstr(config.telescopeNames) && isrow(config.telescopeNames) && ~isempty(config.telescopeNames), ...
    'eo:TelescopeNames','telescopeNames must be a nonempty row cell array of character vectors.');
telescopeCount = numel(config.telescopeNames);
validateattributes(config.telescopeAlignmentMtf,{'numeric'},{'vector','numel',telescopeCount,'real','finite','>=',0,'<=',1});
validateattributes(config.telescopeTransmittance,{'numeric'},{'vector','numel',telescopeCount,'real','finite','>=',0,'<=',1});
validateattributes(config.telescopeFov_deg,{'numeric'},{'vector','numel',telescopeCount,'real','finite','positive','<',90});
validateattributes(config.obstructionRatios,{'numeric'},{'vector','numel',telescopeCount,'real','finite','>=',0,'<',1});
end

function validateDetectors(detectors)
% The existing selection/report schema has three paired detectors per band.
assert(iscellstr(detectors.names) && numel(detectors.names) == 3, ...
    'eo:DetectorCount','The existing downstream model requires three detector entries per band.');
validateattributes(detectors.quantumEfficiency,{'numeric'},{'vector','numel',3,'real','finite','>=',0,'<=',1});
validateattributes(detectors.detectorMtf,{'numeric'},{'vector','numel',3,'real','finite','>=',0,'<=',1});
validateattributes(detectors.pixelCount,{'numeric'},{'vector','numel',3,'integer','finite','positive'});
validateattributes(detectors.pixelSize_m,{'numeric'},{'vector','numel',3,'real','finite','positive'});
end
