function tests = testStudyWorkflow
%TESTSTUDYWORKFLOW Exercise explicit inputs and workflow resource boundaries.
tests = functiontests(localfunctions);
end

function setupOnce(testCase)
root = fileparts(fileparts(mfilename('fullpath')));
testCase.TestData.oldPath = path;
addpath(root,fullfile(root,'lib'));
end

function teardownOnce(testCase)
path(testCase.TestData.oldPath);
end

function testComputeOnlyPreservesInputsAndCallerState(testCase)
config = fixture();
originalConfig = config;
originalDirectory = pwd;
originalPath = path;
folder = tempname;
config.outputRoot = folder;
results = runAnalysis(config,@(varargin) 1);
verifyEqual(testCase,results.config,config);
verifyEqual(testCase,originalConfig.outputRoot,'');
verifyEqual(testCase,pwd,originalDirectory);
verifyEqual(testCase,path,originalPath);
verifyFalse(testCase,isfolder(folder));
verifyEqual(testCase,size(results.optics.mtf{1,1,1}),[2 2]);
verifyEqual(testCase,results.axes.altitudes_km,config.altitudes_km(:));
verifyTrue(testCase,results.coverage.executionSucceeded);
verifyEmpty(testCase,results.configuration.rows);
end

function testDetectorAndMissionInputsReachModels(testCase)
config = fixture();
config.groundSampleDistance_m = 100;
config.wavelengths_m = [1.5e-6 2.1e-6 0.8e-6];
config.co2Detectors.pixelSize_m(1) = 20e-6;
config.co2Detectors.quantumEfficiency(1) = 0.6;
config.co2Detectors.detectorMtf(1) = 0.4;
config.co2Detectors.pixelCount(1) = 1000;
config.satelliteConfigurations = [2 2];
config.mtfRequirement = 0;
config.snrRequirement = 0;
config.swaths_km = 1e6; % Reject all coverage, avoiding calls to the real toolkit.
results = runAnalysis(config);
verifyEqual(testCase,results.config,config);
expectedSnr = computeSnrGrid(1.5e-6,20e-6,0.6,0.8,100,0, ...
    config.altitudes_km,config.pupilDiameters_mm,0);
expectedMtf = computeMtfGrid(1.5e-6,20e-6,0.4,0.9,100,0, ...
    config.altitudes_km,config.pupilDiameters_mm,0);
verifyEqual(testCase,results.optics.snr{1,1,1},expectedSnr,'RelTol',1e-12);
verifyEqual(testCase,results.optics.mtf{1,1,1},expectedMtf,'AbsTol',1e-12);
verifyTrue(testCase,all(isnan(results.selection{1}.viableDiameters_mm(:))));
end

function testCoverageFailuresRemainVisible(testCase)
config = fixture();
results = runAnalysis(config,@failingSolver);
verifyFalse(testCase,results.coverage.executionSucceeded);
verifyEqual(testCase,results.coverage.diagnostics{1,1,1}.solverErrorIdentifiers{1}, ...
    'test:WorkflowSolverFailure');
end

function revisit_days = failingSolver(varargin)
error('test:WorkflowSolverFailure','Deliberate solver failure');
revisit_days = NaN;
end

function testValidationPrecedesOutputCreation(testCase)
config = fixture();
config.writeOutputs = true;
config.outputRoot = tempname;
config.revisitRequirement_days = 7.5;
originalPath = path;
originalDirectory = pwd;
verifyError(testCase,@() runAnalysis(config),'eo:WholeDayHorizon');
verifyFalse(testCase,isfolder(config.outputRoot));
verifyEqual(testCase,path,originalPath);
verifyEqual(testCase,pwd,originalDirectory);
end

function testUnsupportedDetectorCountFailsClearly(testCase)
config = fixture();
config.co2Detectors.names = {'One detector'};
verifyError(testCase,@() runAnalysis(config),'eo:DetectorCount');
end

function config = fixture()
% Small explicit input set; do not execute master or depend on environment overrides.
config = struct('writeOutputs',false,'outputRoot','', ...
    'groundSampleDistance_m',80,'wavelengths_m',[1.61e-6 2.01e-6 0.76e-6], ...
    'altitudes_km',[400 520],'pupilDiameters_mm',[28 90],'swaths_km',[20 500], ...
    'snrRequirement',1e9,'mtfRequirement',1,'revisitRequirement_days',7, ...
    'maxDetectorCount',3,'overlapFraction',0.05,'cloudFraction',1/6, ...
    'localTimeAscendingNode_h',6,'maintenanceLowerFactor',0.98, ...
    'maintenanceUpperFactor',1.02,'satelliteConfigurations',[1 1], ...
    'telescopeNames',{{'Example'}},'telescopeAlignmentMtf',0.9, ...
    'telescopeTransmittance',0.8,'telescopeFov_deg',10,'obstructionRatios',0);
config.co2Detectors = struct('names',{{'CAPYORK','H2RG','SATURN VISIR'}}, ...
    'quantumEfficiency',[0.85 0.7 0.6],'detectorMtf',[0.45 0.5 0.48], ...
    'pixelCount',[1200 2048 1000],'pixelSize_m',[15 18 30]*1e-6);
config.o2Detectors = struct('names',{{'CCD98-50','H2RG','SATURN VISIR'}}, ...
    'quantumEfficiency',[0.65 0.7 0.6],'detectorMtf',[0.5 0.5 0.45], ...
    'pixelCount',[6000 2048 1000],'pixelSize_m',[13 18 30]*1e-6);
end
