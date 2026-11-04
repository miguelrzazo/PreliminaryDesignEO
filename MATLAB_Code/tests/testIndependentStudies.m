function tests=testIndependentStudies
tests=functiontests(localfunctions);
end
function setupOnce(testCase)
root=fileparts(fileparts(mfilename('fullpath')));
testCase.TestData.oldPath=path;
addpath(root,fullfile(root,'lib'),fullfile(root,'tests','reference'));
end
function teardownOnce(testCase)
path(testCase.TestData.oldPath);
end
function testLifetimeAnchorsAndOriginalSegments(testCase)
altitudes_km=[300 400 450 500 600 700 800 900 1200 2000 2100];
expected=arrayfun(@referenceLifetime,altitudes_km);
actual=computeLifetimeApproximation(altitudes_km,[400 500 700 900 2000],[365 3650 36500 365000 3650000]);
verifyEqual(testCase,isnan(actual),isnan(expected)); mask=isfinite(expected);
verifyEqual(testCase,actual(mask),expected(mask),'RelTol',1e-12);
end
function testWeatherPreservesSeededSequenceAndRng(testCase)
config=struct('dayCount',365,'randomSeed',42,'writeOutputs',false, ...
    'conditionNames',{{'Cond. favorables','Lluvia moderada','Lluvia intensa','Nubes densas'}}, ...
    'rateRanges_Mbps',[130 150;120 140;90 120;125 145], ...
    'seasonProbabilities_percent',[76 12 3 9;52 36 7 5;64 22 5 9;82 1 0 17]);
previousRng=rng;
expected=referenceWeather(); actual=runWeatherLinkStudy(config);
currentRng=rng;
verifyEqual(testCase,currentRng,previousRng);
verifyEqual(testCase,actual.dailyConditions,expected.dailyConditions);
verifyEqual(testCase,actual.dailyRates_Mbps,expected.dailyRates_Mbps,'AbsTol',1e-12);
verifyEqual(testCase,actual.weeklyRates_Mbps,expected.weeklyRates_Mbps,'AbsTol',1e-12);
end
function testAcquisitionEnvelopeMatchesScalarPasses(testCase)
eccentricities=[0;0.001;0.002]; angles_rad=[pi/2;pi/3;pi/4]; latitudes_deg=linspace(24,26,100);
[minimum,maximum]=computeAcquisitionHeightEnvelope(eccentricities,angles_rad,latitudes_deg,97.9,7000e3,6378.137e3);
for caseIndex=1:3
    [expectedMinimum,expectedMaximum]=referenceHeightEnvelope(eccentricities(caseIndex),angles_rad(caseIndex), ...
        latitudes_deg,97.9,7000e3,6378.137e3);
    verifyEqual(testCase,minimum(caseIndex),expectedMinimum,'AbsTol',1e-10);
    verifyEqual(testCase,maximum(caseIndex),expectedMaximum,'AbsTol',1e-10);
end
end
function testRgtReference(testCase)
config=struct('earthGravitationalParameter_m3_s2',3.986004418e14,'earthRadius_m',6378.137e3, ...
    'earthJ2',1.08263e-3,'earthJ3',-2.5327e-6,'earthRotation_rad_s',7.292115373194e-5, ...
    'sunNodalRate_rad_s',1.991e-7,'altitudes_km',200:2:1000,'rootBounds_m',[200e3 1000e3], ...
    'referenceAltitude_km',630,'maxRepeatDays',30,'rootTolerance',1e-10,'maxIterations',100,'writeOutputs',false);
expected=referenceRgt(); actual=runRgtStudy(config);
verifyEqual(testCase,actual.solutions,expected.solutions,'AbsTol',1e-7,'RelTol',1e-10);
verifyEqual(testCase,actual.residualDrift_km_day,expected.residualDrift_km_day,'AbsTol',1e-8);
verifyEqual(testCase,actual.inclinations_deg,expected.inclinations_deg,'AbsTol',1e-12);
end
function testMemoryRetainsSequentialBacklogModel(testCase)
bits={8e9*[1 2 0 1];8e9*[2 0 1 2]}; mask=logical([0 0 1 0 1;0 1 0 1 0]);
expected=zeros(2,5); bounded=expected; overflow=zeros(2,1);
for sampleIndex=2:5
    for satelliteIndex=1:2
        requested=max(0,expected(satelliteIndex,sampleIndex-1)+bits{satelliteIndex}(sampleIndex-1)/8e9-mask(satelliteIndex,sampleIndex));
        expected(satelliteIndex,sampleIndex)=requested;
        bounded(satelliteIndex,sampleIndex)=min(2,requested);
        overflow(satelliteIndex)=overflow(satelliteIndex)+max(0,requested-2);
    end
end
actual=computeOnboardMemory(bits,mask,2,1);
verifyEqual(testCase,actual.required_GB,expected); verifyEqual(testCase,actual.bounded_GB,bounded);
verifyEqual(testCase,actual.overflow_GB,overflow);
end
function testSingleMemorySample(testCase)
actual=computeOnboardMemory({[];[]},false(2,1),2,1);
verifyEqual(testCase,actual.required_GB,zeros(2,1));
verifyEqual(testCase,actual.overflow_GB,zeros(2,1));
end
function testSensitivityLeavesInputUnchanged(testCase)
folder=tempname; mkdir(folder);
cleanup=onCleanup(@() rmdir(folder,'s'));
inputFile=fullfile(folder,'samples.csv');
source=table([100;80],[90;28],[520;520],[2;2],[2;2],zeros(2,1),zeros(2,1),zeros(2,1), ...
    'VariableNames',{'GSD_m','Diametro_Pupila_mm','Altura_km','Num_Telescopios', ...
    'Num_Satelites','Masa_Seca_Satelite_kg','Masa_Combustible_kg','Masa_Total_Constelacion_kg'});
writetable(source,inputFile); originalText=fileread(inputFile);
config=struct('inputFile',inputFile,'outputDirectory',fullfile(folder,'output'), ...
    'figureDirectory',fullfile(folder,'figures'),'displaySatelliteCount',2, ...
    'updateInputSummary',false,'writeOutputs',false);
actual=runGsdSensitivityStudy(config);
verifyEqual(testCase,fileread(inputFile),originalText);
verifyFalse(testCase,isfolder(config.outputDirectory));
verifyEqual(testCase,actual.summary.GSD_m,[80;100]);
for rowIndex=1:height(actual.summary)
    [expectedDryMass,~,~,~,~,~,~,~,~,~,expectedArea]= ...
        referenceDryMass(actual.summary.Diametro_Pupila_mm(rowIndex),2);
    [expectedTotalMass,expectedFuelMass]=referenceTotalMass(520,expectedDryMass,expectedArea);
    verifyEqual(testCase,actual.summary.Masa_Seca_Satelite_kg(rowIndex),expectedDryMass,'RelTol',1e-10);
    verifyEqual(testCase,actual.summary.Masa_Combustible_kg(rowIndex),expectedFuelMass,'RelTol',1e-10);
    verifyEqual(testCase,actual.summary.Masa_Total_Constelacion_kg(rowIndex),2*expectedTotalMass,'RelTol',1e-10);
end
end
function testPhasingOriginalSampling(testCase)
config=struct('earthGravitationalParameter_m3_s2',3.986004418e14,'earthRadius_m',6378.137e3, ...
    'finalAltitude_m',630e3,'phaseApogeeAltitude_m',830e3,'timeStep_s',120,'writeOutputs',false);
actual=runPhasingStudy(config);
finalRadius=config.earthRadius_m+config.finalAltitude_m;
apogeeRadius=config.earthRadius_m+config.phaseApogeeAltitude_m;
phaseAxis=(finalRadius+apogeeRadius)/2;
finalRate=2*pi/(2*pi*sqrt(finalRadius^3/config.earthGravitationalParameter_m3_s2));
phaseRate=2*pi/(2*pi*sqrt(phaseAxis^3/config.earthGravitationalParameter_m3_s2));
time_s=0; accumulated=0; lastTime_s=0;
while accumulated<pi
    lastTime_s=time_s; accumulated=(finalRate-phaseRate)*time_s; time_s=time_s+120;
end
verifyEqual(testCase,actual.time_s(end),lastTime_s);
verifyEqual(testCase,actual.displayedDeploymentTime_s,time_s);
end
function testSpectrumComputeOnly(testCase)
config=struct('sBand_Hz',[2e9 4e9],'xBand_Hz',[8e9 12e9], ...
    'detectorWavelengths_m',[0.76e-6 2.3e-6],'writeOutputs',false);
actual=runSpectrumStudy(config);
verifyEqual(testCase,actual.detectorBand_Hz,[3e8/2.3e-6 3e8/0.76e-6]);
end
