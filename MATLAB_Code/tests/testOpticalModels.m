function tests = testOpticalModels
%TESTOPTICALMODELS Regression agreement does not establish scientific validity.
tests = functiontests(localfunctions);
end

function setupOnce(testCase)
root = fileparts(fileparts(mfilename('fullpath')));
testCase.TestData.oldPath = path;
addpath(root, fullfile(root,'lib'), fullfile(root,'tests','reference'));
end

function teardownOnce(testCase)
path(testCase.TestData.oldPath);
end

function testOpticalScalarAgreement(testCase)
% Non-square, singleton and transposed input axes expose orientation errors.
grids = {[200 520 1000], [20 28 90 200]; 520, [20 200]; [200;520], 28};
for g = 1:size(grids,1)
    altitudes_km = grids{g,1}; pupilDiameters_mm = grids{g,2};
    for obstructionRatio = [0 0.2]
        for wavelength_m = [1.61e-6 2.01e-6 0.76e-6]
            for threshold = [0 400]
                expected = referenceSNRfunction(wavelength_m,18e-6,0.7,0.8,80,obstructionRatio,altitudes_km,pupilDiameters_mm,'','',1,threshold);
                actual = computeSnrGrid(wavelength_m,18e-6,0.7,0.8,80,obstructionRatio,altitudes_km,pupilDiameters_mm,threshold);
                assertGrid(testCase,actual,expected);
            end
            for threshold = [0 0.25]
                expected = referenceMTFfunction(wavelength_m,18e-6,0.5,0.9,80,obstructionRatio,altitudes_km,pupilDiameters_mm,'','',1,threshold);
                actual = computeMtfGrid(wavelength_m,18e-6,0.5,0.9,80,obstructionRatio,altitudes_km,pupilDiameters_mm,threshold);
                assertGrid(testCase,actual,expected);
            end
        end
    end
end
end

function testMtfPiecewiseBoundaries(testCase)
% X = altitude_m * lambda / (2*GSD*diameter_m).
for obstructionRatio = [0 0.2]
    normalizedFrequencies = unique([0.01 obstructionRatio (1-obstructionRatio)/2 (1+obstructionRatio)/2 1 1.1]);
    normalizedFrequencies = normalizedFrequencies(normalizedFrequencies > 0);
    pupilDiameters_mm = 520e3 * 2.01e-6 ./ (2*80*normalizedFrequencies) * 1000;
    expected = referenceMTFfunction(2.01e-6,18e-6,0.5,0.9,80,obstructionRatio,520,pupilDiameters_mm,'','',1,0);
    actual = computeMtfGrid(2.01e-6,18e-6,0.5,0.9,80,obstructionRatio,520,pupilDiameters_mm,0);
    assertGrid(testCase,actual,expected);
end
end

function testOpticalExportRoundTrip(testCase)
% Check export metadata and CSV axes.
folder = tempname;
mkdir(folder);
cleanup = onCleanup(@() rmdir(folder,'s'));
snrPrefix = fullfile(folder,'snr');
mtfPrefix = fullfile(folder,'mtf');
altitudes_km = [400 520]; pupilDiameters_mm = [20 28 200];
snr = SNRfunction(2.01e-6,18e-6,0.7,0.8,80,0.2,altitudes_km,pupilDiameters_mm,snrPrefix,'Example',1,400,true);
mtf = MTFfunction(2.01e-6,18e-6,0.5,0.9,80,0.2,altitudes_km,pupilDiameters_mm,mtfPrefix,'Example',1,0.25,true);
assertGrid(testCase,table2array(readtable([snrPrefix '_resultados.csv'],'ReadRowNames',true)),snr);
assertGrid(testCase,table2array(readtable([mtfPrefix '_resultados.csv'],'ReadRowNames',true)),mtf);
verifyTrue(testCase,contains(fileread([snrPrefix '_estadisticas.txt']),'Ancho de banda: 20.0 nm'));
end

function testAxisOrientationDoesNotChangeResults(testCase)
altitudes_km = [200 520]; pupilDiameters_mm = [20 28 90];
snr = computeSnrGrid(2.01e-6,18e-6,0.7,0.8,80,0.2,altitudes_km,pupilDiameters_mm,0);
mtf = computeMtfGrid(2.01e-6,18e-6,0.5,0.9,80,0.2,altitudes_km,pupilDiameters_mm,0);
assertGrid(testCase,snr,computeSnrGrid(2.01e-6,18e-6,0.7,0.8,80,0.2,altitudes_km.',pupilDiameters_mm.',0));
assertGrid(testCase,mtf,computeMtfGrid(2.01e-6,18e-6,0.5,0.9,80,0.2,altitudes_km.',pupilDiameters_mm.',0));
end

function testRequirementBoundaryIsAccepted(testCase)
snr = computeSnrGrid(2.01e-6,18e-6,0.7,0.8,80,0.2,520,200,0);
mtf = computeMtfGrid(2.01e-6,18e-6,0.5,0.9,80,0.2,520,200,0);
verifyEqual(testCase,computeSnrGrid(2.01e-6,18e-6,0.7,0.8,80,0.2,520,200,snr),snr);
verifyEqual(testCase,computeMtfGrid(2.01e-6,18e-6,0.5,0.9,80,0.2,520,200,mtf),mtf);
verifyTrue(testCase,isnan(computeSnrGrid(2.01e-6,18e-6,0.7,0.8,80,0.2,520,200,snr+1)));
verifyTrue(testCase,isnan(computeMtfGrid(2.01e-6,18e-6,0.5,0.9,80,0.2,520,200,mtf+1)));
end

function testCompatibilityCallsCanSkipExports(testCase)
altitudes_km = [400 520]; pupilDiameters_mm = [20 28 200];
folder = tempname;
snr = SNRfunction(2.01e-6,18e-6,0.7,0.8,80,0.2,altitudes_km,pupilDiameters_mm,fullfile(folder,'snr'),'Example',1,400,false);
mtf = MTFfunction(2.01e-6,18e-6,0.5,0.9,80,0.2,altitudes_km,pupilDiameters_mm,fullfile(folder,'mtf'),'Example',1,0.25,false);
verifyFalse(testCase,isfolder(folder));
assertGrid(testCase,snr,computeSnrGrid(2.01e-6,18e-6,0.7,0.8,80,0.2,altitudes_km,pupilDiameters_mm,400));
assertGrid(testCase,mtf,computeMtfGrid(2.01e-6,18e-6,0.5,0.9,80,0.2,altitudes_km,pupilDiameters_mm,0.25));
end

function assertGrid(testCase,actual,expected)
if islogical(expected)
    verifyEqual(testCase,actual,expected);
    return;
end
verifyEqual(testCase,size(actual),size(expected));
verifyEqual(testCase,isnan(actual),isnan(expected));
mask = isfinite(expected);
verifyEqual(testCase,actual(mask),expected(mask),'AbsTol',1e-12,'RelTol',1e-12);
end
