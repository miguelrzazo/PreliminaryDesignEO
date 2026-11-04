function tests = testAnalysisHandoff
%TESTANALYSISHANDOFF Verify memory handoffs against the existing CSV path.
tests = functiontests(localfunctions);
end

function setupOnce(testCase)
root = fileparts(fileparts(mfilename('fullpath')));
testCase.TestData.oldPath = path;
addpath(fullfile(root,'lib'));
end

function teardownOnce(testCase)
path(testCase.TestData.oldPath);
end

function testMemoryPathNeedsNoIntermediateFiles(testCase)
[data,params] = fixture();
folder = tempname; mkdir(folder);
oldDirectory = pwd;
cleanup = onCleanup(@() cleanFolder(folder,oldDirectory));
cd(folder);
data.selection = select(data,params);
summary = ConfigAnalysis(params,data,false);
verifyEqual(testCase,data.selection{1}.viableDiameters_mm(:,1,1),[28;28]);
verifyNotEmpty(testCase,summary.rows);
verifyFalse(testCase,isfolder('MTF'));
verifyFalse(testCase,isfolder('SNR'));
verifyFalse(testCase,isfolder('coverage'));
verifyFalse(testCase,isfolder('HvsDmin'));
verifyFalse(testCase,isfolder('OptimumConfigs'));
end

function testCsvAndMemoryPathsAgree(testCase)
[data,params] = fixture();
folder = tempname; mkdir(folder);
oldDirectory = pwd;
cleanup = onCleanup(@() cleanFolder(folder,oldDirectory));
cd(folder);
data.selection = select(data,params);
writeFixture(data,params);
csvSelection = CrossDataFunction(params.alturas_orbitales,params.swaths_km, ...
    params.diametros_pupila,params.telescope_names,params.fov_limit, ...
    params.satellite_configs,params.MTF_Requirement,params.SNR_Requirement, ...
    params.Cov_Requirement,[],false);
verifyEqual(testCase,csvSelection,data.selection);
expected = ConfigAnalysis(params,[],false);
actual = ConfigAnalysis(params,data,false);
assertSummary(testCase,actual,expected);
% Stale files must not participate when current-run data is supplied.
writematrix(nan(2,2),fullfile('coverage',coverageFilename(1)));
verifyEqual(testCase,select(data,params),data.selection);
assertSummary(testCase,ConfigAnalysis(params,data,false),actual);
end

function testCsvCoordinateOrderIsValidated(testCase)
[data,params] = fixture();
folder = tempname; mkdir(folder);
oldDirectory = pwd;
cleanup = onCleanup(@() cleanFolder(folder,oldDirectory));
cd(folder);
data.selection = select(data,params);
writeFixture(data,params);
filename = 'MTF/MTF_Lambda2_Detector1_Telescopio1_resultados.csv';
opticalTable = readtable(filename,'ReadRowNames',true);
opticalTable = opticalTable(end:-1:1,:);
writetable(opticalTable,filename,'WriteRowNames',true);
verifyError(testCase,@() CrossDataFunction(params.alturas_orbitales,params.swaths_km, ...
    params.diametros_pupila,params.telescope_names,params.fov_limit, ...
    params.satellite_configs,params.MTF_Requirement,params.SNR_Requirement, ...
    params.Cov_Requirement,[],false),'eo:ImportedAxisMismatch');
end

function testReorderedAxesAreRejected(testCase)
[data,params] = fixture();
data.axes.altitudes_km = flipud(data.axes.altitudes_km);
verifyError(testCase,@() select(data,params),'eo:ResultAxisMismatch');
verifyError(testCase,@() ConfigAnalysis(params,data,false),'eo:ResultAxisMismatch');
end

function testMalformedMemoryCannotFallBackToCsv(testCase)
[data,params] = fixture();
data.optics.mtf{2,1,1} = ones(1,2);
verifyError(testCase,@() select(data,params),'eo:ResultDimensions');
data.selection = {struct('viableDiameters_mm',28*ones(2,1,3))};
verifyError(testCase,@() ConfigAnalysis(params,data,false),'eo:ResultDimensions');
end

function testNoSolutionsReturnsEmptySummary(testCase)
[data,params] = fixture();
for detectorIndex = 1:3
    data.optics.mtf{2,detectorIndex,1}(:) = NaN;
end
data.selection = select(data,params);
summary = ConfigAnalysis(params,data,false);
verifyEmpty(testCase,summary.rows);
verifyEqual(testCase,summary.best.mass,Inf);
end

function selection = select(data,params)
selection = CrossDataFunction(params.alturas_orbitales,params.swaths_km, ...
    params.diametros_pupila,params.telescope_names,params.fov_limit, ...
    params.satellite_configs,params.MTF_Requirement,params.SNR_Requirement, ...
    params.Cov_Requirement,data,false);
end

function [data,params] = fixture()
params = struct('alturas_orbitales',[400 520],'diametros_pupila',[28 90], ...
    'swaths_km',[20 50],'telescope_names',{{'Example'}},'fov_limit',10, ...
    'satellite_configs',[1 1],'MTF_Requirement',0.25,'SNR_Requirement',400, ...
    'Cov_Requirement',7,'GSD',80,'Npix_detectors',[2048 2048 2048], ...
    'Pixel_size',[18e-6 18e-6 18e-6],'LTAN_hour',6, ...
    'maintenance_lower_factor',0.98,'maintenance_upper_factor',1.02);
params.regionalRevisitSolver = @(varargin) 1;
data.axes = struct('altitudes_km',params.alturas_orbitales(:), ...
    'pupilDiameters_mm',params.diametros_pupila, 'swaths_km',params.swaths_km, ...
    'telescopeNames',{params.telescope_names}, 'satelliteConfigurations',params.satellite_configs);
data.optics.mtf = cell(3,3,1);
data.optics.snr = cell(3,3,1);
data.coverage.days = cell(1,3,1);
for detectorIndex = 1:3
    data.optics.mtf{2,detectorIndex,1} = 0.3*ones(2,2);
    data.optics.snr{3,detectorIndex,1} = 500*ones(2,2);
    data.coverage.days{1,detectorIndex,1} = 3*ones(2,2);
end
end

function writeFixture(data,params)
mkdir('MTF'); mkdir('SNR'); mkdir('coverage'); mkdir('HvsDmin');
for detectorIndex = 1:3
    mtf = array2table(data.optics.mtf{2,detectorIndex,1}, ...
        'RowNames',{'Alt_400','Alt_520'},'VariableNames',{'Diam_28','Diam_90'});
    snr = array2table(data.optics.snr{3,detectorIndex,1}, ...
        'RowNames',{'Alt_400','Alt_520'},'VariableNames',{'Diam_28','Diam_90'});
    writetable(mtf,sprintf('MTF/MTF_Lambda2_Detector%d_Telescopio1_resultados.csv',detectorIndex),'WriteRowNames',true);
    writetable(snr,sprintf('SNR/SNR_Lambda3_Detector%d_Telescopio1_resultados.csv',detectorIndex+3),'WriteRowNames',true);
    writematrix(data.coverage.days{1,detectorIndex,1},fullfile('coverage',coverageFilename(detectorIndex)));
    diameters = data.selection{1}.viableDiameters_mm(:,1,detectorIndex);
    selected = array2table([params.alturas_orbitales(:),diameters], ...
        'VariableNames',{'Altura_km','Dmin_Example_mm'});
    writetable(selected,sprintf('HvsDmin/HvsDmin_1sat_1tel_Detector%d.csv',detectorIndex));
end
end

function filename = coverageFilename(detectorIndex)
filename = sprintf('coverage_1 Satelite(s); 1 Telescopio(s): Example; Detector %d.csv',detectorIndex);
end

function assertSummary(testCase,actual,expected)
verifyEqual(testCase,size(actual.rows),size(expected.rows));
for column = 1:size(expected.rows,2)
    if column == 6 || column == 1
        verifyEqual(testCase,actual.rows(:,column),expected.rows(:,column));
    else
        verifyEqual(testCase,cell2mat(actual.rows(:,column)), ...
            cell2mat(expected.rows(:,column)),'AbsTol',1e-12,'RelTol',1e-12);
    end
end
verifyEqual(testCase,actual.best.mass,expected.best.mass,'RelTol',1e-12);
end

function cleanFolder(folder,oldDirectory)
cd(oldDirectory);
if isfolder(folder), rmdir(folder,'s'); end
end
