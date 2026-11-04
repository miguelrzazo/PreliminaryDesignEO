function tests = testCoverageModel
%TESTCOVERAGEMODEL Check adapter behavior independently of toolkit physics.
% Deterministic solvers validate the refactor, not the scientific revisit model.
tests = functiontests(localfunctions);
end

function setupOnce(testCase)
root = fileparts(fileparts(mfilename('fullpath')));
testCase.TestData.oldPath = path;
addpath(fullfile(root,'lib'), fullfile(root,'tests','reference'));
end

function teardownOnce(testCase)
path(testCase.TestData.oldPath);
end

function testScalarRegression(testCase)
for telescopeCount = [1 2]
    for cloudFraction = [0 1/6]
        inputs = {80,[200 520 100000],[20 50 500],7,2048,10,2,telescopeCount,0.05,cloudFraction,3};
        % Depend on both orbital elements and FOV to detect indexing mistakes.
        solver = @(orbitalElements,varargin) orbitalElements(1)/1e6 + varargin{10}/100;
        expected = referenceCoverage(inputs{:},solver);
        actual = computeCoverageGrid(inputs{:},solver);
        verifyEqual(testCase,size(actual),[3 3]);
        verifyEqual(testCase,isnan(actual),isnan(expected));
        validMask = isfinite(expected);
        verifyEqual(testCase,actual(validMask),expected(validMask),'AbsTol',1e-12,'RelTol',1e-12);
    end
end
end

function testOnlyFeasibleCellsReachSolverInOriginalOrder(testCase)
solverCalls = zeros(0,2);
[grid,diagnostics] = computeCoverageGrid(80,[200 520],[20 50 500],7,2048,10,2,2,0.05,0,3,@recordSolver);
expectedAltitudes_m = 6378.137e3 + [200;200;520;520]*1000;
expectedHalfFov_deg = atand([20;50;20;50]*0.95 ./ (2*[200;200;520;520]));
verifyEqual(testCase,solverCalls(:,1),expectedAltitudes_m);
verifyEqual(testCase,solverCalls(:,2),expectedHalfFov_deg,'AbsTol',1e-12);
verifyEqual(testCase,nnz(diagnostics.geometryFeasible),4);
verifyEqual(testCase,grid(:,1:2),ones(2,2));
verifyTrue(testCase,all(isnan(grid(:,3))));

    function revisit_days = recordSolver(orbitalElements,varargin)
        verifyEqual(testCase,cell2mat(varargin(1:8)),[25 1 0 2 1 0 1 0]);
        verifyEqual(testCase,varargin{9},[1 7]);
        verifyEqual(testCase,orbitalElements([2 4 5]),[0 0 0]);
        solverCalls(end+1,:) = [orbitalElements(1),varargin{10}];
        revisit_days = 1;
    end
end

function testCloudCorrectionAndRequirementBoundary(testCase)
inputs = {80,520,20,7,2048,10,2,2,0.05};
[grid,diagnostics] = computeCoverageGrid(inputs{:},0,3,@(varargin) 7);
verifyEqual(testCase,grid,7);
verifyFalse(testCase,diagnostics.requirementExceeded);
[grid,diagnostics] = computeCoverageGrid(inputs{:},1/6,3,@(varargin) 7);
verifyTrue(testCase,isnan(grid));
verifyTrue(testCase,diagnostics.requirementExceeded);
grid = computeCoverageGrid(inputs{:},1/6,3,@(varargin) 2);
verifyEqual(testCase,grid,2/(1-1/6));
end

function testSolverFailureDoesNotStopOtherCells(testCase)
[grid,diagnostics] = computeCoverageGrid(80,[200 520],20,7,2048,10,2,2,0.05,0,3,@sometimesFailingSolver);
verifyTrue(testCase,isnan(grid(1)));
verifyEqual(testCase,grid(2),1);
verifyEqual(testCase,diagnostics.solverErrorIdentifiers{1},'test:CoverageFailure');
verifyNotEmpty(testCase,diagnostics.solverErrors{1});
verifyEmpty(testCase,diagnostics.solverErrors{2});
end

function revisit_days = sometimesFailingSolver(orbitalElements,varargin)
if orbitalElements(1) < 6800e3
    error('test:CoverageFailure','Deliberate solver failure');
end
revisit_days = 1;
end

function testNonfiniteResultIsVisible(testCase)
[grid,diagnostics] = computeCoverageGrid(80,520,20,7,2048,10,2,2,0.05,0,3,@(varargin) NaN);
verifyTrue(testCase,isnan(grid));
verifyTrue(testCase,diagnostics.nonfiniteSolverResult);
verifyEmpty(testCase,diagnostics.solverErrors{1});
end

function testSingletonAndColumnAxes(testCase)
inputs = {80,[200;520],[20;50;500],7,2048,10,2,2,0.05,0,3,@(varargin) 1};
grid = computeCoverageGrid(inputs{:});
verifyEqual(testCase,size(grid),[2 3]);
grid = computeCoverageGrid(80,520,[20 50],7,2048,10,2,2,0.05,0,3,@(varargin) 1);
verifyEqual(testCase,size(grid),[1 2]);
end

function testFullyRejectedGeometrySkipsSolver(testCase)
[grid,diagnostics] = computeCoverageGrid(80,100000,[20 50],7,2048,10,2,2,0.05,0,3,@unexpectedSolverCall);
verifyTrue(testCase,all(isnan(grid)));
verifyFalse(testCase,any(diagnostics.geometryFeasible));
verifyTrue(testCase,all(cellfun(@isempty,diagnostics.solverErrors)));
end

function result = unexpectedSolverCall(varargin)
error('test:UnexpectedSolverCall','Rejected cells must not call the solver');
result = NaN;
end

function testCompatibilityCallCanSkipExports(testCase)
% The swath fails geometry, so the external toolkit is unnecessary.
oldDirectory = pwd; oldWarnings = warning;
[grid,diagnostics] = CoverageRevisitCalc(80,520,1e6,7,2048,10,2,2,1,'Example',0.05,0,3,false);
verifyTrue(testCase,isnan(grid));
verifyFalse(testCase,diagnostics.geometryFeasible);
verifyEqual(testCase,pwd,oldDirectory);
currentWarnings = warning;
verifyEqual(testCase,currentWarnings,oldWarnings);
end

function testExportRoundTripAndFigureCleanup(testCase)
folder = tempname; mkdir(folder);
oldDirectory = pwd;
cleanup = onCleanup(@() cleanExportFolder(folder,oldDirectory));
cd(folder);
previousFigures = numel(findall(groot,'Type','figure'));
[grid,diagnostics] = computeCoverageGrid(80,[200 520],[20 50 500],7,2048,10,1,1,0.05,0,3,@(varargin) 1);
configurationName = '1 Satelite(s); 1 Telescopio(s): Example; Detector 1';
writeCoverageResult(grid,[200 520],[20 50 500],7,diagnostics,configurationName);
actual = readmatrix(fullfile('coverage',['coverage_' configurationName '.csv']));
verifyEqual(testCase,isnan(actual),isnan(grid));
validMask = isfinite(grid);
verifyEqual(testCase,actual(validMask),grid(validMask),'AbsTol',1e-12);
verifyTrue(testCase,isfile(fullfile('coverage',['heatmap_' configurationName '.jpg'])));
verifyTrue(testCase,isfile(fullfile('coverage',['heatmap_' configurationName '.pdf'])));
verifyEqual(testCase,numel(findall(groot,'Type','figure')),previousFigures);
end

function cleanExportFolder(folder,oldDirectory)
cd(oldDirectory);
if isfolder(folder), rmdir(folder,'s'); end
end

function testInvalidHorizonFailsBeforeSolver(testCase)
% The external solver compares whole-day iterations against the horizon.
verifyError(testCase,@() computeCoverageGrid(80,520,20,7.5,2048,10,2,2,0.05,0,3,@unexpectedSolverCall), ...
    'eo:WholeDayHorizon');
end
