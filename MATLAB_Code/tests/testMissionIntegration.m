function tests = testMissionIntegration
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

function testSunSynchronousRateBalance(testCase)
computeCoverageGrid(80,630,190,7,2048,10,2,2,.05,1/6,3,@checkRate);
    function revisit = checkRate(elements,varargin)
        a = elements(1); inclination = elements(3);
        nodalRate = -1.5*sqrt(3.986004418e14/a^3)*1.082629989052e-3* ...
            (6378.137e3/a)^2*cos(inclination);
        verifyEqual(testCase,nodalRate,2*pi/(365.2421897*86400),'AbsTol',1e-18);
        revisit = 1;
    end
end

function testRegionalWorstLatitudeControlsWidth(testCase)
params = struct('GSD',80,'Cov_Requirement',7,'solapamiento',.05, ...
    'cobertura_nubes',0,'regionalRevisitSolver',@worstAt41);
cache = containers.Map('KeyType','char','ValueType','double');
[width,revisit] = selectRegionalSwath(690,2,2,[190 240],params,cache);
verifyEqual(testCase,width,240);
verifyEqual(testCase,revisit,6);
    function days = worstAt41(elements,latitude,varargin)
        acquired = 2*(elements(1)/1000-6378.137)*tand(varargin{9})/.95;
        if latitude==41 && acquired<240-1e-8, days=8; else, days=6; end
    end
end

function testFrozenArcsAndJointEnvelopeRemainDistinct(testCase)
design = struct('Altura_km',690,'Num_Satelites',2,'Num_Telescopios',2, ...
    'GSD',80,'Swath_Adquirido_km',240,'Masa_Seca_Satelite_kg',3, ...
    'Area_Satelite_m2',.1,'Masa_Combustible_kg',.05);
inputs = missionStudyInputs(design,tempname);
inputs.frozen.writeOutputs = false;
actual = runFrozenOrbitStudy(inputs.frozen);
verifyEqual(testCase,actual.summary.MaxAscendingArcVariation_km(1), ...
    actual.summary.MaxDescendingArcVariation_km(1),'AbsTol',1e-9);
verifyTrue(testCase,all(actual.circular.AscendingArcVariation_km <= actual.circular.CombinedEnvelope_km+1e-10));
verifyLessThan(testCase,max(actual.circular.AscendingArcVariation_km),max(actual.circular.CombinedEnvelope_km));
verifyFalse(testCase,isfolder(inputs.frozen.outputDirectory));
end
