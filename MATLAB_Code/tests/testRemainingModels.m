function tests=testRemainingModels
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
function testDryMassScalarReference(testCase)
fields={'dryMass_kg','thematicMapperDryMass_kg','seosatDryMass_kg','instrumentLength_m', ...
    'instrumentArea_m2','instrumentVolume_m3','instrumentMass_kg','instrumentPower_W', ...
    'satelliteVolume_m3','satelliteLength_m','satelliteArea_m2','cubesatUnits'};
for telescopeCount=[1 2]
    for columnAxis=[false true]
        diameters_mm=[28 90 125 200 250];
        if columnAxis, diameters_mm=diameters_mm.'; end
        expected=cell(1,12); [expected{:}]=referenceDryMass(diameters_mm,telescopeCount);
        actual=computeDryMass(diameters_mm,telescopeCount);
        for fieldIndex=1:numel(fields)
            assertArray(testCase,actual.(fields{fieldIndex}),expected{fieldIndex});
        end
    end
end
end
function testTotalMassScalarReference(testCase)
fields={'totalMass_kg','fuelMass_kg','impulseCounts','totalDeltaV_m_s', ...
    'cycleTimes_s','cycleDeltaV_m_s','unmarginedFuelMass_kg'};
for upperFactor=[1 1.02]
    expected=cell(1,7);
    [expected{:}]=referenceTotalMass([400;520],[2;2.5],[0.03;0.04],0.98,upperFactor);
    actual=computeTotalMass([400;520],[2;2.5],[0.03;0.04],0.98,upperFactor);
    for fieldIndex=1:numel(fields)
        assertArray(testCase,actual.(fields{fieldIndex}),expected{fieldIndex});
    end
end
% Scalar mass/area broadcast preserves the altitude shape.
actual=computeTotalMass([400 520],2,0);
verifyEqual(testCase,actual.totalMass_kg,[2 2]);
verifyEqual(testCase,actual.cycleTimes_s,[Inf Inf]);
verifyError(testCase,@() computeTotalMass([400 520],[1 2 3],0.02),'eo:MassInputSize');
end
function testSelectionAndConfigurationReference(testCase)
[data,params]=fixture();
expected=referenceSelection(params.alturas_orbitales,params.swaths_km,params.diametros_pupila, ...
    params.telescope_names,params.fov_limit,params.satellite_configs, ...
    params.MTF_Requirement,params.SNR_Requirement,params.Cov_Requirement,data,false);
actual=CrossDataFunction(params.alturas_orbitales,params.swaths_km,params.diametros_pupila, ...
    params.telescope_names,params.fov_limit,params.satellite_configs, ...
    params.MTF_Requirement,params.SNR_Requirement,params.Cov_Requirement,data,false);
verifyEqual(testCase,actual,expected);
data.selection=actual;
% The regional filter intentionally changes the old single-latitude schema.
params.regionalRevisitSolver = @(varargin) 1;
actual=ConfigAnalysis(params,data,false);
viable=actual.table(logical(actual.table.Es_Viable),:);
verifyNotEmpty(testCase,viable);
verifyEqual(testCase,actual.best.mass,min(viable.Masa_Total_Constelacion_kg));
verifyTrue(testCase,all(viable.Swath_Adquirido_km>=viable.Swath_Requerido_km));
verifyTrue(testCase,all(viable.Swath_Adquirido_km<=viable.Swath_Real_Sistema_km));
verifyTrue(testCase,all(viable.Detectores_Necesarios>=viable.Num_Telescopios));
verifyTrue(testCase,all(viable.Dias_Cobertura<=params.Cov_Requirement));
verifyNotEmpty(testCase,actual.maintenance);
end
function testFirstInputOrderAndEmptyRows(testCase)
[viable,theoretical]=selectMinimumDiameters([0.3 0.4;NaN NaN],[500 600;500 600], ...
    [7 2;2 3],[200 520],[90 28],[20 40],6,1,1,0.25,400,7);
verifyEqual(testCase,viable(1),90); verifyTrue(testCase,isnan(viable(2)));
verifyEqual(testCase,theoretical(1),90);
end
function [data,params]=fixture()
params=struct('alturas_orbitales',[400 520],'diametros_pupila',[28 90], ...
    'swaths_km',[20 50],'telescope_names',{{'Refractivo','Korsch'}},'fov_limit',[10 3], ...
    'satellite_configs',[1 1;2 2;1 1],'MTF_Requirement',0.25,'SNR_Requirement',400, ...
    'Cov_Requirement',7,'GSD',80,'Npix_detectors',[2048 2048 2048], ...
    'Pixel_size',[18e-6 18e-6 18e-6],'LTAN_hour',6, ...
    'maintenance_lower_factor',0.98,'maintenance_upper_factor',1.02);
data.axes=struct('altitudes_km',params.alturas_orbitales(:), ...
    'pupilDiameters_mm',params.diametros_pupila,'swaths_km',params.swaths_km, ...
    'telescopeNames',{params.telescope_names},'satelliteConfigurations',params.satellite_configs);
data.optics.mtf=cell(3,3,2); data.optics.snr=cell(3,3,2); data.coverage.days=cell(3,3,2);
for detectorIndex=1:3
    for telescopeIndex=1:2
        data.optics.mtf{2,detectorIndex,telescopeIndex}=[0.3 0.4;0.1 0.3];
        data.optics.snr{3,detectorIndex,telescopeIndex}=[500 600;500 600];
        for configurationIndex=1:3
            data.coverage.days{configurationIndex,detectorIndex,telescopeIndex}=[3 2;NaN NaN];
        end
        data.coverage.days{3,detectorIndex,telescopeIndex}(:)=NaN;
    end
end
end
function assertArray(testCase,actual,expected)
verifyEqual(testCase,size(actual),size(expected));
verifyEqual(testCase,isnan(actual),isnan(expected));
verifyEqual(testCase,isinf(actual),isinf(expected));
mask=isfinite(expected);
verifyEqual(testCase,actual(mask),expected(mask),'AbsTol',1e-10,'RelTol',1e-10);
end
