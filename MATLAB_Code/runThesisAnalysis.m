function report = runThesisAnalysis()
%RUNTHESISANALYSIS Reproduce the design, mission budgets and thesis figures.
scriptDirectory = fileparts(mfilename('fullpath'));
root = fileparts(scriptDirectory);
previousPath = path; pathCleanup = onCleanup(@() path(previousPath));
addpath(scriptDirectory,fullfile(scriptDirectory,'lib'));
stage = tempname; mkdir(stage);
stageCleanup = onCleanup(@() rmdir(stage,'s'));
pipelineDirectory = fullfile(stage,'Pipeline');
previousOutput = getenv('TFG_OUTPUT_ROOT');
environmentCleanup = onCleanup(@() setenv('TFG_OUTPUT_ROOT',previousOutput));
setenv('TFG_OUTPUT_ROOT',pipelineDirectory);
previousVisibility = get(groot,'defaultFigureVisible');
graphicsCleanup = onCleanup(@() set(groot,'defaultFigureVisible',previousVisibility));
set(groot,'defaultFigureVisible','off');
analysis = runPipeline(scriptDirectory);
assert(analysis.coverage.executionSucceeded,'eo:CoverageExecution','Coverage solver failures occurred.');
solutions = analysis.configuration.table;
viable = sortrows(solutions(logical(solutions.Es_Viable),:),'Masa_Total_Constelacion_kg');
assert(~isempty(viable),'eo:NoViableDesign','No design satisfies the mission requirements.');
best = viable(1,:); best.GSD = analysis.config.groundSampleDistance_m;
design = table2struct(best);
missionDirectory = fullfile(stage,'Mission'); mkdir(missionDirectory);
writetable(best,fullfile(missionDirectory,'selected_design.csv'));
writetable(viable,fullfile(missionDirectory,'viable_designs.csv'));
inputs = missionStudyInputs(design,missionDirectory);
frozen = runFrozenOrbitStudy(inputs.frozen);
contacts = runContactStudy(inputs.contact);
keplerianInputs = inputs.contact;
keplerianInputs.propagator = 'two-body-keplerian'; keplerianInputs.writeOutputs = false;
keplerian = runContactStudy(keplerianInputs);
writetable(keplerian.summary,fullfile(missionDirectory,'Contacts','contacts_two-body-keplerian.csv'));

assert(all(contacts.minimumRate_Mbps <= inputs.contact.referenceRate_Mbps), ...
    'eo:DownlinkRequirement','The assumed station rate cannot download the acquired data.');
assert(all(max(contacts.memory.required_GB,[],2) <= inputs.contact.memoryCapacity_GB), ...
    'eo:MemoryRequirement','The selected memory capacity is insufficient.');
endOfLife = runDeorbitStudy(inputs.deorbit);
rgt = runRgtStudy(inputs.rgt);
latitudes = (25:2:49).';
revisit = zeros(size(latitudes));
a = 6378.137e3+design.Altura_km*1000;
inc = deg2rad(inputs.contact.inclination_deg);
psi = atand((1-analysis.config.overlapFraction)*design.Swath_Adquirido_km/(2*design.Altura_km));
for latitudeIndex = 1:numel(latitudes)
    revisit(latitudeIndex) = RevisitCalc([a,0,inc,0,0],latitudes(latitudeIndex), ...
        1,0,design.Num_Satelites,1,0,1,0,[1,analysis.config.revisitRequirement_days],psi) / ...
        (1-analysis.config.cloudFraction);
end
assert(all(isfinite(revisit) & revisit>0 & revisit<=analysis.config.revisitRequirement_days), ...
    'eo:RegionalRequirement','The selected design fails the sampled regional requirement.');
writetable(table(latitudes,revisit,'VariableNames',{'Latitude_deg','EffectiveRevisit_days'}), ...
    fullfile(missionDirectory,'latitude_revisit.csv'));
report = struct('matlab_version',version,'selected',design, ...
    'telescope_index',find(strcmp(analysis.config.telescopeNames,design.Tipo_Telescopio),1), ...
    'revisit_requirement_days',analysis.config.revisitRequirement_days, ...
    'cloud_fraction',analysis.config.cloudFraction,'overlap_fraction',analysis.config.overlapFraction, ...
    'latitudes_deg',latitudes,'revisit_days',revisit, ...
    'selected_frozen',table2struct(frozen.summary), ...
    'keplerian',table2struct(keplerian.summary), ...
    'sgp4',table2struct(contacts.summary),'end_of_life',table2struct(endOfLife.summary), ...
    'rgt',rgt.nearestSolution);
fid = fopen(fullfile(missionDirectory,'summary.json'),'w');
assert(fid>=0,'eo:OutputFile','Cannot open mission summary.');
fprintf(fid,'%s\n',jsonencode(report)); fclose(fid);
publishThesisFigures(stage);
gsdInputs = gsdSensitivityInputs(design);
gsdInputs.outputDirectory = fullfile(root,'Compiled data');
gsdInputs.figureDirectory = fullfile(root,'Latex_Code','4.Payload');
runGsdSensitivityStudy(gsdInputs);
dataDirectory = fullfile(root,'Compiled data');
copyfile(pipelineDirectory,fullfile(dataDirectory,'Pipeline'),'f');
copyfile(missionDirectory,fullfile(dataDirectory,'Mission'),'f');
fprintf('THESIS_ANALYSIS_COMPLETE\n'); disp(best);
end

function analysis = runPipeline(scriptDirectory)
run(fullfile(scriptDirectory,'master.m'));
analysis = results;
end
