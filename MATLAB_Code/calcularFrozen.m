%% Inputs from the selected design; edit config below for an auxiliary study.
studies = missionStudyInputs();
config = studies.frozen;
results = runFrozenOrbitStudy(config);
