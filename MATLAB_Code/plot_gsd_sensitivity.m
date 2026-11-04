%% Editable GSD sensitivity study inputs
projectDirectory=fileparts(fileparts(mfilename('fullpath')));
config=struct();
config.inputFile=fullfile(projectDirectory,'Compiled data','GSD_sensitivity_summary.csv');
config.outputDirectory=fullfile(projectDirectory,'Compiled data');
config.figureDirectory=fullfile(projectDirectory,'Latex_Code','4.Payload');
config.displaySatelliteCount=2;
config.updateInputSummary=false;
config.writeOutputs=true;
results=runGsdSensitivityStudy(config);
