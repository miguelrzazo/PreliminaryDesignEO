%% GSD sensitivity around the adopted mission design
projectDirectory=fileparts(fileparts(mfilename('fullpath')));
report=jsondecode(fileread(fullfile(projectDirectory,'Compiled data','Mission','summary.json')));
config=gsdSensitivityInputs(report.selected);
config.outputDirectory=fullfile(projectDirectory,'Compiled data');
config.figureDirectory=fullfile(projectDirectory,'Latex_Code','4.Payload');
results=runGsdSensitivityStudy(config);
