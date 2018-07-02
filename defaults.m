addpath('./exportfig/');
%addpath('./distributionPlot/');
set(0,'defaultAxesFontName', 'Arial')
set(0,'defaultTextFontName', 'Arial')
set(0,'DefaultAxesFontSize', 16);
set(0,'DefaultTextFontSize', 16);


npyrs=240;
nneurons=300;
stimduration=4000;
nbranches=20;
ninputs=10;

npyrbranches=npyrs*nbranches;

nruns=10;
CUTOFF=10; % Hz
%results = containers.Map();
