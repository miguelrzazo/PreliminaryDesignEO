function timings = benchmarkOpticalRefactor()
%BENCHMARKOPTICALREFACTOR Compare equivalent compute-only SNR workloads.
% Run numerical tests before interpreting timings. Results depend on hardware.
root = fileparts(fileparts(mfilename('fullpath')));
oldPath = path;
cleanup = onCleanup(@() path(oldPath));
addpath(root, fullfile(root,'lib'), fullfile(root,'tests','reference'));
args = {2.01e-6,18e-6,0.7,0.8,80,0.2,200:10:1000,20:200};
timings.scalarSeconds = timeit(@() referenceSNRfunction(args{:},'','',1,400));
timings.vectorSeconds = timeit(@() computeSnrGrid(args{:},400));
timings.speedRatio = timings.scalarSeconds / timings.vectorSeconds;
end
