if ismac
  datapath = '/Users/kloosterman/Library/CloudStorage/Dropbox/PROJECTS/LRaudio/data';
  toolspath = '/Users/kloosterman/Documents/GitHub/';
else
  datapath = '/home/mpib/kloosterman/projectdata/LRaudio/data';
  toolspath = '/home/mpib/kloosterman/GitHub/';
end

restoredefaultpath
addpath(fullfile(toolspath, 'fieldtrip')); ft_defaults
addpath(fullfile(toolspath, 'mMSE'));
addpath(fullfile(toolspath, 'LRaudio'))
addpath(fullfile(toolspath, 'qsub-tardis')) 

outputpath = fullfile(fileparts(datapath), 'outputdata');
mkdir(outputpath)

%% define subjects
SUBJ={};
for i=1:36
  SUBJ{end+1} = sprintf('%d', i);
end
SUBJbool = true(size(SUBJ));
SUBJbool([10, 12, 15, 17]) = false; % exclude 10 and 15, 17 missing, 12 conditions missing
SUBJ = SUBJ(SUBJbool);
disp(SUBJ)

%% compute entropy
% make list of files to analyze on tardis
cfglist = {};
cfg=[];
mkdir(fileparts(datapath), 'mse')
for isub = 1:length(SUBJ)
  cfg.SUBJ = SUBJ{isub};
  cfg.datafile = fullfile(datapath, SUBJ{isub}, sprintf('clean_SUB%s', [SUBJ{isub} '.mat']));
  for icond = 1:2
    cfg.icond = icond;
    cfg.outpath = fullfile(fileparts(datapath), 'mse', sprintf('SUB%s_cond%d.mat', SUBJ{isub}, icond));
    cfglist{isub,icond} = cfg;
  end      
end
%
% submit to tardis, or run locally
fun2run = @computemMSE;
if ismac
  mse = cellfun(fun2run, cfglist, 'Uni', 0);
else
  mse = qsubcellfun(fun2run, cfglist, 'memreq', 40e9, 'timreq', 12*60*60, 'stack', 1, ...
    'StopOnError', false, 'backend', 'slurm', 'options', ' --cpus-per-task=4 ');  
end
mse

%%

