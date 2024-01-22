if ismac
%   datapath = '/Users/kloosterman/Library/CloudStorage/Dropbox/PROJECTS/LRaudio/data';
  datapath = '/Users/kloosterman/gridmaster2012/projectdata/LRaudio/data';
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
overwrite = 0;
mkdir(fileparts(datapath), 'mse')
for isub = 1:length(SUBJ)
  cfg.SUBJ = SUBJ{isub};
  cfg.datafile = fullfile(datapath, SUBJ{isub}, sprintf('clean_SUB%s', [SUBJ{isub} '.mat']));
  for icond = 1:2
    cfg.icond = icond;
    cfg.outpath = fullfile(fileparts(datapath), 'mse', sprintf('SUB%s_cond%d.mat', SUBJ{isub}, icond));
    if ~exist(cfg.outpath, 'file')
      cfglist{end+1} = cfg;
    else
      disp('File exists, skipping')
    end
  end      
end
%
% submit to tardis, or run locally
fun2run = @computemMSE;
if ismac
  cellfun(fun2run, cfglist, 'Uni', 0);
else
  qsubcellfun(fun2run, cfglist, 'memreq', 100e9, 'timreq', 23*60*60, 'stack', 1, ...
    'StopOnError', false, 'backend', 'slurm', 'options', ' --cpus-per-task=4 ');  
end

%% Merge mse files across subjects and cond
disp 'Merge mse'
mse_tmp = {};
trialinfo = [];
for isub = 1:length(SUBJ)
  cfg.SUBJ = SUBJ{isub};
  cfg.datafile = fullfile(datapath, SUBJ{isub}, sprintf('clean_SUB%s', [SUBJ{isub} '.mat']));
  for icond = 1:2
    cfg.icond = icond;
    cfg.outpath = fullfile(fileparts(datapath), 'mse', sprintf('SUB%s_cond%d.mat', SUBJ{isub}, icond));
    disp(cfg.outpath)
    if exist(cfg.outpath, 'file')      
      load(cfg.outpath)
      mse.freq = mse.timescales;
      mse.powspctrm = mse.sampen;
      mse.dimord = 'chan_freq_time';
      mse_tmp{isub,icond} = mse;
      trialinfo(isub,:,icond) = mean(mse.trialinfo);
    else
      disp('File not found, skipping')
    end
  end  
end
cfg=[];
cfg.keepindividual = 'yes';
mse_merged{1} = ft_freqgrandaverage(cfg, mse_tmp{:,1});
mse_merged{2} = ft_freqgrandaverage(cfg, mse_tmp{:,2});
cfg.keepindividual = 'no';
mse_merged{3} = ft_freqgrandaverage(cfg, mse_merged{1:2});
cfg=[];
cfg.operation = 'subtract';
cfg.parameter = 'powspctrm';
mse_merged{4} = ft_math(cfg, mse_merged{1}, mse_merged{2});
% add behavior
mse_merged{1}.trialinfo = trialinfo(:,:,1);
mse_merged{2}.trialinfo = trialinfo(:,:,2);
mse_merged{3}.trialinfo = mean(trialinfo, 3);
mse_merged{4}.trialinfo = trialinfo(:,:,1) - trialinfo(:,:,2);

%% plot mse
cfg=[];
cfg.layout = 'acticap-64ch-standard2.mat';
%     cfg.zlim = [0.8 1.2];
cfg.colorbar = 'yes';
% cfg.baseline = [-0.5 0];
% cfg.baselinetype = 'relchange';
cfg.zlim = 'maxabs';
% cfg.zlim = [1.17 1.23];
%     cfg.xlim = [-0.5 1.5];
ft_multiplotTFR(cfg,mse_merged{4})

% cfg.frequency = 80;
% ft_multiplotER(cfg,mse_merged{1:2})
%% run stats: correlation mMSE vs behavior


