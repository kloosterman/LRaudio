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

%% make eeg layout
load('acticap-64ch-standard2.mat')
lay.label(find(contains(lay.label, 'Ref'))) = {'FCz'};

%% define subjects
SUBJ={};
for i=1:36
  SUBJ{end+1} = sprintf('%d', i);
end
SUBJbool = true(size(SUBJ));
SUBJbool([10, 12, 15, 17]) = false; % exclude 10 and 15, 17 missing, 12 conditions missing
SUBJ = SUBJ(SUBJbool);
disp(SUBJ)

%% compute entropy, freqanalysis and ERPs
% make list of files to analyze on tardis
cfglist = {};
cfg=[];
cfg.analysis = 'mse'; % freq, mse, or erp
cfg.evoked = 'regress'; % empty, regress, or subtract
cfg.csd = 'csd'; % empty or csd
mkdir(fullfile(fileparts(datapath), cfg.analysis, cfg.evoked, cfg.csd))
overwrite = 1;
for isub = 1:length(SUBJ)
  cfg.SUBJ = SUBJ{isub};
  cfg.datafile = fullfile(datapath, SUBJ{isub}, sprintf('clean_SUB%s', [SUBJ{isub} '.mat']));
  for icond = 1:2
    cfg.icond = icond;
    cfg.outpath = fullfile(fileparts(datapath), cfg.analysis, cfg.evoked, cfg.csd, sprintf('SUB%s_cond%d.mat', SUBJ{isub}, icond));
    if overwrite || ~exist(cfg.outpath, 'file')  
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
else % mse: 'memreq', 100e9, 'timreq', 23*60*60, 'options', ' --cpus-per-task=4 '
  qsubcellfun(fun2run, cfglist, 'memreq', 100e9, 'timreq', 23*60*60, 'stack', 1, ...
    'StopOnError', false, 'backend', 'slurm', 'options', ' --cpus-per-task=4 ');  
  return
end

%% Merge mse files across subjects and cond
disp 'Merge analyses'
mse_tmp = {}; freq_tmp = {}; timelock_tmp = {};
trialinfo = [];
for isub = 1:length(SUBJ)
%   datafile = fullfile(datapath, SUBJ{isub}, sprintf('clean_SUB%s', [SUBJ{isub} '.mat']));
  for icond = 1:2
    path = fullfile(fileparts(datapath), 'mse', 'regress', sprintf('SUB%s_cond%d.mat', SUBJ{isub}, icond));
    disp(path)
    if exist(path, 'file')
      load(path)
      mse.freq = mse.timescales;
      mse.powspctrm = mse.sampen;
      mse.dimord = 'chan_freq_time';
      mse_tmp{isub,icond} = mse;
      trialinfo(isub,:,icond) = mean(mse.trialinfo);
    else
      disp('File not found, skipping')
    end
    path = fullfile(fileparts(datapath), 'freq', sprintf('SUB%s_cond%d.mat', SUBJ{isub}, icond));
    disp(path)
    if exist(path, 'file')
      load(path)
      cfg=[];
      cfg.baseline = [-0.5 0];
      cfg.baselinetype = 'relchange';
      freq = ft_freqbaseline(cfg, freq);
      freq_tmp{isub,icond} = freq;
    else
      disp('File not found, skipping')
    end
    path = fullfile(fileparts(datapath), 'erp', sprintf('SUB%s_cond%d.mat', SUBJ{isub}, icond));
    disp(path)
    if exist(path, 'file')
      load(path)
      timelock_tmp{isub,icond} = timelock;
    else
      disp('File not found, skipping')
    end
  end
end
%% combine mse 
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

%% combine freq 
cfg=[];
cfg.keepindividual = 'yes';
freq_merged{1} = ft_freqgrandaverage(cfg, freq_tmp{:,1});
freq_merged{2} = ft_freqgrandaverage(cfg, freq_tmp{:,2});
cfg.keepindividual = 'no';
freq_merged{3} = ft_freqgrandaverage(cfg, freq_merged{1:2});
cfg=[];
cfg.operation = 'subtract';
cfg.parameter = 'powspctrm';
freq_merged{4} = ft_math(cfg, freq_merged{1}, freq_merged{2});
% add behavior
freq_merged{1}.trialinfo = trialinfo(:,:,1);
freq_merged{2}.trialinfo = trialinfo(:,:,2);
freq_merged{3}.trialinfo = mean(trialinfo, 3);
freq_merged{4}.trialinfo = trialinfo(:,:,1) - trialinfo(:,:,2);
%% combine timelock 
cfg=[];
cfg.keepindividual = 'yes';
timelock_merged{1} = ft_timelockgrandaverage(cfg, timelock_tmp{:,1});
timelock_merged{2} = ft_timelockgrandaverage(cfg, timelock_tmp{:,2});
cfg.keepindividual = 'no';
cfg.parameter = 'individual';
timelock_merged{3} = ft_timelockgrandaverage(cfg, timelock_merged{1:2});
cfg=[];
cfg.operation = 'subtract';
cfg.parameter = 'individual';
timelock_merged{4} = ft_math(cfg, timelock_merged{1}, timelock_merged{2});
% add behavior
timelock_merged{1}.trialinfo = trialinfo(:,:,1);
timelock_merged{2}.trialinfo = trialinfo(:,:,2);
timelock_merged{3}.trialinfo = mean(trialinfo, 3);
timelock_merged{4}.trialinfo = trialinfo(:,:,1) - trialinfo(:,:,2);

%% plot TFR
cfg=[];
cfg.layout = lay;
%     cfg.zlim = [0.8 1.2];
cfg.colorbar = 'yes';
% cfg.baseline = [-0.5 0];
% cfg.baselinetype = 'relchange';
cfg.zlim = 'maxabs';
% cfg.zlim = [1.17 1.23];
%     cfg.xlim = [-0.5 1.5];
% ft_multiplotTFR(cfg,mse_merged{3})
ft_multiplotTFR(cfg,freq_merged{3})
%% plot time courses
cfg=[];
cfg.layout = lay;
cfg.frequency = [10 20 ];
% cfg.baseline = [-0.5 0];
% cfg.baselinetype = 'relchange';
% ft_multiplotER(cfg,mse_merged{4})
ft_multiplotER(cfg,freq_merged{1:2})

%% plot freq
cfg=[];
cfg.layout = lay;
cfg.colorbar = 'yes';
% cfg.baseline = [-0.5 0];
% cfg.baselinetype = 'relchange';
cfg.zlim = 'maxabs';
ft_multiplotTFR(cfg,freq_merged{2})
% % plot freq time courses
% cfg=[];
% cfg.layout = lay;
% cfg.frequency = 6;
% % cfg.baseline = [-0.5 0];
% % cfg.baselinetype = 'relchange';
% cfg.zlim = 'maxabs';
% ft_multiplotER(cfg,freq_merged{1})
%% run stats: correlation mMSE vs behavior
cfg = []; 
cfg.method    = 'triangulation';
% cfg.layout = 'acticap-64ch-standard2.mat';
cfg.layout = lay;
% cfg.template = 'elec1010_neighb.mat';
neighbours       = ft_prepare_neighbours(cfg);
% cfg = [];
% % cfg.layout = 'acticap-64ch-standard2.mat';
% % cfg.layout = 'EEG1010.lay';
% cfg.layout = lay;
% cfg.neighbours = neighbours;
% ft_neighbourplot(cfg)

cfg = []; 
cfg.design = freq_merged{4}.trialinfo(:,3); % 3 is accuracy
% cfg.frequency = 80;
cfg.latency = [-0.5 1.5];
cfg.uvar     = [];
cfg.ivar     = 1;
cfg.method           = 'montecarlo';
cfg.statistic        = 'ft_statfun_correlationT';  %depsamplesT ft_statfun_correlationT_corrcol
% cfg.statistic        = 'ft_statfun_partialcorrelationT';  %depsamplesT ft_statfun_correlationT_corrcol
cfg.type             = 'Pearson'; % Spearman Pearson
cfg.correctm         = 'cluster';  %'no'
cfg.clusteralpha     = 0.05;
cfg.clusterstatistic = 'maxsum';
cfg.tail             = 0;
cfg.clustertail      = 0;
cfg.alpha            = 0.025;
cfg.numrandomization = 100;
cfg.neighbours       = neighbours;
% cfg.minnbchan        = 0;
cfg.spmversion = 'spm12';
corrstat = ft_freqstatistics(cfg, mse_merged{4}); % incong - cong

%% plot 
cfg=[];
% cfg.layout = 'acticap-64ch-standard2.mat';
% cfg.layout = 'EEG1010.lay';
cfg.layout = lay;
cfg.parameter = 'rho';
%     cfg.zlim = [0.8 1.2];
cfg.colorbar = 'yes';
% cfg.baseline = [-0.5 0];
% cfg.baselinetype = 'relchange';
cfg.zlim = 'maxabs';
% cfg.zlim = [1.17 1.23];
%     cfg.xlim = [-0.5 1.5];
% cfg.maskparameter = 'mask';
% corrstat.mask = double(corrstat.posclusterslabelmat == 1);
ft_multiplotTFR(cfg, corrstat)
% ft_multiplotER(cfg, corrstat)

%% get data for scatter
cfg=[];
% cfg.latency = [0.25 0.4];
% cfg.frequency = [70.027211 142.952381];
cfg.latency = [0.5];
cfg.frequency = [50 100];
% cfg.channel = {'FC1' 'FC3' 'FCz'};
cfg.channel = {'FC3', 'FC1', 'FCz', 'FC2', 'C3', 'C1', 'Cz', 'C2'};
cfg.avgoverchan = 'yes';
cfg.avgovertime = 'yes';
cfg.avgoverfreq = 'yes';
corrdat = ft_selectdata(cfg, mse_merged{4});
figure; scatter(corrdat.powspctrm , mse_merged{4}.trialinfo(:,3), 100, 'MarkerEdgeColor',[1 1 1],...
  'MarkerFaceColor',[0 0 0], 'LineWidth',1.5);
xlabel('Incongruent–congruent mMSE');    ylabel('Incongruent–congruent Accuracy')
[rho, p] = corr(corrdat.powspctrm , mse_merged{4}.trialinfo(:,3));
title(sprintf('r = %.2f, p = %g', rho, p));
lsline; box on; axis square
