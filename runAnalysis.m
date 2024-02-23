if ismac
%   datapath = '/Users/kloosterman/Library/CloudStorage/Dropbox/PROJECTS/LRaudio/data';
%   datapath = '/Users/kloosterman/gridmaster2012/projectdata/LRaudio/data';
  datapath = '/Users/kloosterman/gridmaster2012/projectdata/LRaudio/source';
  toolspath = '/Users/kloosterman/Documents/GitHub/';
else
%   datapath = '/home/mpib/kloosterman/projectdata/LRaudio/data';
  datapath = '/home/mpib/kloosterman/projectdata/LRaudio/source';
  toolspath = '/home/mpib/kloosterman/GitHub/';
end

restoredefaultpath
addpath(fullfile(toolspath, 'fieldtrip')); ft_defaults
addpath(fullfile(toolspath, 'mMSE'));
addpath(fullfile(toolspath, 'LRaudio'))
addpath(fullfile(toolspath, 'qsub-tardis')) 

outputpath = fullfile(fileparts(datapath), 'outputdata');
mkdir(outputpath)

% make eeg layout
load('acticap-64ch-standard2.mat')
lay.label(find(contains(lay.label, 'Ref'))) = {'FCz'};

% define subjects
SUBJ={};
for i=1%:36
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
cfg.evoked = ''; % empty, regress, or subtract
cfg.csd = ''; % empty or csd
cfg.sensor_or_source = 'source';
overwrite = 1;
for isub = 1:length(SUBJ)
  cfg.SUBJ = SUBJ{isub};
  for icond = 2 %1:2
    cfg.icond = icond;
    switch cfg.sensor_or_source
      case 'sensor'
        cfg.datafile = fullfile(datapath, SUBJ{isub}, sprintf('clean_SUB%s', [SUBJ{isub} '.mat']));
        cfg.outpath = fullfile(fileparts(datapath), cfg.analysis, cfg.evoked, cfg.csd, sprintf('SUB%s_cond%d.mat', SUBJ{isub}, icond));
      case 'source'
        cfg.datafile = fullfile(datapath, SUBJ{isub}, sprintf('SourceTimeSeries_BW_1-100Hz_ParcelSpace_Block*.mat'));
        cfg.outpath = fullfile(datapath, cfg.analysis, cfg.evoked, cfg.csd, sprintf('SUB%s_cond%d.mat', SUBJ{isub}, icond));
    end    
    if overwrite || ~exist(cfg.outpath, 'file')
      mkdir(fileparts(cfg.outpath))
      cfglist{end+1} = cfg;
    else
      disp('File exists, skipping')
    end
  end
end

if ismac % submit to tardis, or run locally
  cellfun(@computemMSE, cfglist, 'Uni', 0);
else % mse: 'memreq', 100e9, 'timreq', 23*60*60, 'options', ' --cpus-per-task=4 '
  if strcmp(cfg.analysis, 'mse')
    qsubcellfun(@computemMSE, cfglist, 'memreq', 100e9, 'timreq', 46*60*60, 'stack', 1, ...
      'StopOnError', false, 'backend', 'slurm', 'options', ' --cpus-per-task=4 --partition long');
  else
    qsubcellfun(@computemMSE, cfglist, 'memreq', 5e9, 'timreq', 1*60*60, 'stack', 1, ...
      'StopOnError', false, 'backend', 'slurm', 'options', ' --cpus-per-task=4 ');
  end
  return
end

%% Merge mse files across subjects and cond
disp 'Merge analyses'
mse_tmp = {}; freq_tmp = {}; timelock_tmp = {};
trialinfo = [];
for isub = 1:length(SUBJ)
%   datafile = fullfile(datapath, SUBJ{isub}, sprintf('clean_SUB%s', [SUBJ{isub} '.mat']));
  for icond = 1:2
%     path = fullfile(fileparts(datapath), 'mse', 'regress', 'csd', sprintf('SUB%s_cond%d.mat', SUBJ{isub}, icond));
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
%     path = fullfile(fileparts(datapath), 'freq', 'regress', 'csd', sprintf('SUB%s_cond%d.mat', SUBJ{isub}, icond));
    path = fullfile(fileparts(datapath), 'freq', 'regress', sprintf('SUB%s_cond%d.mat', SUBJ{isub}, icond));
    disp(path)
    if exist(path, 'file')
      load(path)
      freq_tmp{isub,icond} = freq; % freq_bl or raw: freq
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
% combine mse 
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

% combine freq 
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

% combine timelock 
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
cfg=[];   cfg.layout = lay;
cfg.zlim = 'maxabs';
% cfg.zlim = [1.17 1.23];
% cfg.xlim = [-0.5 1.5];
% cfg.baseline = [-0.5 0];    cfg.baselinetype = 'relchange';
ft_multiplotTFR(cfg, freq_merged{4}); colorbar
%% plot time courses
cfg=[];   cfg.layout = lay;
cfg.frequency = [4 8 ];
ft_multiplotER(cfg,freq_merged{1:2})

%% plot Central pooling mse 
xlim = [0.75 1.25];
channel = {'FC3', 'FC1', 'FCz', 'FC2', 'C3', 'C1', 'Cz', 'C2'};
% channel = {'FC3', 'FC1', 'FCz', 'FC2', 'C2', 'Cz', 'C1', 'C3', 'CP3', 'CP1', 'CPz', 'CP2'};
% channel = {'C3', 'CP3', 'P3'};

f = figure; f.Position = [680         520        800         800*0.5]; 
cfg=[];   cfg.layout = lay;   cfg.figure = 'gcf';
cfg.channel = channel;
cfg.zlim = [1.17 1.23];
% cfg.zlim = [-0.15 0.15]; 
subplot(2,3,1);     ft_singleplotTFR(cfg, mse_merged{3});

cfg=[];   cfg.layout = lay;   cfg.figure = 'gcf';
cfg.xlim = xlim;   cfg.ylim = [50 100];   
% cfg.zlim = [1.17 1.23]; cfg.highlightchannel = channel;
subplot(2,3,2);     ft_topoplotTFR(cfg, mse_merged{3}); colorbar

cfg=[];    cfg.figure = 'gcf';  cfg.channel = channel;
subplot(2,3,3);     ft_singleplotER(cfg, mse_merged{1:2});

% plot Central pooling freq 
cfg=[];   cfg.layout = lay;   cfg.figure = 'gcf';   cfg.channel = channel;
cfg.zlim = [-0.15 0.15]; %'maxabs';      
subplot(2,3,4);     ft_singleplotTFR(cfg, freq_merged{3});

cfg=[];   cfg.layout = lay;   cfg.figure = 'gcf';
cfg.xlim = xlim;   cfg.ylim = [2 8];    cfg.zlim = 'maxabs';
subplot(2,3,5);     ft_topoplotTFR(cfg, freq_merged{3}); colorbar

cfg=[];    cfg.figure = 'gcf';          cfg.channel = channel;
subplot(2,3,6);     ft_singleplotER(cfg, freq_merged{1:2}); shg

plotpath = '/Users/kloosterman/Library/CloudStorage/Dropbox/PROJECTS/LRaudio/plots';
orient landscape
saveas(f, fullfile(plotpath, ['msevsfreq_' [channel{:}]]), 'pdf')

%% run stats: correlation mMSE vs behavior
cfg = []; 
cfg.method    = 'triangulation';
cfg.layout = lay;
neighbours       = ft_prepare_neighbours(cfg);
% cfg.neighbours = neighbours;
% ft_neighbourplot(cfg)

cfg = []; 
cfg.design = mse_merged{4}.trialinfo(:,3); % 3 is accuracy
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

%% plot correlation
cfg=[];    cfg.layout = lay;
cfg.parameter = 'rho';
cfg.colorbar = 'yes';
cfg.zlim = 'maxabs';
% cfg.zlim = [1.17 1.23];
%     cfg.xlim = [-0.5 1.5];
% cfg.maskparameter = 'mask';
% corrstat.mask = double(corrstat.posclusterslabelmat == 1);
ft_multiplotTFR(cfg, corrstat)
% ft_multiplotER(cfg, corrstat)

%% get data for scatter
cfg=[];
cfg.latency = [0.1 0.4];
% cfg.latency = [0.75 1.25];
cfg.frequency = [70.027211 142.952381];
% cfg.latency = [0.5];
% cfg.frequency = [2 8];
% cfg.channel = {'C3', 'CP3', 'P3'};
cfg.channel = {'FC3', 'FC1', 'FCz', 'FC2', 'C3', 'C1', 'Cz', 'C2'};
% cfg.channel = {'FC3', 'FC1', 'FCz', 'FC2', 'C2', 'Cz', 'C1', 'C3', 'CP3', 'CP1', 'CPz', 'CP2'};
cfg.avgoverchan = 'yes'; cfg.avgovertime = 'yes'; cfg.avgoverfreq = 'yes';
corrdat = ft_selectdata(cfg, mse_merged{4});
figure; scatter(corrdat.powspctrm , mse_merged{4}.trialinfo(:,3), 100, 'MarkerEdgeColor',[1 1 1],...
  'MarkerFaceColor',[0 0 0], 'LineWidth',1.5);
xlabel('Incongruent–congruent mMSE');    ylabel('Incongruent–congruent Accuracy')
[rho, p] = corr(corrdat.powspctrm , mse_merged{4}.trialinfo(:,3), 'type', 'Spearman');
title(sprintf('r = %.2f, p = %g', rho, p));
lsline; box on; axis square
