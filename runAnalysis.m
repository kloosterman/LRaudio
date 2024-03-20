if ismac
%   datapath = '/Users/kloosterman/Library/CloudStorage/Dropbox/PROJECTS/LRaudio/data';
  datapath = '/Users/kloosterman/gridmaster2012/projectdata/LRaudio/data';
%   datapath = '/Users/kloosterman/gridmaster2012/projectdata/LRaudio/source';
  toolspath = '/Users/kloosterman/Documents/GitHub/';
else
  datapath = '/home/mpib/kloosterman/projectdata/LRaudio/data';
%   datapath = '/home/mpib/kloosterman/projectdata/LRaudio/source';
  toolspath = '/home/mpib/kloosterman/GitHub/';
end

restoredefaultpath
addpath(fullfile(toolspath, 'fieldtrip')); ft_defaults
addpath(fullfile(toolspath, 'mMSE'));
addpath(fullfile(toolspath, 'LRaudio'))
addpath(fullfile(toolspath, 'qsub-tardis'))
addpath(genpath(fullfile(toolspath, 'BrainSlicer')))

outputpath = fullfile(fileparts(datapath), 'outputdata');
mkdir(outputpath)

% make eeg layout
load('acticap-64ch-standard2.mat')
lay.label(find(contains(lay.label, 'Ref'))) = {'FCz'};

% define subjects
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
cfg.analysis = 'freq'; % freq, mse, or erp
cfg.evoked = 'subtract'; % empty, regress, or subtract
cfg.csd = ''; % empty or csd
cfg.sensor_or_source = 'source';
overwrite = 1;
for isub = 1:length(SUBJ)
  cfg.SUBJ = SUBJ{isub};
  for icond = 1:2
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
    path = fullfile(fileparts(datapath), 'mse', 'subtract', sprintf('SUB%s_cond%d.mat', SUBJ{isub}, icond));
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
    path = fullfile(fileparts(datapath), 'freq', 'subtract', sprintf('SUB%s_cond%d.mat', SUBJ{isub}, icond));
    disp(path)
    if exist(path, 'file')
      load(path)
      freq_tmp{isub,icond} = freq_bl; % freq_bl (baseline corrected) freq (raw power)
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
% close all
xlim = [0.75 1.25];
channel = {'FC5', 'FC3', 'FC1', 'FCz', 'FC2', 'FC4', 'C2', 'Cz', 'C1', 'C3', 'F3', 'F1'};
% channel = {'FC3', 'FC1', 'FCz', 'FC2', 'C3', 'C1', 'Cz', 'C2'};
% channel = {'FC3', 'FC1', 'FCz', 'FC2', 'C2', 'Cz', 'C1', 'C3', 'CP3', 'CP1', 'CPz', 'CP2'};
% channel = {'C3', 'CP3', 'P3'};

f = figure; f.Position = [680         520        800         800*0.5]; 
cfg=[];   cfg.layout = lay;   cfg.figure = 'gcf';
cfg.channel = channel;  cfg.zlim = [1.17 1.23]; cfg.title = 'Entropy';
subplot(2,3,1);     ft_singleplotTFR(cfg, mse_merged{3}); xline(0); xlabel('Time from stim (s)'); ylabel('Time scale (ms)')

cfg=[];   cfg.layout = lay;   cfg.figure = 'gcf';
cfg.xlim = [0.1 0.5];   cfg.ylim = [50 100];   % cfg.zlim = [1.17 1.23]; 
cfg.highlightchannel = channel; cfg.highlight = 'on';
subplot(2,3,2);     ft_topoplotTFR(cfg, mse_merged{3}); %colorbar

cfg=[];    cfg.figure = 'gcf';  cfg.channel = channel; cfg.title = ' ';
subplot(2,3,3);     ft_singleplotER(cfg, mse_merged{1:2}); hold on; xline(0,'HandleVisibility','off'); 
legend({'Incong.', 'Cong.'}, 'Location', 'South'); legend boxoff
xlabel('Time from stim (s)'); ylabel('Entropy')
% ax=gca; YL = ax.YLim; patch([0; 0; 0.46; 0.46], [YL(1); YL(2); YL(2); YL(1);], [0.5 0.5 0.5 ])

% plot Central pooling freq 
cfg=[];   cfg.layout = lay;   cfg.figure = 'gcf';   cfg.channel = channel; 
cfg.ylim = [0 20]; cfg.zlim = [-0.2 0.2]; cfg.title = 'Power';
subplot(2,3,4); ft_singleplotTFR(cfg, freq_merged{3}); xline(0); xlabel('Time from stim (s)'); ylabel('Frequency (Hz)')

cfg=[]; cfg.layout = lay;   cfg.figure = 'gcf'; %cfg.colormap = colormaps(2);
cfg.xlim = [0.1 0.5];   cfg.ylim = [2 8];    cfg.zlim = 'maxabs';
subplot(2,3,5); ft_topoplotTFR(cfg, freq_merged{3}); %colorbar;

cfg=[]; cfg.figure = 'gcf'; cfg.channel = channel; cfg.title = ' ';  cfg.frequency = [2 8]; % cfg.showlegend = 'yes'; cfg.dataname = {'Incong.', 'Cong.'};
subplot(2,3,6); cfg = ft_singleplotER(cfg, freq_merged{1:2}); 
xline(0,'HandleVisibility','off'); % legend({'Incong.', 'Cong.'}, 'Location', 'Southeast'); legend boxoff
xlabel('Time from stim (s)'); ylabel('Power (psc)'); % %shg

plotpath = '/Users/kloosterman/Library/CloudStorage/Dropbox/PROJECTS/LRaudio/plots';
orient landscape
saveas(f, fullfile(plotpath, ['msevsfreq_' [channel{:}]]), 'pdf')
saveas(f, fullfile(plotpath, ['msevsfreq_' [channel{:}]]), 'png')

%% run stats: correlation mMSE vs behavior
cfg = []; cfg.method = 'triangulation'; cfg.layout = lay;
neighbours = ft_prepare_neighbours(cfg); % cfg.neighbours = neighbours; ft_neighbourplot(cfg)

% get 2-8 Hz power to control for
cfg=[]; cfg.frequency = [2 8]; cfg.latency = [0.1 0.5]; cfg.channel = channel; cfg.avgoverchan = 'yes'; cfg.avgovertime = 'yes'; cfg.avgoverfreq = 'yes';
freq_control = ft_selectdata(cfg, freq_merged{4});

cfg=[]; cfg.design = mse_merged{4}.trialinfo(:,3); % 3 is accuracy
cfg.frequency = [60 100]; cfg.avgoverfreq = 'yes';
cfg.channel = channel;    cfg.avgoverchan = 'yes';
cfg.latency = [-0.5 1.5];
cfg.statistic = 'ft_statfun_partialcorrelationT';  %ft_statfun_correlationT depsamplesT ft_statfun_correlationT_corrcol
cfg.type = 'Pearson'; % Spearman Pearson
cfg.correctm = 'cluster';  %'no'
cfg.numrandomization = 1000; cfg.uvar=[]; cfg.ivar = 1; cfg.method = 'montecarlo'; cfg.clusteralpha = 0.05; cfg.clusterstatistic = 'maxsum'; cfg.tail = 0; cfg.clustertail = 0; cfg.alpha = 0.025; cfg.spmversion = 'spm12'; cfg.neighbours       = neighbours;
corrstat_mse = ft_freqstatistics(cfg, mse_merged{4}); % incong - cong

cfg.design = [mse_merged{4}.trialinfo(:,3) freq_control.powspctrm]; % control for theta
corrstat_mse_control = ft_freqstatistics(cfg, mse_merged{4}); % incong - cong

cfg.frequency = [2 8]; % straight corr theta vs behav
cfg.design = mse_merged{4}.trialinfo(:,3);
corrstat_freq = ft_freqstatistics(cfg, freq_merged{4}); % incong - cong

%% plot corr time series and scatter for significant part
clear pl; f=figure; f.Position = [  350   500   814   410];
subplot(1,2,1);hold on; pl(1) = plot(corrstat_mse.time, squeeze(corrstat_mse.rho)); title('Correlation time courses')
plot_sig_bar(corrstat_mse.time, squeeze(corrstat_mse.mask)', -0.25, 5, pl(1).Color)

pl(2)=plot(squeeze(corrstat_mse_control.time), squeeze(corrstat_mse_control.rho)); 
plot_sig_bar(corrstat_mse_control.time, squeeze(corrstat_mse_control.mask)', -0.35, 5, pl(2).Color)

pl(3)=plot(squeeze(corrstat_freq.time), squeeze(corrstat_freq.rho)); xline(0); yline(0);
plot_sig_bar(corrstat_freq.time, squeeze(corrstat_freq.mask)', -0.45, 5, pl(3).Color);
text(0.4, -0.4, 'p<0.05, corrected')
legend(pl, {'Entropy v. behav' 'Entropy v. behav controlled for 2-8Hz' '2-8 Hz power v. behav'}, 'Location', 'southoutside')
xlabel('Time from stim (s)'); ylabel('Pearson correlation')

cfg=[]; cfg.channel = channel; cfg.latency = [0.1 0.5]; cfg.frequency = [60 100]; %[70.027211 142.952381];
cfg.avgoverchan = 'yes'; cfg.avgovertime = 'yes'; cfg.avgoverfreq = 'yes';
corrdat = ft_selectdata(cfg, mse_merged{4});
subplot(1,2,2); scatter(corrdat.powspctrm , mse_merged{4}.trialinfo(:,3), 100, 'MarkerEdgeColor',[1 1 1],...
  'MarkerFaceColor',[0 0 0], 'LineWidth',1.5);
xlabel('Incongruent–congruent mMSE');    ylabel('Incongruent–congruent Accuracy')
[rho, p] = corr(corrdat.powspctrm , mse_merged{4}.trialinfo(:,3), 'type', 'Pearson');
% control for theta
[rho_p, p_p] = partialcorr(corrdat.powspctrm , mse_merged{4}.trialinfo(:,3), freq_control.powspctrm, 'type', 'Pearson');
title(sprintf('r = %.2f, p = %1.3f\ncontrolled for theta: r = %.2f, p = %1.3f', rho, p, rho_p, p_p));
lsline; box on; axis square

orient landscape
saveas(f, fullfile(plotpath, 'Corr_tc_scatter'), 'pdf')
saveas(f, fullfile(plotpath, 'Corr_tc_scatter'), 'png')


% 
% %% plot correlation
% cfg=[];    cfg.layout = lay;
% cfg.parameter = 'rho';
% cfg.colorbar = 'yes';
% cfg.zlim = 'maxabs';
% % cfg.zlim = [1.17 1.23];
% %     cfg.xlim = [-0.5 1.5];
% cfg.maskparameter = 'mask';
% corrstat.mask = double(corrstat.posclusterslabelmat == 1);
% ft_multiplotTFR(cfg, corrstat)
% % ft_multiplotER(cfg, corrstat)
% 
% orient landscape
% saveas(f, fullfile(plotpath, 'Corr_tc_scatter'), 'pdf')
% 
