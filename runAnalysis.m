restoredefaultpath
global datapath
if ismac
  datapath = '/Users/kloosterman/projectdata/LRaudio';
  toolspath = '/Users/kloosterman/Documents/GitHub/';
  addpath(fullfile(toolspath, 'fieldtrip-20231025')); ft_defaults
else
  datapath = '/home/mpib/kloosterman/projectdata/LRaudio';
  %   datapath = '/home/mpib/kloosterman/projectdata/LRaudio/source';
  toolspath = '/home/mpib/kloosterman/GitHub/';
  addpath(fullfile(toolspath, 'fieldtrip')); ft_defaults
end

addpath(fullfile(toolspath, 'mMSE'));
addpath(fullfile(toolspath, 'LRaudio'))
addpath(fullfile(toolspath, 'qsub-tardis'))
addpath(genpath(fullfile(toolspath, 'BrainSlicer')))
% plotpath = '/Users/kloosterman/Library/CloudStorage/Dropbox/PROJECTS/LRaudio/plots';
plotpath = '/Users/kloosterman/Dropbox/PROJECTS/LRaudio/plots';

% make eeg layout
% load('acticap-64ch-standard2.mat')
% lay.label(find(contains(lay.label, 'Ref'))) = {'FCz'};
load Acticap_64_UzL.mat

% define subjects
SUBJ={};
for i=1:36
  SUBJ{end+1} = sprintf('%d', i);
end
SUBJbool = true(size(SUBJ));
% SUBJbool([10, 12, 15, 17 ]) = false; % exclude 10 and 15, 17 missing, 12 conditions missing
% exclude 10 and 15, 17 missing, 12 conditions missing,
% #6 staircase not converging + datamissing, accuracy only 0.65
% #31: staircase at max, accuracy only 0.6006
% #3: staircase at minimum (0.01) in many trials
% 34 sc at max until halfway, then drops, convergence weird?
SUBJbool([10, 12, 15, 17,    31, 6, 3 ]) = false;
% SUBJbool([10, 12, 15, 17 ]) = false;
SUBJ = SUBJ(SUBJbool);
disp(SUBJ)

%% compute entropy, freqanalysis and ERPs
% make list of files to analyze on tardis
overwrite = 1;
cfglist = {}; cfg=[];
cfg.evoked = ''; % empty, regress, or subtract
cfg.csd = ''; % empty or csd
cfg.sensor_or_source = 'sensor';
cfg.runperblock = 'yes'; % empty for all together, or per block.
for isub = 1:length(SUBJ)
  cfg.SUBJ = SUBJ{isub};
  for icond = 1:2
    cfg.icond = icond;
    switch cfg.sensor_or_source
      case 'sensor'
        cfg.datafile = fullfile(datapath, 'data', SUBJ{isub}, sprintf('clean_SUB%s', [SUBJ{isub} '.mat']));
      case 'source'
        cfg.datafile = fullfile(datapath, 'source', SUBJ{isub}, sprintf('SourceTimeSeries_BW_1-100Hz_ParcelSpace_Block*.mat'));
    end
    if overwrite;      cfglist{end+1} = cfg;       end
  end
end

if ismac % submit to tardis, or run locally
  cellfun(@computemMSE, cfglist, 'Uni', 0);
else % mse: 'memreq', 100e9, 'timreq', 23*60*60, 'options', ' --cpus-per-task=4 '
  %   if strcmp(cfg.analysis, 'mse')
  % all scales:
  %     qsubcellfun(@computemMSE, cfglist, 'memreq', 100e9, 'timreq', 46*60*60, 'stack', 1, ...
  %       'StopOnError', false, 'backend', 'slurm', 'options', ' --cpus-per-task=4 --partition long');
  % scales 10-40
  qsubcellfun(@computemMSE, cfglist, 'memreq', 10e9, 'timreq', 8*60*60, 'stack', 1, ...
    'StopOnError', false, 'backend', 'slurm', 'options', ' --cpus-per-task=4 --partition long');
  %   else
  %     qsubcellfun(@computemMSE, cfglist, 'memreq', 5e9, 'timreq', 1*60*60, 'stack', 1, ...
  %       'StopOnError', false, 'backend', 'slurm', 'options', ' --cpus-per-task=4 ');
  % end
  return
end

%% Merge mse files across subjects and cond
evoked = 'subtract'; % subtract_avgref subtract
csd = ''; % csd
runperblock = 'yes';

disp 'Merge analyses'
mse_tmp = {}; freq_tmp = {}; timelock_tmp = {};
trialinfo = []; trialinfo_trl = {}; trialinfo_blocks = {};
for isub = 1:length(SUBJ)
  for icond = 1:2
    path = fullfile(datapath, 'mse', evoked, csd, sprintf('SUB%s_cond%d.mat', SUBJ{isub}, icond));
    disp(path);      load(path)
    mse.freq = mse.timescales;
    mse.powspctrm = mse.sampen;
    mse.dimord = 'chan_freq_time';
    mse_tmp{isub,icond} = mse;
    trialinfo(isub,:,icond) = mean(mse.trialinfo);
    trialinfo_trl{isub,icond} = mse.trialinfo;

    path = fullfile(datapath, 'freq', evoked, csd, sprintf('SUB%s_cond%d.mat', SUBJ{isub}, icond));
    disp(path);      load(path)
    freq_tmp{isub,icond} = freq_bl; % freq_bl (baseline corrected) freq (raw power)

    path = fullfile(datapath, 'timelock', '', csd, sprintf('SUB%s_cond%d.mat', SUBJ{isub}, icond));
    disp(path);      load(path)
    timelock_tmp{isub,icond} = timelock;

    if strcmp(runperblock, 'yes')
      path = fullfile(datapath, 'mse', evoked, csd, sprintf('SUB%s_cond%d_blocks.mat', SUBJ{isub}, icond));
      disp(path);    load(path)
      for iblock=1:length(mse)
        mse{iblock}.freq = mse{iblock}.timescales;
        mse{iblock}.powspctrm = mse{iblock}.sampen;
        mse{iblock}.dimord = 'chan_freq_time';
        mse{iblock}.trialinfo = mean(mse{iblock}.trialinfo);
      end
      cfg=[];
      cfg.appenddim = 'rpt';
      mse_blocks{isub,icond} = ft_appendfreq(cfg, mse{:}); % ft_appendfreq?
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
% cfg.keepindividual = 'no';
% cfg.parameter = 'individual';
% timelock_merged{3} = ft_timelockgrandaverage(cfg, timelock_merged{1:2});
cfg=[];
cfg.operation = 'subtract';
cfg.parameter = 'individual';
timelock_merged{4} = ft_math(cfg, timelock_merged{1}, timelock_merged{2});
% add behavior
timelock_merged{1}.trialinfo = trialinfo(:,:,1);
timelock_merged{2}.trialinfo = trialinfo(:,:,2);
timelock_merged{3}.trialinfo = mean(trialinfo, 3);
timelock_merged{4}.trialinfo = trialinfo(:,:,1) - trialinfo(:,:,2);

%% plot ERPs
cfg=[];   cfg.layout = lay;
% cfg.baseline = [-0.5 0];    cfg.baselinetype = 'relchange';
ft_multiplotER(cfg, timelock_merged{1}, timelock_merged{2}); colorbar
%% plot TFR / topo MSE
cfg=[];   cfg.layout = lay;
% cfg.zlim = 'maxabs';
cfg.zlim = [1.17 1.23];
% cfg.xlim = [-0.5 1.5];
% cfg.baseline = [-0.5 0];    cfg.baselinetype = 'relchange';
ft_multiplotTFR(cfg, mse_merged{3}); colorbar
%% plot TFR / topo freq
cfg=[];   cfg.layout = lay; cfg.zlim = 'maxabs'; % cfg.xlim = [-0.5 1.5];
ft_multiplotTFR(cfg, freq_merged{3}); colorbar
%% plot time courses
cfg=[];   cfg.layout = lay;
cfg.frequency = [4 8 ];
ft_multiplotER(cfg,freq_merged{1:2})

%% plot Central pooling mse
close all
XLIM = [-0.5 1.5];
zlim = [1.15 1.23];
% channel = {'FC5', 'FC3', 'FC1', 'FCz', 'FC2', 'FC4', 'C2', 'Cz', 'C1', 'C3', 'F3', 'F1'};
% channel = {'FC2', 'FCz', 'FC1', 'C3', 'C1', 'Cz', 'C2', 'CP1'}; % CSD
% channel = {'FC2', 'FCz', 'FC1', 'C1', 'Cz', 'C2'}; % avg ref 'C3',
channel = {'FC2', 'FCz', 'FC1', 'C3', 'C1', 'Cz', 'C2', 'CP2', 'CPz', 'CP1'}
channel = {'FCz', 'FC1', 'C1', 'Cz', 'CPz', 'CP1'}
% channel = {'FC2', 'FCz', 'FC1', 'C3', 'C1', 'Cz', 'C2', 'CP2', 'CPz', 'CP1'} % for csd
channel = {'FC2', 'FCz', 'FC1',  'C1', 'Cz', 'C2'} % for avgref  'C3', , 'CP2', 'CPz', 'CP1'

f = figure; f.Position = [680         520        800         800*0.5];
cfg=[];   cfg.layout = lay;   cfg.figure = 'gcf'; cfg.channel = channel;  cfg.zlim = zlim;  cfg.xlim = XLIM; cfg.title = 'Entropy';
subplot(2,3,1);     ft_singleplotTFR(cfg, mse_merged{3}); xline(0); xlabel('Time from stim (s)'); ylabel('Time scale (ms)')

cfg=[];   cfg.layout = lay;   cfg.figure = 'gcf'; cfg.xlim = [0 0.3];   cfg.ylim = [70 90];   cfg.zlim = zlim; cfg.highlightchannel = channel; cfg.highlight = 'on';%cfg.zlim = [1.17 1.22];
subplot(2,3,2);     ft_topoplotTFR(cfg, mse_merged{3}); %colorbar

cfg=[];    cfg.figure = 'gcf';  cfg.channel = channel; cfg.title = ' ';  cfg.xlim = XLIM;
subplot(2,3,3);     ft_singleplotER(cfg, mse_merged{1:2}); hold on; xline(0,'HandleVisibility','off');
legend({'Incong.', 'Cong.'}, 'Location', 'South'); legend boxoff
xlabel('Time from stim (s)'); ylabel('Entropy')
% ax=gca; YL = ax.YLim; patch([0; 0; 0.46; 0.46], [YL(1); YL(2); YL(2); YL(1);], [0.5 0.5 0.5 ])

% plot Central pooling freq
cfg=[];   cfg.layout = lay;   cfg.figure = 'gcf';   cfg.channel = channel;
cfg.xlim = XLIM; cfg.ylim = [0 100]; cfg.zlim = [-0.17 0.17]; cfg.title = 'Power modulation';
subplot(2,3,4); ft_singleplotTFR(cfg, freq_merged{3}); xline(0); xlabel('Time from stim (s)'); ylabel('Frequency (Hz)')

channel_freq = {'FCz',  'Cz', 'CPz', 'CP1'}; %   'CP3' 'C1',
% cfg=[]; cfg.layout = lay;   cfg.figure = 'gcf'; cfg.highlightchannel = channel_freq; cfg.highlight = 'on';
cfg=[]; cfg.layout = lay;   cfg.figure = 'gcf'; cfg.highlightchannel = channel; cfg.highlight = 'on';
cfg.xlim = [0 0.2];   cfg.ylim = [4 8];    cfg.zlim = 'maxabs';
% cfg.xlim = [0.1 1.5];   cfg.ylim = [60 80];    cfg.zlim = 'maxabs';
subplot(2,3,5); ft_topoplotTFR(cfg, freq_merged{3}); %colorbar;

cfg=[]; cfg.figure = 'gcf'; cfg.channel = channel; cfg.title = ' ';  cfg.frequency = [4 8]; cfg.xlim = XLIM;
% cfg=[]; cfg.figure = 'gcf'; cfg.channel = channel_freq; cfg.title = ' ';  cfg.frequency = [60 80];  cfg.xlim = XLIM;
subplot(2,3,6); cfg = ft_singleplotER(cfg, freq_merged{1:2});
xline(0,'HandleVisibility','off'); % legend({'Incong.', 'Cong.'}, 'Location', 'Southeast'); legend boxoff
xlabel('Time from stim (s)'); ylabel('Power (% signal change)'); % %shg

orient landscape
saveas(f, fullfile(plotpath, ['msevsfreq_' [channel{:}]]), 'pdf')
saveas(f, fullfile(plotpath, ['msevsfreq_' [channel{:}]]), 'png')

%% run stats: correlation mMSE vs behavior
behav_col=3;% 3 is accuracy,6 is delta_fsemitones
corrtype='Pearson'; % Spearman Pearson
% corrtype='Spearman'; % Spearman Pearson
cond_leg = {'Incong.', 'Congr.', 'Cong. avg', 'Incong–Congr.'};
colors = {'r' 'b' 'g'};
ctrl_band = [4 8]; % 4 8 1 3
cfg = []; cfg.method = 'triangulation'; cfg.layout = lay;
neighbours = ft_prepare_neighbours(cfg); % cfg.neighbours = neighbours; ft_neighbourplot(cfg)

clear pl; f=figure; f.Position = [  350   250   400   400];
for icond = 4%1:4
  % get 2-8 Hz power to control for
  cfg=[]; cfg.frequency = ctrl_band; cfg.latency = [-0.1 0.1]; cfg.channel = channel; cfg.avgoverchan = 'yes'; cfg.avgovertime = 'yes'; cfg.avgoverfreq = 'yes';
  freq_control = ft_selectdata(cfg, freq_merged{icond});

  cfg=[]; cfg.design = mse_merged{icond}.trialinfo(:,behav_col); %
  cfg.frequency = [60 100]; cfg.avgoverfreq = 'yes';
  cfg.channel = channel;    cfg.avgoverchan = 'yes';
  cfg.latency = [-0.2 1.0];
  cfg.statistic = 'ft_statfun_partialcorrelationT';  %ft_statfun_correlationT depsamplesT ft_statfun_correlationT_corrcol
  cfg.type = corrtype;
  cfg.correctm = 'cluster';  %'no'
  cfg.numrandomization = 10000; cfg.uvar=[]; cfg.ivar = 1; cfg.method = 'montecarlo'; cfg.clusteralpha = 0.05; cfg.clusterstatistic = 'maxsum'; cfg.tail = 0; cfg.clustertail = 0;  cfg.spmversion = 'spm8'; cfg.neighbours       = neighbours;
  cfg.alpha = 0.025;
  corrstat_mse = ft_freqstatistics(cfg, mse_merged{icond}); % incong - cong

  cfg.design = [mse_merged{icond}.trialinfo(:,behav_col) freq_control.powspctrm]; % control for theta
  corrstat_mse_control = ft_freqstatistics(cfg, mse_merged{icond}); % incong - cong

  cfg.frequency = ctrl_band; % straight corr theta vs behav
  cfg.design = mse_merged{icond}.trialinfo(:,behav_col);
  corrstat_freq = ft_freqstatistics(cfg, freq_merged{icond}); % incong - cong

  % plot corr time series and scatter for significant part
  subplot(1,1,1);hold on;
  pl(1) = plot(corrstat_mse.time, squeeze(corrstat_mse.rho));  %'Color', colors{1}
  title(sprintf('%s\nCorrelation time courses', cond_leg{icond}))
  plot_sig_bar(corrstat_mse.time, squeeze(corrstat_mse.mask)', -0.25, 5, pl(1).Color)
  xline(0,'HandleVisibility','off'); yline(0,'HandleVisibility','off');
end
text(0.3, -0.3, 'p<0.05, corrected')
legend( cond_leg)
xlabel('Time from stim (s)'); ylabel([corrtype ' correlation'])
orient landscape
saveas(f, fullfile(plotpath, sprintf('Corr_percond%s.pdf', cond_leg{icond})), 'pdf')
saveas(f, fullfile(plotpath, sprintf('Corr_percond%s.png', cond_leg{icond})), 'png')

% % plot Incong-Cong inc controls + scatter
icond=4;
f=figure; f.Position = [  350   250   800   400]; clear pl
subplot(1,2,1);hold on;
pl(1) = plot(corrstat_mse.time, squeeze(corrstat_mse.rho),'Color', colors{1});  %
title(sprintf('%s\nCorrelation time courses', cond_leg{icond}))
plot_sig_bar(corrstat_mse.time, squeeze(corrstat_mse.mask)', -0.25, 5, pl(1).Color)
xline(0,'HandleVisibility','off'); yline(0,'HandleVisibility','off');
pl(2)=plot(squeeze(corrstat_mse_control.time), squeeze(corrstat_mse_control.rho), 'Color', colors{2});
plot_sig_bar(corrstat_mse_control.time, squeeze(corrstat_mse_control.mask)', -0.35, 5, pl(2).Color)
pl(3)=plot(squeeze(corrstat_freq.time), squeeze(corrstat_freq.rho), 'Color', colors{3}); xline(0); yline(0);
plot_sig_bar(corrstat_freq.time, squeeze(corrstat_freq.mask)', -0.45, 5, pl(3).Color);
text(0.3, -0.3, 'p<0.05, corrected')
legend(pl, {'Entropy v. behav' 'Entropy v. behav controlled for 4-8Hz' '4-8 Hz power v. behav'}, 'Location', 'southoutside')
xlabel('Time from stim (s)'); ylabel([corrtype ' correlation'])
% xlim(cfg.latency)

% scatter
cfg=[]; cfg.channel = channel; cfg.latency = [0.1 0.3]; cfg.frequency = [60 100]; %[70.027211 142.952381];
% cfg=[]; cfg.channel = channel; cfg.latency = [0.5 0.5]; cfg.frequency = [60 100]; %[70.027211 142.952381];
% cfg=[]; cfg.channel = channel; cfg.latency = [1 1.2]; cfg.frequency = [60 100]; %[70.027211 142.952381];
cfg.avgoverchan = 'yes'; cfg.avgovertime = 'yes'; cfg.avgoverfreq = 'yes';
corrdat = ft_selectdata(cfg, mse_merged{icond});
subplot(1,2,2);
scatter(corrdat.powspctrm , mse_merged{icond}.trialinfo(:,behav_col), 100, 'MarkerEdgeColor',[1 1 1],...
  'MarkerFaceColor',[0 0 0], 'LineWidth',1.5);
% scatter(corrdat.powspctrm , mse_merged{icond}.trialinfo(:,behav_col), 100, 'MarkerEdgeColor',[1 1 1],...
%   'MarkerFaceColor',[0 0 0], 'LineWidth',1.5, 'Marker', '.');
% text(corrdat.powspctrm , mse_merged{icond}.trialinfo(:,behav_col), SUBJ)
xlabel('Incongruent–congruent mMSE');    ylabel('Incongruent–congruent Accuracy')
[rho, p] = corr(corrdat.powspctrm , mse_merged{icond}.trialinfo(:,behav_col), 'type', corrtype);
% control for theta
[rho_p, p_p] = partialcorr(corrdat.powspctrm , mse_merged{icond}.trialinfo(:,behav_col), freq_control.powspctrm, 'type', corrtype);
% title(sprintf('r = %.2f, p = %1.3f\ncontrolled for theta: r = %.2f, p = %1.3f', rho, p, rho_p, p_p));
title(sprintf('%g-%g s, r = %.2f p = %.3f\ncontrolled for theta: r = %.2f p = %.3f',  cfg.latency, rho, p, rho_p, p_p));
lsline; box on; axis square

saveas(f, fullfile(plotpath, sprintf('Corr_ctrl%s.pdf', cond_leg{icond})), 'pdf')
saveas(f, fullfile(plotpath, sprintf('Corr_ctrl%s.png', cond_leg{icond})), 'png')

%% plot correlation
cfg=[];    cfg.layout = lay; cfg.parameter = 'rho'; cfg.colorbar = 'yes'; cfg.zlim = 'maxabs';
ft_multiplotTFR(cfg, corrstat_mse)
%% plot delta+f over time, check staircase convergence, use cleaned data
% figure; plot(cleaned.trialinfo(:,6))
f=figure;
for isub = 1:length(SUBJ)
  for icond = 1:2
    subplot(5,7,isub); hold on
    plot(trialinfo_trl{isub,icond}(:,6))
    title(SUBJ{isub})
  end
end

%% repeated measures correlation
% get mse for time scales and chans and times of interest
cfg=[];
cfg.frequency = [60 100];  cfg.avgoverfreq = 'yes';
cfg.channel = channel';    cfg.avgoverchan = 'yes';
cfg.latency = [0.1 0.3];   cfg.avgovertime = 'yes';
mse_sel = cellfun(@(x) ft_selectdata(cfg, x), mse_blocks);
% % only keep subjects > 3 runs
% mse_nruns = arrayfun(@(x) size(x.powspctrm,1), mse_sel(:,1));
% mse_sel = mse_sel(mse_nruns > 3,:);

% take difference for mse
cfg=[];
cfg.operation = 'subtract';
cfg.parameter = 'powspctrm';
mse_sel_diff = arrayfun(@(x,y) ft_math(cfg, x,y), mse_sel(:,1), mse_sel(:,2));
% demean
mse_sel_diff_demean = cellfun(@(x) x-mean(x), {mse_sel_diff.powspctrm}, 'UniformOutput',false);
% mse_sel_diff_demean = cellfun(@(x) x, {mse_sel_diff.powspctrm}, 'UniformOutput',false);
msedat = vertcat(mse_sel_diff_demean{:});

% take difference for accuracy
cfg=[];
cfg.operation = 'subtract';
cfg.parameter = 'trialinfo';
behav_sel_diff = arrayfun(@(x,y) ft_math(cfg, x,y), mse_sel(:,1), mse_sel(:,2));
%demean
behav_sel_diff_demean = cellfun(@(x) x-mean(x), {behav_sel_diff.trialinfo}, 'UniformOutput',false);
% behav_sel_diff_demean = cellfun(@(x) x, {behav_sel_diff.trialinfo}, 'UniformOutput',false);
behavdat = vertcat(behav_sel_diff_demean{:});
behavdat = behavdat(:,3);

[r,p] = corr(msedat, behavdat, 'Type','Pearson')

[r,p] = corr(msedat, behavdat, 'Type','Spearman')
figure; scatter(msedat, behavdat); lsline

