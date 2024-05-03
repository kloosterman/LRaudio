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

age_groups = {'YA' 'OA'};

% define subjects
SUBJ=cell(1,2);
for iage = 1:2
  for i=1:37
    SUBJ{iage}{end+1} = sprintf('%d', i);
  end
  SUBJbool = true(size(SUBJ{iage}));
  if strcmp(age_groups{iage}, 'YA')
    % SUBJbool([10, 12, 15, 17 ]) = false; % exclude 10 and 15, 17 missing, 12 conditions missing
    % exclude 10 and 15, 17 missing, 12 conditions missing,
    % #6 staircase not converging + datamissing, accuracy only 0.65
    % #31: staircase at max, accuracy only 0.6006
    % #3: staircase at minimum (0.01) in many trials
    % 34 sc at max until halfway, then drops, convergence weird?
    SUBJbool([37, 10, 12, 15, 17,    31, 6, 3 ]) = false;
  else
    % 7 and 34 don't exist
    SUBJbool([ 7, 34 ]) = false;
  end
  SUBJ{iage} = SUBJ{iage}(SUBJbool);
end
disp(SUBJ)

%% compute entropy, freqanalysis and ERPs
% make list of files to analyze on tardis TODO loop over age_groups
overwrite = 1;
cfglist = {}; cfg=[];
cfg.evoked = 'subtract'; % empty, regress, or subtract
cfg.csd = ''; % empty or csd
cfg.sensor_or_source = 'sensor';
cfg.runperblock = 'no'; % empty for all together, or per block.
for iage = 2
  for isub = 1:length(SUBJ{iage})
    cfg.SUBJ = SUBJ{iage}{isub};
    for icond = 1:2
      cfg.icond = icond;
      switch cfg.sensor_or_source
        case 'sensor'
          if strcmp(age_groups{iage}, 'YA')
            cfg.datafile = fullfile(datapath,  age_groups{iage}, 'data', SUBJ{iage}{isub}, sprintf('clean_SUB%s.mat', SUBJ{iage}{isub}));
          else
            cfg.datafile = fullfile(datapath,  age_groups{iage}, 'data', SUBJ{iage}{isub}, sprintf('NICOSA_%s_clean_task.mat', SUBJ{iage}{isub}));
          end
        case 'source'
          cfg.datafile = fullfile(datapath, 'source', SUBJ{iage}{isub}, sprintf('SourceTimeSeries_BW_1-100Hz_ParcelSpace_Block*.mat'));
      end
      if overwrite;      cfglist{end+1} = cfg;       end
    end
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
end
return

%% Merge mse files across subjects and cond
evoked = 'subtract'; % subtract_avgref subtract
csd = ''; % csd
runperblock = 'no';

eeg = [];
for iage = 1:2
  mse_tmp = {}; freq_tmp = {}; freq_tmp_bl = {}; timelock_tmp = {};
  trialinfo = []; trialinfo_trl = {}; trialinfo_blocks = {};
  for isub = 1:length(SUBJ{iage})
    for icond = 1:2
      path = fullfile(datapath, age_groups{iage}, 'mse', evoked, csd, sprintf('SUB%s_cond%d.mat', SUBJ{iage}{isub}, icond));
      disp(path);      load(path)
      mse.freq = mse.timescales;
      mse.powspctrm = mse.sampen;
      mse.dimord = 'chan_freq_time';
      mse_tmp{isub,icond} = mse;
      trialinfo(isub,:,icond) = mean(mse.trialinfo);
      trialinfo_trl{isub,icond} = mse.trialinfo;

      path = fullfile(datapath, age_groups{iage}, 'freq', evoked, csd, sprintf('SUB%s_cond%d.mat', SUBJ{iage}{isub}, icond));
      disp(path);      load(path)
      freq_tmp_bl{isub,icond} = freq_bl; % freq_bl (baseline corrected) freq (raw power)
      freq_tmp{isub,icond} = freq; % freq_bl (baseline corrected) freq (raw power)

      path = fullfile(datapath, age_groups{iage}, 'timelock', evoked, csd, sprintf('SUB%s_cond%d.mat', SUBJ{iage}{isub}, icond));
      disp(path);      load(path)
      timelock_tmp{isub,icond} = timelock;

      if strcmp(runperblock, 'yes')
        path = fullfile(datapath, age_groups{iage},  'mse', evoked, csd, sprintf('SUB%s_cond%d_blocks.mat', SUBJ{iage}{isub}, icond));
        disp(path);    load(path)
        for iblock=1:length(mse)
          mse{iblock}.freq = mse{iblock}.timescales;
          mse{iblock}.powspctrm = mse{iblock}.sampen;
          mse{iblock}.dimord = 'chan_freq_time';
          mse{iblock}.trialinfo = mean(mse{iblock}.trialinfo);
        end
        cfg=[];        cfg.appenddim = 'rpt';
        mse_blocks{isub,icond} = ft_appendfreq(cfg, mse{:}); 
      end
    end
  end
  % combine mse
  cfg=[];         cfg.keepindividual = 'yes';
  mse_merged{1} = ft_freqgrandaverage(cfg, mse_tmp{:,1});
  mse_merged{2} = ft_freqgrandaverage(cfg, mse_tmp{:,2});
  cfg.keepindividual = 'no';
  mse_merged{3} = ft_freqgrandaverage(cfg, mse_merged{1:2});

  cfg=[];       cfg.operation = 'subtract';     cfg.parameter = 'powspctrm';
  mse_merged{4} = ft_math(cfg, mse_merged{1}, mse_merged{2});
  % add behavior
  mse_merged{1}.trialinfo = trialinfo(:,:,1);
  mse_merged{2}.trialinfo = trialinfo(:,:,2);
  mse_merged{3}.trialinfo = mean(trialinfo, 3);
  mse_merged{4}.trialinfo = trialinfo(:,:,1) - trialinfo(:,:,2);

  % combine freq
  cfg=[];     cfg.keepindividual = 'yes';
  freq_merged{1} = ft_freqgrandaverage(cfg, freq_tmp{:,1});
  freq_merged{2} = ft_freqgrandaverage(cfg, freq_tmp{:,2});
  freq_merged_bl{1} = ft_freqgrandaverage(cfg, freq_tmp_bl{:,1});
  freq_merged_bl{2} = ft_freqgrandaverage(cfg, freq_tmp_bl{:,2});
  cfg.keepindividual = 'no';
  freq_merged{3} = ft_freqgrandaverage(cfg, freq_merged{1:2});
  freq_merged_bl{3} = ft_freqgrandaverage(cfg, freq_merged_bl{1:2});
  cfg=[];       cfg.operation = 'subtract';         cfg.parameter = 'powspctrm';
  freq_merged{4} = ft_math(cfg, freq_merged{1}, freq_merged{2});
  freq_merged_bl{4} = ft_math(cfg, freq_merged_bl{1}, freq_merged_bl{2});
  % add behavior
  freq_merged{1}.trialinfo = trialinfo(:,:,1);    freq_merged{2}.trialinfo = trialinfo(:,:,2);
  freq_merged{3}.trialinfo = mean(trialinfo, 3);  freq_merged{4}.trialinfo = trialinfo(:,:,1) - trialinfo(:,:,2);
  freq_merged_bl{1}.trialinfo = trialinfo(:,:,1);    freq_merged_bl{2}.trialinfo = trialinfo(:,:,2);
  freq_merged_bl{3}.trialinfo = mean(trialinfo, 3);  freq_merged_bl{4}.trialinfo = trialinfo(:,:,1) - trialinfo(:,:,2);

  % combine timelock
  cfg=[];        cfg.keepindividual = 'yes';
  timelock_merged{1} = ft_timelockgrandaverage(cfg, timelock_tmp{:,1});
  timelock_merged{2} = ft_timelockgrandaverage(cfg, timelock_tmp{:,2});
  % cfg.keepindividual = 'no';
  % cfg.parameter = 'individual';
  % timelock_merged{3} = ft_timelockgrandaverage(cfg, timelock_merged{1:2});
  cfg=[];      cfg.operation = 'subtract';       cfg.parameter = 'individual';
  timelock_merged{4} = ft_math(cfg, timelock_merged{1}, timelock_merged{2});
  % add behavior
  timelock_merged{1}.trialinfo = trialinfo(:,:,1);
  timelock_merged{2}.trialinfo = trialinfo(:,:,2);
  timelock_merged{3}.trialinfo = mean(trialinfo, 3);
  timelock_merged{4}.trialinfo = trialinfo(:,:,1) - trialinfo(:,:,2);

  eeg(iage).timelock = timelock_merged;
  eeg(iage).freq = freq_merged;
  eeg(iage).freq_bl = freq_merged_bl;
  eeg(iage).mse = mse_merged;
end

%% plot ERPs
cfg=[];   cfg.layout = lay;
% cfg.baseline = [-0.5 0];    cfg.baselinetype = 'relchange';
ft_multiplotER(cfg, timelock_merged{1}, timelock_merged{2}); colorbar
%% plot TFR / topo MSE
cfg=[];   cfg.layout = lay;
% cfg.zlim = 'maxabs';
% cfg.zlim = [1.17 1.23];
% cfg.xlim = [-0.5 1.5];
% cfg.baseline = [-0.5 0];    cfg.baselinetype = 'relchange';
cfg.colorbar = 'yes';
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
% zlim = [1.05 1.14];
% channel = {'FC5', 'FC3', 'FC1', 'FCz', 'FC2', 'FC4', 'C2', 'Cz', 'C1', 'C3', 'F3', 'F1'};
% channel = {'FC2', 'FCz', 'FC1', 'C3', 'C1', 'Cz', 'C2', 'CP1'}; % CSD
% channel = {'FC2', 'FCz', 'FC1', 'C1', 'Cz', 'C2'}; % avg ref 'C3',
channel = {'FC2', 'FCz', 'FC1', 'C3', 'C1', 'Cz', 'C2', 'CP2', 'CPz', 'CP1'}
channel = {'FCz', 'FC1', 'C1', 'Cz', 'CPz', 'CP1'}
% channel = {'FC2', 'FCz', 'FC1', 'C3', 'C1', 'Cz', 'C2', 'CP2', 'CPz', 'CP1'} % for csd
channel = {'FC2', 'FCz', 'FC1',  'C1', 'Cz', 'C2'} % for avgref  'C3', , 'CP2', 'CPz', 'CP1'
channel_freq = {'FCz',  'Cz', 'CPz', 'CP1'}; %   'CP3' 'C1',

f = figure; f.Position = [680         520        800         922];
% TFR mse 
cfg=[];   cfg.layout = lay;   cfg.figure = 'gcf'; cfg.channel = channel;  cfg.zlim = [1.15 1.23];  cfg.xlim = XLIM; cfg.title = ['Entropy YA'];
subplot(6,3,1);     ft_singleplotTFR(cfg, eeg(1).mse{3}); xline(0); ylabel('Time scale (ms)')
cfg.zlim = [1.05 1.14]; cfg.title = ['Entropy OA'];
subplot(6,3,4);     ft_singleplotTFR(cfg, eeg(2).mse{3}); xline(0); xlabel('Time from stim (s)'); ylabel('Time scale (ms)')
% topo
cfg=[];   cfg.layout = lay;   cfg.figure = 'gcf'; cfg.xlim = [0 0.3];   cfg.ylim = [70 90];   cfg.zlim = zlim; cfg.highlightchannel = channel; cfg.highlight = 'on';%cfg.zlim = [1.17 1.22];
subplot(6,3,2);     ft_topoplotTFR(cfg, eeg(1).mse{3}); %colorbar
cfg.zlim = [1.04 1.13];
subplot(6,3,5);     ft_topoplotTFR(cfg, eeg(2).mse{3}); %colorbar
% time series
cfg=[];    cfg.figure = 'gcf';  cfg.channel = channel; cfg.title = ' ';  cfg.xlim = XLIM;
subplot(6,3,3);     ft_singleplotER(cfg, eeg(1).mse{1:2}); hold on; xline(0,'HandleVisibility','off');
subplot(6,3,6);     ft_singleplotER(cfg, eeg(2).mse{1:2}); hold on; xline(0,'HandleVisibility','off');
legend({'Incong.', 'Cong.'}, 'Location', 'South'); legend boxoff
xlabel('Time from stim (s)'); ylabel('Entropy')
% ax=gca; YL = ax.YLim; patch([0; 0; 0.46; 0.46], [YL(1); YL(2); YL(2); YL(1);], [0.5 0.5 0.5 ])

% plot Central pooling freq
cfg=[];   cfg.layout = lay;   cfg.figure = 'gcf';   cfg.channel = channel;
cfg.xlim = XLIM; cfg.ylim = [0 100]; cfg.zlim = [-0.17 0.17]; cfg.title = 'Power modulation YA';
subplot(6,3,7); ft_singleplotTFR(cfg, eeg(1).freq_bl{3}); xline(0); xlabel('Time from stim (s)'); ylabel('Frequency (Hz)')
cfg.title = 'Power modulation OA';
subplot(6,3,10); ft_singleplotTFR(cfg, eeg(2).freq_bl{3}); xline(0); xlabel('Time from stim (s)'); ylabel('Frequency (Hz)')

cfg=[]; cfg.layout = lay;   cfg.figure = 'gcf'; cfg.highlightchannel = channel; cfg.highlight = 'on';
cfg.xlim = [0.2 0.4];   cfg.ylim = [4 8];    cfg.zlim = 'maxabs';
% cfg.xlim = [0.1 1.5];   cfg.ylim = [60 80];    cfg.zlim = 'maxabs';
subplot(6,3,8); ft_topoplotTFR(cfg, eeg(1).freq_bl{3}); %colorbar;
subplot(6,3,11); ft_topoplotTFR(cfg, eeg(2).freq_bl{3}); %colorbar;

cfg=[]; cfg.figure = 'gcf'; cfg.channel = channel; cfg.title = ' ';  cfg.frequency = [4 8]; cfg.xlim = XLIM;
% cfg=[]; cfg.figure = 'gcf'; cfg.channel = channel_freq; cfg.title = ' ';  cfg.frequency = [60 80];  cfg.xlim = XLIM;
subplot(6,3,9); ft_singleplotER(cfg, eeg(1).freq_bl{1:2});
subplot(6,3,12); ft_singleplotER(cfg, eeg(2).freq_bl{1:2});
xline(0,'HandleVisibility','off'); % legend({'Incong.', 'Cong.'}, 'Location', 'Southeast'); legend boxoff
xlabel('Time from stim (s)'); ylabel('Power (% signal change)'); % %shg

orient portrait
saveas(f, fullfile(plotpath, ['msevsfreq_' [channel{:}] ]), 'pdf') % '_' age_groups{iage}
saveas(f, fullfile(plotpath, ['msevsfreq_' [channel{:}] ]), 'png') % '_' age_groups{iage}

%% run stats: correlation mMSE vs behavior TODO correlate overall entropy vs delta_f!
behav_col=3;% 3 is accuracy,6 is delta_fsemitones
% corrtype='Pearson'; % Spearman Pearson
corrtype='Spearman'; % Spearman Pearson
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
    title(SUBJ{iage}{isub})
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

