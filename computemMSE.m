function [mse] = computemMSE(cfg)
global datapath

plotit = 0;

datafile = cfg.datafile;
SUBJ = cfg.SUBJ;
icond = cfg.icond;
evoked = cfg.evoked;
csd = cfg.csd;
sensor_or_source = cfg.sensor_or_source;

switch sensor_or_source
  case 'sensor'
    disp(datafile)
    load(datafile)
    data = cleaned; clear cleaned
  case 'source'
    data = {};
    list = dir(datafile);
    for i = 1:length(list)
      disp(fullfile(list(i).folder, list(i).name))
      load(fullfile(list(i).folder, list(i).name))
      ntrials = size(sourcedata.trialinfo,1);
      sourcedata.trial = mat2cell(sourcedata.trial, ones(ntrials,1));
      sourcedata.trial = cellfun(@squeeze, sourcedata.trial, 'Uni', 0)';
      removefields(sourcedata, {'dimord'}); %  'sampleinfo'
      data{i} = sourcedata;   clear sourcedata
    end
    cfg = [];
    data = ft_appenddata(cfg, data{:});
end

if plotit && ismac
  ft_databrowser([],data)
  load('Acticap_64_UzL.mat')
%   lay.label(find(contains(lay.label, 'Ref'))) = {'FCz'};
  cfg = [];
  cfg.layout = lay;
  ft_multiplotER(cfg, data)
end

disp 'select valid trials...'
cfg=[];
cfg.trials = ~isnan(data.trialinfo(:,7));
data = ft_selectdata(cfg, data);

if length(data.trial)<2
  disp('Too few trials found, abort')
  mse=[];
  return
end

disp 'make stim-locked'
% data_ori=data;
cfg=[];
% cfg.offset = -round(data.trialinfo(:,7)/4);
compression_factor = 44100/48000;
cfg.offset = -round(((data.trialinfo(:,7)/1000)+0.5)*data.fsample*compression_factor);
%figure; histogram(cfg.offset)
data = ft_redefinetrial(cfg, data);

disp 'crop trials'
cfg=[];
cfg.latency = [-1 2.5];
data = ft_selectdata(cfg, data);

disp 'reject trials with large variance'
par = [];   par.badtrs = []; par.method = 'zscorecut';
to_plot=0; if ismac; to_plot=1; end
keeptrls = EM_ft_varcut3(data, par, to_plot);
cfg=[];
cfg.trials = keeptrls;
data=ft_selectdata(cfg, data);

if ismac && plotit
  ft_rejectvisual([],data)
  cfg=[];  cfg.viewmode = 'vertical'; cfg.channel = 'C*'; ft_databrowser(cfg, data)  % 
end

% introduce implicit reference DONE PER TRIAL? run after trial rejection to
% be sure
cfg=[];
cfg.channel=1:63;
cfg.reref='no';
cfg.implicitref='TP9';
data = ft_preprocessing(cfg, data);
% make average reference
cfg=[];
cfg.reref='yes';
cfg.refchannel='all';
% cfg.refchannel={'TP9', 'TP10'};
data = ft_preprocessing(cfg, data);

% if ismac
%   ft_databrowser([], ft_timelockanalysis([],data))
%   timelock=ft_timelockanalysis([],data);
%   timelock_ori=ft_timelockanalysis([],data_ori);
%   figure; hold on; plot(timelock.time, mean(timelock.avg)); xlim([-0.5 1]); title('stim-locked'); xline(0)
%   plot(timelock_ori.time, mean(timelock_ori.avg)); xlim([-0.5 1]); title('stim-locked')
%   legend({'stim-locked' 'cue-locked'}); ax=gca; ax.XTick=[-0.5:0.1:1.5 ];
% end

disp 'select trials for cond'
cfg=[];
cfg.trials = data.trialinfo(:,8) == icond-1;
data = ft_selectdata(cfg, data);

% if plotit
%   ft_rejectvisual([],data)
%   cfg=[];  cfg.viewmode = 'vertical';  ft_databrowser(cfg, data)
% end

switch csd
  case 'csd'
    cfg=[];
    cfg.method = 'spline';
    cfg.elec = 'standard_1020.elc';
    data = ft_scalpcurrentdensity(cfg, data);
end

% erp analysis
timelock= ft_timelockanalysis([], data);
outpath = fullfile(fileparts(datapath), 'timelock', evoked, csd, sprintf('SUB%s_cond%d.mat', SUBJ, icond));
mkdir(fileparts(outpath))
disp(outpath)
save(outpath, 'timelock')

switch evoked
  case 'regress'
    timelock= ft_timelockanalysis([], data);
    timelock.avg(:, timelock.time < 0) = 0;
    timelock.avg(:, timelock.time > 1.5) = 0;
    disp 'regress ERP from single trials'
    for itrial = 1:size(data.trial,2)
      for ichan = 1:63 % TODO weigh ERP timepoints
        [~,~,res] = regress(data.trial{itrial}(ichan,:)', [timelock.avg(ichan,:)' ones(size(timelock.avg(ichan,:)'))] );
        %         figure; hold on;
        %         plot(timelock.time, data.trial{itrial}(ichan,:)');
        %         plot(timelock.time, timelock.avg(ichan,:)');
        %         plot(timelock.time, res);
        %         legend({'raw' 'erp' 'regressed'})
        data.trial{itrial}(ichan,:) = res;
      end
    end
  case 'subtract'
    timelock= ft_timelockanalysis([], data);
    for itrial = 1:size(data.trial,2)
      data.trial{itrial} = data.trial{itrial} - timelock.avg;
    end
  otherwise
    disp 'ERP not removed'
end

% cfg=[];
% cfg.bpfilter      = 'yes';
% cfg.bpfreq        = [4 8];
% data_theta = ft_preprocessing(cfg, data);
% data_thetaremoved = data;
% data_thetaremoved.trial = cellfun(@(x,y) x-y, data.trial, data_theta.trial, 'uni',0);
% if ismac & plotit
%   cfg=[];
%   cfg.viewmode = 'vertical';
%   ft_databrowser(cfg, data_theta)
%   ft_databrowser(cfg, data)
%   ft_databrowser(cfg, data_thetaremoved)
% end


if ismac && plotit
  disp 'average trials'
  cfg=[];
  %       cfg.covariance = 'yes';
  timelock= ft_timelockanalysis(cfg, data);
  cfg=[];
  cfg.parameter = 'var';
  cfg.layout = 'EEG1005.lay';
  ft_multiplotER(cfg,timelock)
  
  disp 'subtract ERP from single trials'
  data_noERP = data;
  data_noERPregress = data;
  for itrial = 1:size(data_noERP.trial,2)
    data_noERP.trial{itrial} = data_noERP.trial{itrial} - timelock.avg;
    for ichan = 1:63
      [~,~,res] = regress(data_noERP.trial{itrial}(ichan,:)', timelock.avg(ichan,:)');
      data_noERPregress.trial{itrial}(ichan,:) = res;
    end
  end
  timelock_noERP = ft_timelockanalysis([], data_noERP);
  timelock_noERPregress = ft_timelockanalysis([], data_noERPregress);
  cfg=[];
  %     cfg.parameter = 'var';
  cfg.layout = 'EEG1005.lay';
  ft_multiplotER(cfg,  timelock_noERPregress)
end

%freqanalysis
cfg              = [];
cfg.output       = 'pow';
cfg.channel      = 'EEG';
cfg.method       = 'mtmconvol';
cfg.taper        = 'hanning';
cfg.foi          = logspace(0, 1.5, 30); %2:1:30;                         % analysis 2 to 30 Hz in steps of 2 Hz
cfg.t_ftimwin    = ones(length(cfg.foi),1).*0.5;   % length of time window = 0.5 sec
cfg.toi          = -0.5:0.05:2;                  % time window "slides" from -0.5 to 1.5 sec in steps of 0.05 sec (50 ms)
freq{1} = ft_freqanalysis(cfg, data);

cfg              = [];
cfg.output       = 'pow';
cfg.channel      = 'EEG';
cfg.method       = 'mtmconvol';
cfg.taper        = 'dpss';
cfg.foi          = logspace( 1.5, 2, 10); % 30:4:100;
cfg.tapsmofrq    = ones(length(cfg.foi),1).*4;   % length of time window = 0.5 sec
cfg.t_ftimwin    = ones(length(cfg.foi),1).*0.5;   % length of time window = 0.5 sec
cfg.toi          = -0.5:0.05:2;                  % time window "slides" from -0.5 to 1.5 sec in steps of 0.05 sec (50 ms)
freq{2} = ft_freqanalysis(cfg, data);

cfg=[];
cfg.parameter = 'powspctrm';
cfg.appenddim = 'freq';
freq = ft_appendfreq(cfg, freq{:});

cfg=[];
cfg.baseline = [-0.5 0];
cfg.baselinetype = 'relchange';
freq_bl = ft_freqbaseline(cfg, freq);

if ismac & plotit
  load('acticap-64ch-standard2.mat')
  lay.label(find(contains(lay.label, 'Ref'))) = {'FCz'};
  cfg=[];
  cfg.layout = lay;
  cfg.zlim = 'maxabs';
  ft_multiplotTFR(cfg,freq_bl); colorbar
end

outpath = fullfile(fileparts(datapath), 'freq', evoked, csd, sprintf('SUB%s_cond%d.mat', SUBJ, icond));
mkdir(fileparts(outpath))
disp(outpath)
save(outpath, 'freq', 'freq_bl')

%mse analysis
cfg = [];
cfg.m = 2;
cfg.r = 0.5;
cfg.timwin = 0.5;
cfg.toi = [-0.5 -0.4 -0.3 -0.2 -0.1 0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1 1.1 1.2 1.3 1.4 1.5 1.6 1.7 1.8 1.9 2];
% cfg.toi = [-0.5];
%     cfg.timescales = 1:40;
cfg.timescales = 10:30;
%     cfg.timescales = 40;
cfg.recompute_r = 'perscale_toi_sp';
cfg.coarsegrainmethod = 'filtskip';
cfg.filtmethod = 'lp';
cfg.mem_available = 40e+09;
cfg.allowgpu = true;

% cfg.trials = 1:50;
%     mse{isub,icond} = ft_entropyanalysis(cfg, data);
mse = ft_entropyanalysis(cfg, data);
%     mse_noerp = ft_entropyanalysis(cfg, data_noERP);
% mse_noerpregress = ft_entropyanalysis(cfg, data_noERPregress);
%
if ismac & plotit
  mse.freq = mse.timescales;
  mse.powspctrm = mse.sampen;
  mse.dimord =  'chan_freq_time';
  %   mse_noerp.freq = mse_noerp.timescales;
  %   mse_noerp.powspctrm = mse_noerp.sampen;
  %   mse_noerp.dimord =  'chan_freq_time';
  %   mse_noerpregress.freq = mse_noerpregress.timescales;
  %   mse_noerpregress.powspctrm = mse_noerpregress.sampen;
  %   mse_noerpregress.dimord =  'chan_freq_time';
  %   mse_diff = mse_noerp;
  %   mse_diff.powspctrm = mse.powspctrm - mse_noerpregress.powspctrm;
  cfg=[];
  cfg.layout = 'acticap-64ch-standard2.mat';
  %     cfg.zlim = [0.8 1.2];
  cfg.colorbar = 'yes';
  %   cfg.baseline = [-0.5 0];
  %       cfg.baselinetype = 'relchange';
  %   cfg.zlim = 'maxabs';
  %     cfg.xlim = [-0.5 1.5];
  ft_multiplotTFR(cfg,mse)
end

outpath = fullfile(fileparts(datapath), 'mse', evoked, csd, sprintf('SUB%s_cond%d.mat', SUBJ, icond));
mkdir(fileparts(outpath))
disp(outpath)
save(outpath, 'mse')


