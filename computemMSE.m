function [mse] = computemMSE(cfg)

plotit = 0;

datafile = cfg.datafile;
SUBJ = cfg.SUBJ;
icond = cfg.icond;
outpath = cfg.outpath;
mse = [];

if exist(outpath)
  disp('Output file exists')
  load(outpath)
  return
end

disp(datafile)
load(datafile)
data = cleaned; clear cleaned
%   cfg=[];
%     cfg.viewmode = 'vertical';
%     ft_databrowser(cfg, cleaned)

disp 'select valid trials...'
cfg=[];
cfg.trials = ~isnan(data.trialinfo(:,7));
data = ft_selectdata(cfg, data);

if length(data.trial)<2
  disp('Too few trials found, abort')
  mse=[];
  return
end
% make stim-locked
cfg=[];
cfg.offset = - round(data.trialinfo(:,7)/4);
data = ft_redefinetrial(cfg, data);

disp crop
cfg=[];
cfg.latency = [-1 2.5];
data = ft_selectdata(cfg, data);

cfg=[];
cfg.trials = data.trialinfo(:,8) == icond-1;
data = ft_selectdata(cfg, data);

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


cfg = [];
cfg.m = 2;
cfg.r = 0.5;
cfg.timwin = 0.5;
cfg.toi = [-0.5 -0.4 -0.3 -0.2 -0.1 0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1 1.1 1.2 1.3 1.4 1.5 1.6 1.7 1.8 1.9 2];
cfg.toi = [-0.5];
% cfg.timescales = 1:40;
    cfg.timescales = 20;
cfg.recompute_r = 'perscale_toi_sp';
cfg.coarsegrainmethod = 'filtskip';
cfg.filtmethod = 'lp';
cfg.mem_available = 20e+09;
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
  mse_noerp.freq = mse_noerp.timescales;
  mse_noerp.powspctrm = mse_noerp.sampen;
  mse_noerp.dimord =  'chan_freq_time';
  mse_noerpregress.freq = mse_noerpregress.timescales;
  mse_noerpregress.powspctrm = mse_noerpregress.sampen;
  mse_noerpregress.dimord =  'chan_freq_time';
  mse_diff = mse_noerp;
  mse_diff.powspctrm = mse.powspctrm - mse_noerpregress.powspctrm;
  cfg=[];
  cfg.layout = 'acticap-64ch-standard2.mat';
  %     cfg.zlim = [0.8 1.2];
  cfg.colorbar = 'yes';
  cfg.baseline = [-0.5 0];
  %     cfg.baselinetype = 'relchange';
  cfg.zlim = 'maxabs';
  %     cfg.xlim = [-0.5 1.5];
  ft_multiplotTFR(cfg,mse_diff)
end

disp(outpath)
save(outpath, 'mse')