function [mse] = computemMSE(cfg)

plotit = 0;

datafile = cfg.datafile;
SUBJ = cfg.SUBJ;
icond = cfg.icond;
outpath = cfg.outpath;
analysis = cfg.analysis;
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
%   cfg=[];
%     cfg.viewmode = 'vertical';
%     ft_databrowser(cfg, data)

switch csd
  case 'csd'
    cfg=[];
    cfg.method = 'spline';
    cfg.elec = 'standard_1020.elc';
    data = ft_scalpcurrentdensity(cfg, data);
end

if plotit
  load('acticap-64ch-standard2.mat')
  lay.label(find(contains(lay.label, 'Ref'))) = {'FCz'};
  cfg = [];
  cfg.layout = lay;
  ft_multiplotER(cfg, data, scd)
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

switch analysis
  
  case 'mse'
    mse = [];
    cfg = [];
    cfg.m = 2;
    cfg.r = 0.5;
    cfg.timwin = 0.5;
    cfg.toi = [-0.5 -0.4 -0.3 -0.2 -0.1 0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1 1.1 1.2 1.3 1.4 1.5 1.6 1.7 1.8 1.9 2];
    % cfg.toi = [-0.5];
    %     cfg.timescales = 1:40;
%     cfg.timescales = 10:30;
    cfg.timescales = 40;
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
    
    disp(outpath)
    save(outpath, 'mse')

  case 'freq'
    
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
    
    disp(outpath)
    save(outpath, 'freq', 'freq_bl')

  case 'erp'
    disp 'average trials'
    cfg=[];
    timelock= ft_timelockanalysis(cfg, data);
    disp(outpath)
    save(outpath, 'timelock')

end
