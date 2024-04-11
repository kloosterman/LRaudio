%% Usage   : LCMV beamformer pipeline for neuromag dataset using FieldTrip
%% Scripted: Amit Jaiswal @ MEGIN, Helsinki, Finland 
%% Created : 15 Dec. 2018, Version: 1.0.3
%% File No : 
%% ***********************************************************************
function [selecttrials, par] = EM_ft_varcut3(data_summary, par, to_plot)
    fprintf('\nDetecting the trials with excessive max/min variance...\n')
    clear trl_var
    trlindx = 1:size(data_summary.trial,2);
     for trl = trlindx
         trl_var(trl,:)    = max(var(data_summary.trial{1,trl}(:,:)'));
     end
     rest_idx = setdiff(trlindx, par.badtrs);
     if strcmp(par.method, 'maxmin_perct')
       percentiles = prctile(trl_var(rest_idx), par.cov_cut);
       outlr_idx = trl_var(rest_idx) < percentiles(1) | trl_var(rest_idx) > percentiles(2);
     elseif  strcmp(par.method, 'zscorecut') % NKEDIT, remove zscore > 3
       zmaxvar = zscore(log(trl_var(rest_idx)));
       % % TOEGEVOEGD DOOR JOHANNES OM OP DE abs(min(Z)+cutoff) TE THRESHOLDEN
       threshold = median(zmaxvar) + abs(min(zmaxvar)-median(zmaxvar)) + abs(1); % add 1 std to make lenient
       outlr_idx = abs(zscore(log(trl_var(rest_idx)))) > threshold;
     end
     bd_trl_var = trlindx(outlr_idx);
     bd_trls = bd_trl_var;
     % disabled zscore:     bd_trls = union(bd_trl_var, bd_trl_zscore);
     par.bad_trials = sort([par.badtrs, bd_trls]);
     if to_plot
         figure
         subplot(111),
         scatter(trlindx, trl_var, 400, 'bD'); title('Max. variance'), hold on
         scatter(par.bad_trials, trl_var(par.bad_trials), 400, 'ro', 'linewidth',2); 
         xlim([trlindx(1) trlindx(end)]), title('Max. variance')
         legend({'all trials','bad trials'}, 'fontsize', 15)
     end
     fprintf('Found total %d trials to remove ...\n',length(par.bad_trials))
     selecttrials = setdiff(trlindx, par.bad_trials);
end
