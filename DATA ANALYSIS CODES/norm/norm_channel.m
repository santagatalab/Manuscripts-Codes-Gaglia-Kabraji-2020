% data - full dataset to normalize on
% cell_num - number of cells to subset on (so algorithm doesn't take too long)
% reps - number of reps to run of the normalization algorithm
% prior - default 0, if 1 --> positive tail, if -1 --> negative tail
% overexpr - positive population percentage

function [best_cutoff, best_mult] = norm_channel(data, cell_num, reps, prior, overexpr, final_est, label, FigSettings)

    good_idcs = find(data > prctile(data,0.1) & data < prctile(data,99.9));  % remove extreme outliers
    cutoffs = zeros(reps,1);
    initial_estimates = zeros(reps,1);
    multipliers = zeros(reps,1);
    aics = zeros(reps,1);
    npoints0 = 100;
    
    if FigSettings.FigFlag
        mkdir([FigSettings.Folder label])
    end
    
    for r=1:reps
        subset_idcs = good_idcs(round((length(good_idcs)-1).*rand(cell_num,1) + 1)); % define subset of cells to normalize on
        data_subset = data(subset_idcs);
        
        if final_est > 0
            cutoff = prctile(data_subset, (1-final_est)*100);
            mult = 1;
            aic = 0;
            initial_est = 0;
        
        elseif overexpr > 0
            [cutoff, mult, aic] = fit_overexpr(data_subset, overexpr, label, r, FigSettings);
            initial_est = 0;
        else
            [f0,x0] = ksdensity(data_subset, 'npoints', npoints0, 'Bandwidth', 0.2);
            [pks,locs,~,~] = findpeaks(f0,'MinPeakDistance',5,'MinPeakHeight',0.05,'MinPeakWidth',3);

            if length(pks) == 1  % unimodal
                peak_idx = locs(find(pks == max(pks)));
                peak_loc = x0(peak_idx);
                [initial_est, cutoff, mult, aic] = norm_unimodal(data_subset, peak_loc, prior, label, r, FigSettings);

            elseif length(pks) == 2  % bimodal distribution
                peak_locs = sort(x0(locs));
                [cutoff, mult, aic] = norm_bimodal(data_subset, peak_locs, label, r, FigSettings);
                initial_est = 0;

            else
                peak_xlow = x0(locs(find(locs == min(locs))));
                peak_xhigh = x0(locs(find(locs == max(locs))));
                peak_locs = [peak_xlow, peak_xhigh];
                disp([num2str(peak_xlow) ', ' num2str(peak_xhigh)])
                [cutoff, mult, aic] = norm_bimodal(data_subset, peak_locs, label, r, FigSettings);
                initial_est = 0;
            end
        end
        
        cutoffs(r) = cutoff;
        initial_estimates(r) = initial_est;
        multipliers(r) = mult;
        aics(r) = aic;
        
        disp([num2str(r) ' - ' label ': ' num2str(cutoff) ', ' num2str(initial_est) ', ' ...
            num2str(mult) ', ' num2str(aic)])
    end 
    
    % filter out outliers and undefined cutoffs
    bad_reps = find(cutoffs == 0 | isnan(cutoffs));  
    cutoffs(bad_reps) = [];
    initial_estimates(bad_reps) = [];
    multipliers(bad_reps) = [];
    aics(bad_reps) = [];
    
    if isempty(cutoffs)
        best_cutoff = 0;
        best_mult = 1;
    elseif (std(cutoffs) > 1 && std(cutoffs) > std(initial_estimates))
        best_cutoff = median(initial_estimates);
        [~,best_idx] = min(abs(initial_estimates - best_cutoff));
        best_mult = multipliers(best_idx);
    else
        best_idx = find(aics == max(aics));
        best_cutoff = cutoffs(best_idx(1));
        [~,best_idx] = min(abs(cutoffs - best_cutoff));
        best_mult = multipliers(best_idx);
    end
    
    best_mult = max([best_mult, 0.1]); % multiplier can't be less than 0.1
    
    disp(['BEST: ' num2str(best_cutoff) ', ' num2str(best_mult)])
    disp('--------------------------------------------------')
end

