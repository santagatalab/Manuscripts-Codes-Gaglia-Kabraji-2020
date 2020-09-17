function [initial_est, cutoff, multiplier, aic] = ...
        norm_unimodal(data_subset, peak_loc, prior, label, rep, FigSettings)

    % determine whether the data has a positive of negative tail
    if prior ~= 0
        tail = prior;
    else
        tail = 1;
        if sum(data_subset >= peak_loc) < sum(data_subset < peak_loc)  % more data points left of the peak
            if peak_loc-prctile(data_subset,2) > prctile(data_subset,98)-peak_loc 
                tail = -1;
            end
        end
    end

    if tail == 1
        sort_direction = 'descend';
        data_filtered = data_subset(data_subset > prctile(data_subset, 5) | data_subset > 1);
    else
        sort_direction = 'ascend';
        data_filtered = data_subset(data_subset < prctile(data_subset, 95));
    end

    mu_goal = peak_loc;

    % find cutoff point that maximizes the fit of a normal curve centered at the peak
    [~, sorted_idcs] = sort(data_filtered, sort_direction);
    initial_est = find_init_est(data_filtered, sorted_idcs, mu_goal, tail, label, rep, FigSettings);

    % this is a preliminary estimation of the positive population
    pos_perc = length(find(data_subset > initial_est))/length(data_subset);
    
    % set ceiling for positive percentage
    pos_perc = min([pos_perc, 0.98]);
        
    % GMM parameter setup
    S = {};
    if pos_perc < 0.02   % if positive population is very small, use k=3 for GMM
        pos_perc = 0.02;
        k=3;
        S.Sigma = zeros(1,1,k);
        S.Sigma(:,:,2) = cov(data_subset)/2;
        S.Sigma(:,:,[1 3]) = cov(data_subset)/10;
        lower_mu = prctile(data_subset, pos_perc);
        S.mu = [lower_mu mu_goal initial_est]';
        S.ComponentProportion = [pos_perc 1-(2*pos_perc) pos_perc];
    else
        neg_perc = 1-pos_perc;
        k=2;
        S.Sigma = zeros(1,1,k);
        S.Sigma(:,:,1) = cov(data_subset)/2;
        S.Sigma(:,:,2) = cov(data_subset)/2;
        if tail == 1
            S.mu = [mu_goal initial_est]';
            S.ComponentProportion = [1-pos_perc pos_perc];
        elseif tail == -1
            S.mu = [initial_est mu_goal]';
            S.ComponentProportion = [neg_perc 1-neg_perc];
        end
    end    

    [cutoff,multiplier,aic,pos_perc_updated] = fit_gmm(data_subset, tail, S, label, rep, FigSettings);

    % overexpression
    if k<3 && (pos_perc_updated > 0.7 || pos_perc_updated < 0.3)
        pos_perc_updated = max([pos_perc_updated, 0.05]);
        [cutoff,multiplier,aic] = fit_overexpr(data_subset,pos_perc_updated,label, rep, FigSettings);

        % unable to fit model --> default to original estimate
        if cutoff == 0 || isnan(cutoff)
            cutoff = initial_est;
            multiplier =  1;            % ????
            aic = 0;                    % todo: ??
        end
    end
end

