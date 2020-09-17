
function [cutoff, multiplier, aic] = norm_bimodal(data_subset, peak_locs, label, rep, FigSettings)

    midpt = mean(peak_locs);

    subset_left = data_subset(find(data_subset < midpt));
    subset_right = data_subset(find(data_subset >= midpt));
    k=2;
    S = {};
    S.Sigma = zeros(1,1,k);
    S.mu = [peak_locs(1) peak_locs(2)]';
    S.ComponentProportion = [length(subset_left) length(subset_right)]/length(data_subset);
    
    if min(S.ComponentProportion) < 0.3  % use the overexpression model
        [cutoff,multiplier,aic] = fit_overexpr(data_subset, S.ComponentProportion(2), label, rep, FigSettings);
    else    
        S.Sigma(:,:,1) = cov(subset_left);
        S.Sigma(:,:,2) = cov(subset_right);
        tail = 0;
        [cutoff,multiplier,aic] = fit_gmm(data_subset,tail,S,label, rep, FigSettings);
    end
end

