function [cutoff] = find_init_est(data_orig, sorted_idcs, mu_goal, tail, label, rep, FigSettings)

    stop = 0;
    filtered_idcs = [];
    idx = 1;
    interval = 20;
    data_sorted = data_orig(sorted_idcs);
    
    pd = fitdist(data_sorted,'Normal'); 

    while stop == 0 && ~isempty(data_sorted) && idx < length(data_orig)  
        if (tail == 1 && pd.mu < mu_goal) || (tail == -1 && pd.mu > mu_goal)
            stop = 1;
        else
            stop_idx = min([interval, length(data_sorted)]);
            data_sorted(1:stop_idx) = [];   
            pd = fitdist(data_sorted,'Normal');
            filtered_idcs = [filtered_idcs; sorted_idcs(idx:idx+stop_idx-1)];
            idx = idx + interval;
        end 
    end
    
    if tail == -1
        cutoff = max(data_orig(filtered_idcs));
        if isempty(cutoff)
            cutoff = min(data_orig);
        end
    else
        cutoff = min(data_orig(filtered_idcs));
        if isempty(cutoff)
            cutoff = max(data_orig);
        end
    end
    
    if FigSettings.FigFlag
        f = figure; %('visible',fig_display);
        [y,x] = ksdensity(data_orig);
        ax = axis;
        plot(x, y);
        hold on
        plot([mu_goal mu_goal], [ax(3) ax(4)], '--b')
        hold on
        plot([cutoff cutoff], [ax(3) ax(4)], '--r')
        title([label ' - initial cutoff estimate']) 
        filename = [FigSettings.Folder label '\' num2str(rep) '_init_est.fig'];
        saveas(f, filename, 'fig')
        if ~FigSettings.Debug
            close(f)
        end
    end
end

