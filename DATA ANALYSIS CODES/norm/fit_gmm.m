function [add_fact,mult_fact,aic,pos_perc] = fit_gmm(data_subset,tail,S,label, rep, FigSettings)

    try
        x = linspace(min(data_subset),max(data_subset),1000);
        k = length(S.mu);
        options = statset('MaxIter',2000,'TolFun',1e-8);  %'Display','final',
        model1 = fitgmdist(data_subset,k,'Start',S,'Options',options,'RegularizationValue', 0.0001);
        aic = model1.AIC;
        pos_perc = model1.ComponentProportion(k);
        
        t = linspace(min(data_subset),max(data_subset),100);
        [n,h] = hist(data_subset,t);
        if FigSettings.FigFlag
            f = figure; %('visible',fig_display);
            plot(h,n/sum(n)/(t(2)-t(1)))
        end

        all_probs = [];
        all_comps = [];
        all_means = [];
        
        for clst=1:k
            clst_mean = model1.mu(clst);
            clst_std = sqrt(model1.Sigma(:,:,clst));  % TODO: does this make sense?
            clst_comp = model1.ComponentProportion(clst);
            clst_prob = normpdf(x,clst_mean,clst_std)*clst_comp;
            
            all_probs = [all_probs; clst_prob];
            all_comps = [all_comps; clst_comp];
            all_means = [all_means; clst_mean];
            if FigSettings.FigFlag
                hold on
                plot(x,clst_prob)
                ax = axis;
                hold on
                plot([S.mu(clst) S.mu(clst)], [ax(3) ax(4)], '--r');
            end
        end
        [~, order] = sort(all_means);
        all_probs = all_probs(order,:);
        all_comps = all_comps(order);
        all_means = all_means(order);


        if tail == 1
            posprob = all_probs(size(all_probs,1),:);
            negprob = all_probs(size(all_probs,1)-1,:);
            pos_comp = all_comps(size(all_probs,1),:);
            neg_comp = all_comps(size(all_probs,1)-1,:);
            pos_mean = all_means(size(all_probs,1),:);
            neg_mean = all_means(size(all_probs,1)-1,:);
        else %if tail == -1
            posprob = all_probs(2,:);
            negprob = all_probs(1,:);
            pos_comp = all_comps(2,:);
            neg_comp = all_comps(1,:);
            pos_mean = all_means(2,:);
            neg_mean = all_means(1,:);
        end


        % check that the intersection of the function exists
        close_to_zero = 0.01;
        ind_1 = find(abs(posprob-negprob)./((posprob+negprob)/2) < close_to_zero);
        % make sure there are point and if not increase the tolerance
        count2 = 0;
        tolerance = 0.02;  
        maxcount = 10;
        while length(ind_1)<1 && count2 < maxcount
            count2 = count2 + 1;
            ind_1 = find(abs(posprob-negprob)./((posprob+negprob)/2)<(0.1+tolerance*count2));
        end

        if count2 == maxcount
            disp('The fitting is not working')
        end

        [~,ind_check] = max([neg_comp pos_comp]);

        temp = x(ind_1);
        if ind_check == 1
            ind_1 = ind_1(temp > neg_mean);
        elseif ind_check == 2
            ind_1 = ind_1(temp < pos_mean);
        end


        % now find the exact minimum
        [~,ind_2] = max(posprob(ind_1).*negprob(ind_1));
        % and find the location of the intersection
        inter_point = ind_1(ind_2);
        x_inter = x(inter_point);

        % check the the intersection of the two curves is not in a meaningless
        % point of the distribution
        [~,ind_check] = max([neg_comp pos_comp]);
        if ~isempty(x_inter)
            if (ind_check == 1 && x_inter < neg_mean) || (ind_check == 2 && x_inter > pos_mean)
                % if true the intersection happens in the wrong place, try again
                x_inter = [];
            end
        end
            
        if ~isempty(x_inter)
            % set the output parameters
            add_fact = x_inter;
            mult_fact = abs(pos_mean-neg_mean);
            
            if FigSettings.FigFlag
                hold on
                ax = axis;
                plot([x_inter x_inter], [ax(3) ax(4)], '--g');
                title([label ' - GMM']) 
                filename = [FigSettings.Folder label '\' num2str(rep) '_gmm.fig'];
                saveas(f, filename, 'fig')
                if ~FigSettings.Debug
                    close(f)
                end
            end
        else
            add_fact = 0;
            mult_fact = 0;
            if FigSettings.FigFlag
                title([label ' - GMM']) 
                filename = [FigSettings.Folder label '\' num2str(rep) '_gmm.fig'];
                saveas(f, filename, 'fig')
                if ~FigSettings.Debug
                    close(f)
                end
            end
        end

    catch exception
        disp('There was an error fitting the Gaussian mixture model')
        disp(exception.message)
        
        add_fact = 0;
        mult_fact = 0;
        aic = 0;
    end
end