function [add_fact,mult_fact,aic] = fit_overexpr(data, estimate, label, rep, FigSettings)

    count1 = 0;
    count2 = 0;
    x_inter = [];
    options = statset('MaxIter',2000,'TolFun',1e-8);

    x = linspace(min(data),max(data),1000);

    while length(x_inter) == 0 && count1 < 10
        count1 = count1 + 1;
        % if the population is likely to be monotype plus rare populations then
        % we fit a single gaussian first   

        if estimate <= 0.5
            filt_data = data(data < prctile(data,(1-estimate)*100));
            model1 = fitgmdist(filt_data,1,'Options',options);
            lim_mean = model1.mu + sqrt(model1.Sigma(:,:,1))*2;
            s = -1;
        else
            filt_data = data(data > prctile(data,(1-estimate)*100));
            model1 = fitgmdist(filt_data,1,'Options',options);
            lim_mean = model1.mu - sqrt(model1.Sigma(:,:,1))*2;
            s = +1;
        end
        
        lim_std = sqrt(model1.Sigma(:,:,1))*2;
        aic = model1.AIC;

        obs = data; 

        fun = @(a) -sum( log(  a(1)*normpdf(obs,model1.mu,model1.Sigma(:,:,1)^(1/2)) + (1-a(1))*normpdf(obs,a(2),a(3))   )   );   

        % constaints on the 3rd component - a gaussian model for overexpression
        A = [-1  0  0  ;
              0 -1  0  ;
              0  0 -1  ;
              1  0  0  ;
              0  s  0  ;
              0  0  1  ];  
        b = [0 0 0 1 s*lim_mean lim_std ];
        options_fmincon = optimoptions('fmincon','Display','off');
        par = fmincon(fun,[0.9 lim_mean-0.05 lim_std-0.05],A,b,[],[],[],[],[],options_fmincon);
%         [par,~,~,~] = fmincon(fun,[0.9 lim_mean-0.05 lim_std-0.05],A,b);

        if estimate < 0.5
            neg_mean = model1.mu(1);
            pos_mean = par(2);
            neg_std = model1.Sigma(:,:,1)^(1/2);
            pos_std = par(3);
            neg_comp = par(1);
            pos_comp = (1-par(1));
        else
            neg_mean = par(2);
            pos_mean = model1.mu(1);
            neg_std = par(3);
            pos_std = model1.Sigma(:,:,1)^(1/2);
            neg_comp = (1-par(1));
            pos_comp = par(1);
        end
        negprob = neg_comp*normpdf(x,neg_mean,neg_std);
        posprob = pos_comp*normpdf(x,pos_mean,pos_std);


        % check that the intersection of the function exists
        ind_1 = find(abs(posprob-negprob)./((posprob+negprob)/2)<0.1);
            % make sure there are point and if not increase the tolerance
        while length(ind_1)<1 && count2 < 10
            count2 = count2 + 1;
            ind_1 = find(abs(posprob-negprob)./((posprob+negprob)/2)<(0.1+0.02*count2));
        end

        if count2 == 10
            disp('The fitting is not working')
%             normdata = zeros(size(data_in))+NaN;
            add_fact = NaN;
            mult_fact = NaN;
%             data_prob = zeros(size(data_in))+NaN;
            return
        end

        [~,ind_check] = max([neg_comp pos_comp]);

        temp = x(ind_1);
        if ind_check == 1
            ind_1 = ind_1(temp > neg_mean);
        elseif ind_check == 2
            ind_1 = ind_1(temp < pos_mean);
        end


        % now find the exact minimum
    %     [~,ind_2] = min(abs(posprob(ind_1)-negprob(ind_1)));
        [~,ind_2] = max(posprob(ind_1).*negprob(ind_1));
        % and find the location of the intersection
        inter_point = ind_1(ind_2);
        x_inter = x(inter_point);

        % check the the intersection of the two curves is not in a meaningless
        % point of the distribution
        [~,ind_check] = max([neg_comp pos_comp]);
        if length(x_inter) > 0
            if (ind_check == 1 && x_inter < neg_mean) || (ind_check == 2 && x_inter > pos_mean)
                % if true the intersection happens in the wrong place, try again
                x_inter = [];
            end
        end
    end

    if count1 == 10
        disp('The fitting is not working')
        add_fact = NaN;
        mult_fact = NaN;
        return
    end


        %%% calculate the "positive" posterior probability
    if length(x_inter) > 0
        % set the output parameters
        add_fact = x_inter;
        mult_fact = abs(pos_mean-neg_mean);

    else
        % if the gmm could not find a solution where the two gaussians 
        % intersected than we use the tail of the max

        disp('The fitting is not working')
        add_fact = NaN;
        mult_fact = NaN;
        return
    end
    
    if FigSettings.FigFlag == 1
        f = figure; %('visible',fig_display);
        
        t = linspace(min(data),max(data),60);
        [n,h] = hist(data,t);
        plot(h,n/sum(n)/(t(2)-t(1)))
        hold on
        plot(x,posprob)
        plot(x,negprob)
        hold on
        ax = axis;
        plot([add_fact add_fact], [ax(3) ax(4)], '--r')
        
        title([label ' - overexpression fit']) 
        filename = [FigSettings.Folder label '\' num2str(rep) '_overexpr.fig'];
        saveas(f, filename, 'fig')
        if ~FigSettings.Debug
            close(f)
        end
    end
end
