function [Y, data_filtered, filtered_idcs, explained] = dim_reduction(data, dim, filterOn, PCAon)

    %%% Define distance matrix 
    spearman_corrs = abs(corr(data', data', 'Type','Pearson')); %, 'Type','Spearman'
    dist_matrix = 1-(spearman_corrs/2);
    disp('DONE - distance matrix')


    %%% Get 10 nearest neighbors for each cell, then define "outlierness" of each cell (average of spearman rank correlation of 10 nearest neighbors)
    [min_dists, min_idcs] = mink(dist_matrix,11,2);  % 11 closest neighbors because it include self

    outlierness = zeros(length(data), 1);
    for i=1:length(data)
        nn10_idcs = min_idcs(i,2:11);
        outlierness(i) = mean(dist_matrix(i,nn10_idcs));
    end
    disp('DONE - 10 nn')

    if filterOn
        disp('filtering...')
        %%% iteratively remove outlier values until log likelihood of normal curve fit is maximized
        [outlierness_sorted, sorted_idcs] = sort(outlierness, 'descend'); 
        pd = fitdist(outlierness_sorted,'Normal');
        stop = 0;
        prev_logL = pd.NLogL;   % negative log likelihood of probability distribution
        filtered_idcs = [];
        i = 1;

        while stop == 0 && i <= length(outlierness)   
            outlierness_sorted(1) = [];  % remove largest element
            pd = fitdist(outlierness_sorted,'Normal');
            if pd.NLogL >= prev_logL
                stop = 1;
            else
                prev_logL = pd.NLogL;
                filtered_idcs = [filtered_idcs sorted_idcs(i)];
                if length(outlierness_sorted) <= 2
                    stop = 1;
                    filtered_idcs = [];
                end
            end
            i = i + 1;
        end
        disp(['DONE - removed outliers (' num2str(length(filtered_idcs)) ')'])
    else
        filtered_idcs = [];
    end

    %%% remove outliers from data
    data_filtered = data;
    data_filtered(filtered_idcs,:) = [];

    good_idcs = setdiff(1:length(data), filtered_idcs);
    dist_matrix_filtered = dist_matrix(good_idcs, good_idcs);

    %classical multi-dimensional scaling
    dist_matrix_formatted = dist_matrix_filtered - 0.5; % the cdmscale function expects the distance matrix to contain 0 on the diagonal
    [Y_all,e_all] = cmdscale(dist_matrix_formatted);
    disp('DONE - multidimensional scaling')

    
    %size(Y_all,2) is the dimension of the smallest space in which the n points whose inter-point distances are given by D can be embedded.
    %the reduction to size(Y_all,2) or fewer dimensions provides a reasonable approximation to D only if the negative elements of e are small in magnitude.
    e_filtered = e_all(1:size(Y_all,2));
    
    Y = Y_all(:,1:dim);
%     explained = sum(e_filtered(1:dim))/sum(e_filtered);
%     neg_proportion = abs(sum(e_all(e_all < 0))) / abs(sum(e_all));
    
    if PCAon
        [~,~,~,~,explained_all] = pca(dist_matrix_formatted); 
        explained = sum(explained_all(1:dim)) / sum(explained_all);
    else
        explained = NaN;
    end
end