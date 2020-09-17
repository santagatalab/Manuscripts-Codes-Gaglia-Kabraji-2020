function [hier_clus,hier_idx,hier_C,hier_D,order] = cluster_cycIFdata_extradata(data,numclust,labels,figflag,iters,rep,ax_lim,cmap,image_title,extradata,extralabels)

% data = dataclus_1;
% labels= ColLabelsValue;
% figflag = 1;
% iters = 100;
% rep = 1;
% ax_lim = [-1.5 1.5];
% cmap = 'jet';
tic
opts = statset('Display','final');
[idx,C,~,D] = kmeans(data,numclust,'Distance','cityblock','MaxIter',iters, 'Replicates',rep,'Options',opts);


clus = zeros(max(idx),size(data,2));
clus_extra = zeros(max(idx),size(extradata,2));
toc
disp('Done with kmeans clustering')
for i = 1:max(idx)
    clus(i,:) = mean(data(find(idx==i),:),1);
    if size(extradata,1) ~= 0
%         clus_extra(i,:) = mean(extradata(find(idx==i),:),1);
%     else
%         clus_extra = [];
    end
end

T = clustergram(clus);
toc
disp('Done with hier clustering')

hier_idx = zeros(size(idx));
hier_C = zeros(size(C));
hier_D = zeros(size(D));

mapping_idx = str2num(cell2mat(T.RowLabels));
count = 0;
for i = 1:max(idx)
    hier_idx(idx==i,1)= find(mapping_idx == i);
    hier_C(find(mapping_idx == i),:) = C(i,:);
    hier_D(:,find(mapping_idx == i)) = D(:,i);
end

cellnum = length(idx);

hier_clus = zeros(max(idx),size(data,2));
hier_clus_extra = zeros(max(idx),size(extradata,2));
for i = 1:max(hier_idx)
    hier_clus(i,:) = mean(data(find(hier_idx==i),:),1);
    hier_clus_len(i) = sum(hier_idx==i);%/cellnum;
    if size(extradata,1) ~= 0
        hier_clus_extra(i,:) = mean(extradata(find(hier_idx==i),:),1);
%     else
%         hier_clus_extra = [];
    end 
end
toc
disp('Done with calculating cluster means')

median_dist= [];
dist = zeros(length(idx),2);

for i = 1:max(hier_idx)
    index_cells = hier_idx==i;
    median_dist(i,1) = median(hier_D(index_cells,i));
    median_dist(i,2) = prctile(hier_D(index_cells,i),75)-prctile(hier_D(index_cells,i),25);
    dist(index_cells,1) = hier_idx(index_cells);
    dist(index_cells,2) =  (hier_D(hier_idx(index_cells),i)-median_dist(i,1))/median_dist(i,2);
end


% dist = [];
% for i = 1:length(data(:,1))
%     dist(i,:) = [hier_idx(i) (hier_D(i,hier_idx(i))-median_dist(hier_idx(i),1))/median_dist(hier_idx(i),2)];
% end

toc
disp('Done with calculating dist')

% data_clus_l_ord = [dist(:,1)/numclust-(0.1*numclust/2)+dist(:,2)*0.000000001 dist(:,1)/numclust-(0.1*numclust/2) dist(:,2) data];
data_clus_l_ord = [dist(:,1) (hier_idx/numclust-0.5)*(ax_lim(2)-ax_lim(1))  data extradata];
[data_clus_l_ord,order] = sortrows(data_clus_l_ord,[2 1]);
data_clus_l_ord = data_clus_l_ord(:,2:end);

toc
disp('Done with reordering')

if figflag == 1
    figure
    subplot(1,3,1)
    imagesc(hier_clus_len')
    colormap(cmap)
    caxis([0 max(hier_clus_len)])
    for i = 1:max(hier_idx)
        text(1,i,num2str(hier_clus_len(i)),'VerticalAlignment','middle','HorizontalAlignment','center')
    end
    subplot(1,3,2)
    imagesc([hier_clus hier_clus_extra])
    colormap(cmap)
    xticks = linspace(1, numel([labels extralabels]), numel([labels extralabels]));
    set(gca, 'XTick', xticks, 'XTickLabel', [labels extralabels])
    caxis(ax_lim)    
    if class(image_title) == 'char'
        title(image_title)
    end
    subplot(1,3,3)
    imagesc(data_clus_l_ord)
    colormap(cmap)
    xticks = linspace(1, numel([labels extralabels]), numel([labels extralabels]));
    set(gca, 'XTick', xticks+1, 'XTickLabel', [labels extralabels])
    caxis(ax_lim)
end

ha=get(gcf,'children');
set(ha(1),'position',[0.52 0.05 0.47 0.9])
set(ha(2),'position',[0.06 0.05 0.45 0.9])
set(ha(3),'position',[0.01 0.05 0.03 0.9])

delete(T)