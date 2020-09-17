function [cutoffs] = visualize_cycle(dim, Y, data, all_channels, clust_channels, markers, cmap, condensedView)
    cutoffs = zeros(length(all_channels),1);
    if condensedView
        figure
    end
    if dim == 3
        for i=1:length(all_channels)
            ch = all_channels(i);
            channel_data = data(:,i);
            if ~condensedView
                figure
            else
                hold on
                subplot(4, ceil(length(all_channels)/4), i)
            end
            s1 = scatter3(Y(:,1), Y(:,2), Y(:,3), 5, channel_data, 'filled');
            s1.MarkerFaceAlpha = .3;

            cutoff = prctile(channel_data,95);
            if cutoff < 0
                cutoff = max(channel_data);
            end
            cutoffs(i) = cutoff;
            caxis([-cutoff cutoff])
            colormap(cmap)
            
            if ~condensedView
                colorbar
            end
            if ismember(ch, clust_channels)
                title(markers{ch}, 'Color', 'm')
            else
                title(markers{ch})
            end
        end
    elseif dim == 2
        for i=1:length(all_channels)
            ch = all_channels(i);
            channel_data = data(:,i);
            pos_idcs = find(channel_data > 0);
            neg_idcs = find(channel_data <= 0);

            if ~condensedView
                figure
            else
                hold on
                subplot(4, ceil(length(all_channels)/4), i)
            end
            s1 = scatter(Y(neg_idcs,1), Y(neg_idcs,2), 3, channel_data(neg_idcs), 'filled');
            s1.MarkerFaceAlpha = 1;
            hold on
            s2 = scatter(Y(pos_idcs,1), Y(pos_idcs,2), 3, channel_data(pos_idcs), 'filled');
            s2.MarkerFaceAlpha = 1;

            cutoff = prctile(channel_data,95);
            cutoffs(i) = cutoff;
            caxis([-abs(cutoff) abs(cutoff)])    
            colormap(cmap)
            if ~condensedView
                colorbar
            end

            if ismember(ch, clust_channels)
                title(markers{ch}, 'Color', 'm')
            else
                title(markers{ch})
            end
        end
    end
end