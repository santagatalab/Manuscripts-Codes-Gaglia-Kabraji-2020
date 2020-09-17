% Data - MxN matrix of log2 intensity values
% Filter - MxN matrix of binaries
% channels - vector of channel indices we want to normalize ({1,...,N})
% channel_labels - 1xN cell array of channel names
% Params - object with the following attributes:
%   Params.Priors - 1xN vector, default 0, -1 = negative tail, +1 = positive tail
%   Params.CellNum - integer, number of cells to subsample (50,000 usually good)
%   Params.Reps - integer, number of times to repeat normalization (usu 5 is good)
%   Params.Overexpr - 1xN vector of doubles; predefines positive population (ignored if 0)
% FigSettings - object with the following attributes:
%   FigSettings.FigFlag - binary (1 = save figure)
%   FigSettings.Folder - string; filepath to save figures to
%   FigSettings.Debug - binary (1 = print figure)

function [DataNorm, cutoffs, multipliers] = norm_main(Data, Filter, channel_labels, Params)
    
    DataNorm = Data;
    cutoffs = zeros(size(Data,2),1);
    multipliers = zeros(size(Data,2),1);
    
    channels = Params.Channels;
    FigSettings = Params.FigSettings;
    mkdir(FigSettings.Folder)
    
    for i=1:length(channels)
        ch=channels(i);
        
        overexpr = Params.OverExpr(ch);
        final_est = 0;
        if isfield(Params, 'FinalEst')
            if length(Params.FinalEst) >= ch
                final_est = Params.FinalEst(ch);
            end
        end
        
        filter = Filter(:,ch) == 1;
        data_ch = Data(filter,ch);
        
        if Params.IsNuc
            label = [channel_labels{ch} '_nuc'];
        else
            label = [channel_labels{ch} '_cyt'];
        end
        
        [cutoff, mult] = norm_channel(data_ch, Params.CellNum, Params.Reps, Params.Priors(ch), ...
                                        overexpr, final_est, label, FigSettings);
        
        DataNorm(:,ch) = (Data(:,ch)-cutoff)/mult;
        cutoffs(ch) = cutoff;
        multipliers(ch) = mult;
    end
end