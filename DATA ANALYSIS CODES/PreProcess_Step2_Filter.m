function [outFilt, report] = PreProcess_Step2_Filter(Results, options, figOpt)

% Results must include:
    % Area = nuclear area (nx1)
    % CytArea = cytoplasmic area (nx1)
    % Solidity = solidity of the nucleus (nx1)
    % MedianNucSign = median nuclear signal in all the channels (nxc)
    % MedianCytSign = median cytoplasmic signal in all channels (nxc)
    % Indexes = tissue / ROI from which the cell comes from (nx1)

% options:
    % maxround = max number of channels (*** only required option ***)
    % folder = path to save output images (default is current directory)
    % Index_Names = names of the tissues / ROIs, in order
% options.thresholds can include (eg options.thresholds.solidity):
    % foldDAPI_th = ratio of cytoplasmic to nuclear DAPI
    % absDAPI_th = max DAPI signal (log2)
    % solidity = nuclear integrity (0-1 scale)
	% area_low = min nuclear area
	% area_high = max nuclear area
    % cytarea_low = min cytoplasmic area
    % cytarea_high = max cytoplasmic area
% any option.thresholds left blank or omitted will not be filtered with

% outputs:
    % outFilt = structure with 1 logical filter for each option.thresholds,
              % 1 overall logical filter
    % report = cell array with percentages of cells lost to a particular
              % filter in each round
    % figures = (if figOpt, these are kept open, saved in options.folder) 
              % Cell Loss (red) visualization in subplots:
              % 1, (50k subset) lost to N vs C DAPI
              % 2, cells lost in each cycle (due to abs_DAPI if
                   % applicable, or overall if not)
              % 3, cells lost to solidity
              % 4, (50k subset) lost to area constraints

% step 1: determine which options will be filtered on
filters = zeros(1,7);
filt_names = {'all','DAPI_Fold','DAPI_Thresh','Solidity','NucArea_Low',...
    'NucArea_High', 'CytArea_Low', 'CytArea_High'};
rej = []; % this vector signified which options will *not* be used to filter
try filters(1) = options.thresholds.foldDAPI_th;
catch filters(1) = 100; rej = [rej,1];                          %#ok<*SEPEX>
end
try filters(2) = options.thresholds.absDAPI_th;
catch filters(2) = 100; rej = [rej,2];
end
try filters(3) = options.thresholds.solidity;
catch filters(3) = 100; rej = [rej,3];
end
try filters(4) = options.thresholds.area_low;
catch filters(4) = 100; rej = [rej,4];
end
try filters(5) = options.thresholds.area_high;
catch filters(5) = 100; rej = [rej,5];
end
try filters(6) = options.thresholds.cytarea_low;
catch filters(6) = 100; rej = [rej,6];
end
try filters(7) = options.thresholds.cytarea_high;
catch filters(7) = 100; rej = [rej,7];
end
optNum = 7-length(rej);
filt_names(rej+1) = [];
if length(rej) == 7
    disp('Warning: No filters set')
end

% step 2: determine number of rounds / samples for looping
try cyc = options.maxround./4;
catch
    disp('Warning: Non-integer cycle number')
    cyc = ceil(options.maxround./4);
end
samp = unique(Results.Indexes);
try names = length(options.Index_Names) == length(samp);
catch names = false;
end

% step 3: filtering
Filter = ones(size(Results.MedianNucSign),'logical');
outCell = repmat({Filter},1,optNum+1);
outFilt = cell2struct(outCell,filt_names,2);
report = cell(10,cyc+1);
for c = 1:cyc 
    DAPI_col = (c-1)*4+1;

    reject{2} = Results.MedianNucSign(:,DAPI_col)./Results.MedianCytSign(:,DAPI_col)<filters(1);
    reject{3} = Results.MedianNucSign(:,DAPI_col)<(2^filters(2));
    reject{4} = Results.Solidity(:) < filters(3);
    reject{5} = Results.Area(:) < filters(4);
    reject{6} = Results.Area(:) > filters(5);
    reject{7} = Results.CytArea(:) < filters(6);
    reject{8} = Results.CytArea(:) > filters(7);
    reject{1} = zeros(size(Filter,1),1);
        
    for f = 2:8
       if ~ismember(f-1,rej) 
           reject{1} = reject{1} | reject{f};
           tempFilt = outFilt.(filt_names{f});
           tempFilt(reject{f},DAPI_col:(DAPI_col+3)) = false;  
           outFilt.(filt_names{f}) = tempFilt;
       end
    end
    tempFilt = outFilt.(filt_names{1});
    tempFilt(reject{1},DAPI_col:(DAPI_col+3)) = false;  
    outFilt.(filt_names{1}) = tempFilt;
    
    totcells = size(Results.MedianNucSign,1);    
    report(:,1) = {[],'cell number', 'all rej', 'DAPI fold rej', ...
        'DAPI absolute rej','solidity rej', 'nuc area high rej', ...
        'nuc area low rej', 'cyt area high rej', 'cyt area high rej'} ;
    report(:,c+1) = {['Cycle ' num2str(c)], ...
                     totcells, ...
                     sum(reject{1})/totcells,  ...           
                     sum(reject{2})/totcells , ...      
                     sum(reject{3})/totcells , ...    
                     sum(reject{4})/totcells , ...     
                     sum(reject{5})/totcells , ...    
                     sum(reject{6})/totcells , ...
                     sum(reject{7})/totcells , ...
                     sum(reject{8})/totcells};     
    report(rej+3,c+1) = {NaN};
end

% Step 4: Visualizing
if figOpt
    view = 'on';
else
    view = 'off';
end

    try outfol = [options.folder 'FilterPrints\'];
        mkdir(outfol)
    catch outfol = [pwd '\FilterPrints\'];
        mkdir(outfol)
    end
    try savefile = [outfol options.date 'report.mat'];
    catch savefile = [outfol 'report.mat'];
    end
    save(savefile,'report')
    
    %'fold','abs','sol','n_area_l','n_area_h','c_area_l','c_area_h'
    canfilt = logical([1 1 1 1 1 1 1]); canfilt(rej) = false;
    plot_filter = zeros(size(outCell{1},1),7,'logical');
    for f = 1:7
        if canfilt(f)
            plot_filter(:,f) = outFilt.(filt_names{1+sum(canfilt(1:f))})(:,1);
        else
            plot_filter(:,f) = ones(size(outCell{1},1),1,'logical');
        end
    end
    if canfilt(2)
        thisfilt = {true, outFilt.(filt_names{1+sum(canfilt(1:2))})};
    else
        thisfilt = {false, outFilt.(filt_names{1})};
    end
    Ind = ones(size(plot_filter,1),1,'logical');
    f = figure(1);
    f.Visible = view;
    plothelp(Results, plot_filter, thisfilt, Ind, cyc)
    suptitle('Cell Loss overall')
    saveas(f,[outfol, 'Overall'],'jpeg');
    if ~figOpt
        close all
    end

    % now split up by Indexes
    if length(samp) > 1
    t = 1;
    for i = 1:length(samp)
        indices = Results.Indexes == samp(i);
        f1 = figure(100+i);
        f1.Visible = view;
        plothelp(Results, plot_filter, thisfilt, indices, cyc)
        if names
            title = ['Cell Loss ' options.Index_Names{t}];
            suptitle(title)
            t = t+1;
        else
            title = ['Cell Loss ' num2str(samp(i))];
            suptitle(title)
        end
        saveas(f1,[outfol, title],'jpeg');
        if ~figOpt
            close all
        end
    end
    end
    
end

function plothelp(Results, plot_filter, thisfilt, Ind, cyc)

if sum(Ind) > 50000
    y = randsample(length(Ind),50000);
    reduce = zeros(size(Ind),'logical'); 
    reduce(y) = true;
else
    reduce = ones(size(Ind),'logical');
end

    subplot(2,2,1)
    scatter(log2(double(Results.MedianNucSign(Ind&reduce&plot_filter(:,1)&plot_filter(:,2),1)+1)),...
        log2(double(Results.MedianCytSign(Ind&reduce&plot_filter(:,1)&plot_filter(:,2),1)+1)),...
        5,'filled','MarkerFaceColor','b','MarkerFaceAlpha',0.5)
    hold on
    scatter(log2(double(Results.MedianNucSign(Ind&reduce&~plot_filter(:,1),1)+1)),...
        log2(double(Results.MedianCytSign(Ind&reduce&~plot_filter(:,1),1)+1)),...
        5,'filled','MarkerFaceColor','r','MarkerFaceAlpha',0.75)
    scatter(log2(double(Results.MedianNucSign(Ind&reduce&~plot_filter(:,2),1)+1)),...
        log2(double(Results.MedianCytSign(Ind&reduce&~plot_filter(:,2),1)+1)),...
        5,'filled','MarkerFaceColor','r','MarkerFaceAlpha',0.5)
    xlabel('Nuclear DAPI (log2)') 
    ylabel('Cytoplasmic DAPI (log2)')
    axis tight 
    
    subplot(2,2,2)
    X = 1:cyc;
    tempY = (X-1).*4+1;
    Y1 = sum(thisfilt{2}(Ind,tempY),1)';
    Y2 = sum(~thisfilt{2}(Ind,tempY),1)';
    Y = [Y1,Y2];
    h = bar(X,Y,'stacked','b');
    h(1).FaceColor = 'b';
    h(2).FaceColor = 'r';
    xlabel('Round')
    ylabel('Cell Number')
    if thisfilt{1}
        title('Cell Loss Due To DAPI Threshold')
    else
        title('Cell Loss Overall')
    end
    
    subplot(2,2,3)
    X = [sum(Results.Solidity(Ind) > 0.8), sum(Results.Solidity(Ind) > 0.6 & Results.Solidity(Ind) <= 0.8),...
        sum(Results.Solidity(Ind) > 0.4 & Results.Solidity(Ind) <= 0.6), ...
        sum(Results.Solidity(Ind) > 0.2 & Results.Solidity(Ind) <= 0.4),...
        sum(Results.Solidity(Ind) < 0.2)];
    mask = X < 10;
    X(mask) = [];
    labels = {'0.8 < S', '0.6 < S <= 0.8','0.4 < S <= 0.6','0.2 < S <= 0.4', 'S < 0.2'};
    labels(mask) = [];
    p = pie(X);
    legend(labels,'Location','best')
    for i = 1:2:length(X)*2
        choices = {'b', 'b','g', 'g', 'y','y','r', 'r','m', 'm'};
        p(i).FaceColor = choices{i};
    end
    ylabel('Solidity')
    
    subplot(2,2,4)
    scatter(Results.Area(Ind&reduce&plot_filter(:,4)&plot_filter(:,5)&plot_filter(:,6)&plot_filter(:,7)),...
        Results.CytArea(Ind&reduce&plot_filter(:,4)&plot_filter(:,5)&plot_filter(:,6)&plot_filter(:,7)),...
        5,'filled','MarkerFaceColor','b','MarkerFaceAlpha',0.5)
        % ^ plots all cells that are not filtered out
        hold on
    scatter(Results.Area(Ind&reduce&plot_filter(:,5)&~plot_filter(:,7)),...
        Results.CytArea(Ind&reduce&plot_filter(:,5)&~plot_filter(:,7)),...
        5,'filled','MarkerFaceColor','r','MarkerFaceAlpha',0.5)
        % ^ plots cells that are above the cyt thresh
    scatter(Results.Area(Ind&reduce&~plot_filter(:,4)&plot_filter(:,7)),...
        Results.CytArea(Ind&reduce&~plot_filter(:,4)&plot_filter(:,7)),...
        5,'filled','MarkerFaceColor','r','MarkerFaceAlpha',0.5)
        % ^ plots cells that are below the nuc thresh
    scatter(Results.Area(Ind&reduce&plot_filter(:,4)&~plot_filter(:,6)),...
        Results.CytArea(Ind&reduce&plot_filter(:,4)&~plot_filter(:,6)),...
        5,'filled','MarkerFaceColor','r','MarkerFaceAlpha',0.5)
        % ^ plots cells that are below the cyt thresh
    scatter(Results.Area(Ind&reduce&~plot_filter(:,5)&plot_filter(:,6)),...
        Results.CytArea(Ind&reduce&~plot_filter(:,5)&plot_filter(:,6)),...
        5,'filled','MarkerFaceColor','r','MarkerFaceAlpha',0.5)
        % ^ plots cells that are above the nuc thresh
    ylabel('Cytoplasmic Area')
    xlabel('Nuclear Area')
    axis tight 
end

