function [AggrResults, MorpResults] = PreProcess_Step1_Aggregation_v2(filename, options, figOpt)

% Inputs:
% filename.analfolder = string with path to analysis folder
% filename.resufolder = string with results subfolder if present
% filename.folders    = cell struct with string with tissues names used in naming the result files
% filename.tissues    = cell struct with string with tissues names usable for plots and such
% filename.suffix     = '.mat' usually
% options.maxround    = number of rounds in the experiments - not always cycles*4
% options.date        = used for name the results files uniquely if needed
% options.Markers     = cell struct with marker names for plot and such
% options.magnification = 20 if 20X, 4- if 40X etc


% This could/should be all contained in an excel file with 2 sheets called:
% - Folders, contains all the information to retrieve the input files and
% save the output files, column A is the name of the folder, basefolder,
% analfolder, resufolder, montfolder etc, column B is the string of folder
% - Settings, contains all the parameters used in the data analysis, with
% the space to input a value, a default value if the input value is left
% blank

Area = [];
CytArea = [];
Solidity = [];
MedianNucSign = [];
MedianCytSign = [];
MeanNucSign = [];
MeanCytSign = [];
Indexes = [];
X = [];
Y = [];
nuc = [];
cyt = [];

for fol = 1:length(filename.folders)
    clear Results  
    result_file = [filename.analfolder filename.resufolder filename.folders{fol} filename.suffix ];
    load(result_file)
    
    fol_num = uint8(zeros(length(Results.Area),1) + fol);
    
    Area            = [Area; Results.Area];
    CytArea         = [CytArea; Results.CytArea];
    Solidity        = [Solidity; Results.Solidity];
    
    MedianNucTemp = padarray(Results.MedianNucSign,[0 options.maxround-size(Results.MedianNucSign,2)],NaN,'post');
    MedianCytTemp = padarray(Results.MedianCytSign,[0 options.maxround-size(Results.MedianCytSign,2)],NaN,'post');
    
    MedianNucSign   = [MedianNucSign; MedianNucTemp];
    MedianCytSign   = [MedianCytSign; MedianCytTemp];
   
    MeanNucTemp = padarray(Results.MeanNucSign,[0 options.maxround-size(Results.MeanNucSign,2)],NaN,'post');
    MeanCytTemp = padarray(Results.MeanCytSign,[0 options.maxround-size(Results.MeanCytSign,2)],NaN,'post');
    
    MeanNucSign     = [MeanNucSign; MeanNucTemp];
    MeanCytSign     = [MeanCytSign; MeanCytTemp];

    Indexes         = [Indexes; fol_num];
    X               = [X; round(Results.CentroidX,0)];
    Y               = [Y; round(Results.CentroidY,0)];
    

    for j = 1:size(MedianNucTemp,2)
        [nuc{j}(fol,:),h_n{j}(fol,:)] = hist(log2(double(MedianNucTemp(:,j))+1),linspace(1,16,200));
        [cyt{j}(fol,:),h_c{j}(fol,:)] = hist(log2(double(MedianCytTemp(:,j))+1),linspace(1,16,200));
        nuc{j}(fol,:) = nuc{j}(fol,:)/sum(nuc{j}(fol,:));
        cyt{j}(fol,:) = cyt{j}(fol,:)/sum(cyt{j}(fol,:));
    end
    [n_morph{1}(fol,:),h_morph{1}(fol,:)] = hist(Results.Area,linspace(0,options.magnification*100,200));
    [n_morph{2}(fol,:),h_morph{2}(fol,:)] = hist(Results.CytArea,linspace(0,options.magnification*100,200));
    [n_morph{3}(fol,:),h_morph{3}(fol,:)] = hist(Results.Solidity,linspace(0,1,200));
    for k = 1:3
        n_morph{k}(fol,:) = n_morph{k}(fol,:)/sum(n_morph{k}(fol,:));
    end

    disp(['Done with ' filename.tissues{fol}])
end
disp('Saving files')
%%% make sure the variable are the correct data type to avoid space waste
AggrResults.MedianNucSign = uint16(MedianNucSign);
AggrResults.MedianCytSign = uint16(MedianCytSign);
AggrResults.MeanNucSign = uint16(MeanNucSign);
AggrResults.MeanCytSign = uint16(MeanCytSign);
AggrResults.Indexes = uint8(Indexes);
AggrResults.Area = uint16(Area);
AggrResults.CytArea = uint16(CytArea);
AggrResults.Solidity = double(Solidity);

MorpResults.Area = AggrResults.Area;
MorpResults.CytArea = AggrResults.CytArea;
MorpResults.Solidity = AggrResults.Solidity;
MorpResults.Indexes = AggrResults.Indexes;
MorpResults.X = uint16(X);
MorpResults.Y = uint16(Y);

%%% now plot the outputs
% first estimate how many subplots are needed
if length(filename.folders) < 6
    sub_rows = 2;
    sub_cols = 1;
    plotperfig = length(filename.folders);
    figpos = [0 0 0.5 0.5];
elseif length(filename.folders) < 26
    sub_rows = 2;
    sub_cols = ceil(length(filename.folders)/5);
    plotperfig = 5;
    figpos = [0 0 1 0.5];
elseif length(filename.folders) < 51
    sub_rows = 4;
    sub_cols = ceil(length(filename.folders)/10);
    plotperfig = 5;
    figpos = [0 0 1 1];
else
    sub_rows = 4;
    sub_cols = 5;
    plotperfig = ceil(length(filename.folders)/10);
    figpos = [0 0 1 1];
end

outplotdir= [filename.analfolder filename.resufolder 'Step1_Aggr_Plots\'];
mkdir(outplotdir)
disp('Plotting and saving outputs')
for i = 1:length(nuc)
    ax_lim = [min([h_n{i}(:); h_c{i}(:)]) max([h_n{i}(:); h_c{i}(:)]) min([nuc{i}(:); cyt{i}(:)]) max([nuc{i}(:); cyt{i}(:)])];
    figure(i)
    for j = 1:ceil(length(filename.folders)/plotperfig)
        start_fol = 1+(j-1)*plotperfig;
        end_fol = min(start_fol+plotperfig-1,length(filename.folders));
        % plot the nuclear signal
        subplot(sub_rows,sub_cols,j)
        plot(h_n{i}(start_fol:end_fol,:)',nuc{i}(start_fol:end_fol,:)')
        axis(ax_lim)
        title(['Nuclear ' options.Markers{i}])
        legend(filename.tissues{start_fol:end_fol},'Location','northeast')
        % plot the cytoplasmic signal
        subplot(sub_rows,sub_cols,j+sub_rows*sub_cols/2)
        plot(h_c{i}(start_fol:end_fol,:)',cyt{i}(start_fol:end_fol,:)')
        axis(ax_lim)  
        title(['Cytoplasmic ' options.Markers{i}])
        legend(filename.tissues{start_fol:end_fol},'Location','northeast')
    end
    set(gcf,'Units','Normalized','OuterPosition',figpos)
    % DA added 20200413 after Cyclin A1/2 errored on saving
    if contains(options.Markers{i},{'\','/',':','*','?','"','<','>','|'})
        saveas(figure(i),[outplotdir num2str(i) '_' ...
            options.Markers{i}(isstrprop(options.Markers{i}),'alpha')],'jpeg');
    else
        saveas(figure(i),[outplotdir, [num2str(i) '_' options.Markers{i}]],'jpeg');
    end
    %
    if figOpt == 0
        close all
    end
end

temp_title = {'Area','CytArea','Solidity'};
for i = 1:3
    figure(1000+i)
    if i ~= 3; xUP = prctile(n_morph{i}(:),90) + 1000;
    else xUP = 1; end
    for j = 1:ceil(length(filename.folders)/plotperfig)
        start_fol = 1+(j-1)*plotperfig;
        end_fol = min(start_fol+plotperfig-1,length(filename.folders));
        % plot the nuclear signal
        subplot(sub_rows/2,sub_cols,j)
        plot(h_morph{i}(start_fol:end_fol,:)',n_morph{i}(start_fol:end_fol,:)')
        title(temp_title{i})
        try
            legend(filename.tissues{start_fol:end_fol},'Location',options.Aggr.legendloc)
        catch
            legend(filename.tissues{start_fol:end_fol},'Location','best')
        end
        xlim([0 xUP]);
    end
set(gcf,'Units','Normalized','OuterPosition',[figpos(1:3) figpos(4)/2])
saveas(figure(1000+i),[outplotdir, 'CellShape_', temp_title{i}],'jpeg');
end

if figOpt == 0
    close all
end