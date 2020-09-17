clear all
%%%
% codedir = 'Z:\data\IN_Cell_Analyzer_6000\Giorgio\CycIF Codes\Utility Functions';
% addpath(codedir)
% addpath('...\normalization\functions\')
    % ^ where the normalization codes are found
%%%
filename.basefolder = 'Y:\add_analysys_folder_here\';
filename.suffix = '_Results.mat';
filename.analfolder = [filename.basefolder 'ANALYSIS\'];
filename.resufolder = 'Results_singletissue\'; 
filename.roifolder  = 'ROIs\';
filename.montfolder = 'TumorMasks\MontageforROI_Lv3\';
filename.montsuffix = '_montage.tif';

filename.folders = {'MLB_266_1',      'MLB_267_6',      'MLB_268_2',     ...
                    'MLB_269_5_ROI1', 'MLB_269_5_ROI2', 'MLB_270_1',     ...
                    'MLB_271_3',      'MLB_272_2_ROI1', 'MLB_272_2_ROI2',...
                    'MLB_273_4_ROI1', 'MLB_273_4_ROI2', 'MLB_274_3',     ...    
                    'MLB_275_1',      'MLB_276_2',      'MLB_277_1',     ...
                    'MLB_278_4_ROI1', 'MLB_278_4_ROI2', 'MLB_279_4_ROI1',...
                    'MLB_279_4_ROI2', 'MLB_280_2_ROI1', 'MLB_280_2_ROI2',...
                    'MLB_281_3_ROI1', 'MLB_281_3_ROI2', 'MLB_282_4_ROI1',...
                    'MLB_282_4_ROI2', 'MLB_283_1',      'MLB_284_2',     ...
                    'MLB_285_6' };   
                      
for i = 1:length(filename.folders)
    options.MouseNum(i) = str2num(filename.folders{i}(5:7));
    options.MouseGroup(i) = ceil((options.MouseNum(i)-265)/5);
    if length(filename.folders{i}) == 9
        filename.tissues{i} = [filename.folders{i}(1:3) filename.folders{i}(5:7)];
    else
        filename.tissues{i} = [filename.folders{i}(1:3) filename.folders{i}(5:7) filename.folders{i}(11:14)];
    end
end

options.Markers =  {  'DAPI 0','bk 488'  ,'bk 555' ,'Tim3-1', ...
                      'DAPI 1','TTF1'    ,'B220'   ,'CD45', ...
                      'DAPI 2','FOXP3'   ,'CD4'    ,'CD8-1',...
                      'DAPI 3','CD103'  ,'CD11c'   ,'CD11b',...
                      'DAPI 4','NKp46'  ,'CD3e'    ,'Ki67', ...
                      'DAPI 5','PD-L1'  ,'PD-1'    ,'CD8-2', ...
                      'DAPI 6','GranzB' ,'Perforin','HSF1', ...
                      'DAPI 7','Tim3-2'  ,'Ly6g'    ,'bkgd' ...
                        }; % list of mutliplex markers in order 
                    
options.maxround = length(options.Markers);
options.magnification = 20; % either 20X or 40X
options.FigOpt = 0;
options.date = '20200408';

% STEP 2: additional parameters for filtering          
options.Filtering.folder = [filename.analfolder filename.resufolder 'Step2_'];
options.Filtering.Index_Names = filename.tissues;
options.Filtering.thresholds.foldDAPI_th = 0.5;
options.Filtering.thresholds.absDAPI_th = 9;
options.Filtering.thresholds.solidity = 0.7;
options.Filtering.thresholds.area_low = 20;    % was 50   in previous analysis
options.Filtering.thresholds.area_high = 500;  % was 2000 in previous analysis..
options.Filtering.maxround = options.maxround;

% STEP 3: additional parameters for normalization
options.Norm.Reps = 5;
options.Norm.FigSettings.FigFlag = 1; % to save
options.Norm.FigSettings.Folder = [filename.analfolder filename.resufolder 'Step3_NormPrints\' ];
options.Norm.FigSettings.Debug = 1; % to view
options.Norm.Channels = setdiff(1:length(options.Markers), 1:4:length(options.Markers));  % non-dapi channels
% default to zeros
options.Norm.Priors = zeros(length(options.Markers),1); 
options.Norm.OverExpr = zeros(length(options.Markers),1);
options.Norm.CellNum = min([50000, length(AggrResults.MedianNucSign)]);

% STEP 4: roi extracting
options.tifnames.Aiforia = {'TumorGrade1.tif','TumorGrade2.tif','TumorGrade3.tif'};
options.tifnames.Path = {'Tumors.tif','Blood.tif','Spleen.tif'};
options.ROI.pyramidlevel = 4;

% STEP 5: ROI Masking
options.ROI.pyramidlevel = 2; % where full size is 0
options.ROI.FigOpt = true;
options.ROI.date = '20200424';                 

save([filename.analfolder filename.resufolder 'Results_Settings_' options.date '.mat'],'filename','options')

%% Step 1 aggregate tissues together
clearvars -except filename options
[AggrResults, MorpResults] = PreProcess_Step1_Aggregation_v2(filename, options, options.FigOpt);

save([filename.analfolder filename.resufolder 'Results_Aggr_' options.date '.mat'],'AggrResults');
save([filename.analfolder filename.resufolder 'Results_Morp_' options.date '.mat'],'MorpResults');

%% Step 2 - Filter the data

clearvars -except filename options
load([filename.analfolder filename.resufolder 'Results_Morp_' options.date '.mat'])
load([filename.analfolder filename.resufolder 'Results_Aggr_' options.date '.mat'])
options.FigOpt = 0;

[Filter,report] = PreProcess_Step2_Filter(AggrResults, options.Filtering, options.FigOpt);
save([filename.analfolder filename.resufolder 'Results_Filt_' options.date '.mat'],'Filter','report')

%%
% close all
clearvars -except filename options
load([filename.analfolder filename.resufolder 'Results_Filt_' options.date '.mat'])
load([filename.analfolder filename.resufolder 'Results_Aggr_' options.date '.mat'])
load([filename.analfolder filename.resufolder 'Results_Morp_' options.date '.mat'])

rng(11)
close all

filter = Filter.all;
data_nuc = log2(double(AggrResults.MedianNucSign)+1);
data_cyt = log2(double(AggrResults.MedianCytSign)+1);

% normalize
options.Norm.IsNuc = 1;
[NormResults.MedianNucNorm, cutoffs_nuc, mults_nuc] = norm_main(data_nuc, filter, options.Markers, options.Norm);
options.Norm.IsNuc = 0;
[NormResults.MedianCytNorm, cutoffs_cyt, mults_cyt] = norm_main(data_cyt, filter, options.Markers, options.Norm);

%save
redmap = [linspace(0,255,128) zeros(1,128)+255 ]/255;
blumap = [zeros(1,128)+255 flip(linspace(0,255,128))]/255;
gremap = [linspace(128,255,128) flip(linspace(128,255,128))]/255;
NormResults.colorMap = [redmap' gremap' blumap'];
NormResults.CellID = (1:size(NormResults.MedianNucNorm,1))';

NormResults.MedianNucNorm = int16(round(NormResults.MedianNucNorm*1000,0));
NormResults.MedianCytNorm = int16(round(NormResults.MedianCytNorm*1000,0));
NormResults.nuc_add_fact = int16(1000.*cutoffs_nuc); 
NormResults.nuc_mult_fact = int16(1000.*mults_nuc);
NormResults.cyt_add_fact = int16(1000.*cutoffs_cyt); 
NormResults.cyt_mult_fact = int16(1000.*mults_cyt);

save([filename.analfolder filename.resufolder 'Results_Norm_' options.date '.mat'],'NormResults','options','filename');
save([filename.analfolder filename.resufolder 'Results_Settings_' options.date '.mat'],'filename','options');

%% Step 4: use ROI montages to extract annotations from other sources
% this step depends on where the masks come from and will need to output a
% structure with all the mask types (ie tumor, blood, spleen) containing
% binary or bwlabeled masks that fit perfectly over the original tissue
clearvars -except filename options

for i = 1:length(filename.folders)
    ROI.Aiforia{i} = [];
    for grade = 1:length(options.tifnames.Aiforia)
        Grade{grade} = uint16(imread([filename.analfolder filename.roifolder filename.folders{i}(5:end) '_' options.tifnames.Aiforia{grade}]));
      
        if sum(ROI.Aiforia{i}) == 0
            ROI.Aiforia{i} = Grade{grade}*grade;
        else
            ROI.Aiforia{i} = ROI.Aiforia{i} + Grade{grade}*grade;
        end
    end
    if options.FigOpt == 1
        figure
        imshow(ROI.Aiforia{i},[])
        title(filename.folders{i}(5:end))
    end
    ROI.Tumors{i} = [];
    ROI.Spleen{i} = [];
    ROI.Blood{i} = [];
    for sc = 1:length(options.tifnames.Path)
        Score{sc} = int16(imread([filename.analfolder filename.roifolder filename.folders{i}(5:end) '_' options.tifnames.Path{sc}]));
        if strcmp(options.tifnames.Path{sc}(1:end-4),'Tumors')
            ROI.Tumors{i} = uint16(Score{sc});
        elseif strcmp(options.tifnames.Path{sc}(1:end-4),'Blood')
            ROI.Blood{i} = uint16(Score{sc});
        elseif strcmp(options.tifnames.Path{sc}(1:end-4),'Spleen')
            ROI.Spleen{i} = uint16(Score{sc});
        end
    end
    if options.FigOpt == 1
        figure
        imshow(ROI.Tumors{i}-ROI.Blood{i}-ROI.Spleen{i}*2,[-2 10])
        title(filename.folders{i}(5:end))
    end
end
save([filename.analfolder filename.resufolder 'Results_ROI_' options.date '.mat'],'ROI')

%% Step 5: find whether cells belong to ROIs
clearvars -except filename options
load([filename.analfolder filename.resufolder 'Results_Morp_' options.date '.mat'])
load([filename.analfolder filename.resufolder 'Results_ROI_' options.date '.mat'])
load([filename.analfolder filename.resufolder 'Results_Settings_' options.date '.mat'])

MorpInput.X = MorpResults.X; MorpInput.Y = MorpResults.Y;
MorpInput.Indexes = MorpResults.Indexes; clearvars MorpResults

ROIResults = PreProcess_Step5_ROIMasking(filename, options.ROI, ROI, MorpInput);

save([filename.analfolder filename.resufolder 'Results_ROI_' options.date '.mat'],'ROIResults','-append')













