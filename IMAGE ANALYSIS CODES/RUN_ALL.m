clear all
%
% t-CycIF ASHLAR ANALYSIS PIPELINE 
% Step 0: Stitch and register the fields into an "ome.tif" by using Ashlar 
% before starting this workflow  
% Step 1:(RUN_Step1_new_fields_preilastik_crops.m) Cut the "ome.tif" into 
% fields of a size specified by the user and cut crops out of each field 
% for Ilastik training. Create full stacks of all cycles and channels. 
% Omits fields that have no cells.   
% Step 2: Create segmentation probabilities using Ilastik 
% Step 3:(RUN_Step3_segmentfromilastik.m) Segment based on segmentation
% probabilities produce by Ilastik
% Step 4: (RUN_Step4_CycIF_measurments_ilastik.m) Makes measurements of
% signal and foci segmentation 
% Step 5: (RUN_Step5_ROI.m) creates small montage for annotations (e.g. in
% ImageJ)
% 
% In order to begin the analysis a set of parameters are needed to be
% defined by the user 
% 1) where the files are located and how the names are formatted
% 2) parameters to choose the size of the field and number of crops per
% field 

%%% OUTPUTS FILES AND LOCATIONS
% Starting Inputs - ome.tif file folder = filename.folders.ashlared
% choose the parameters below to define folder location and image analysis
% settings

% Step 1: outputs the following folders and files
% FullStacks: ANALYSIS\FullStacks\
% CroppedStacks: AllTheRawData\OmeTifName\CroppedData\OmeTifName_Field_row_column_#ofCrop.tif
% Coordinates: AllTheRawData\Coordinates_FullStacks\OmeTifName.mat 
% y= # of pixels in 1 column; x = # of pixels in 1 row ; t = size of field desired 
% Coordinates Matrix: Coordinates.Field(row,column) = [keep field (0= no, 1=yes), x1, x2, y1, y2]
% Step 2b (optional): eliminates crops tiles that do not have cells
% Step 3:
% Segmented Images: ANALYSIS\Ilastik_Segmentation\OmeTifName_Field_row_column_Seg.tif
% Check Segmented Images: ANALYSIS\Ilastik_Segmentation\OmeTifName_Field_row_column_checkseg.tif
% Step 4: 
% Nucleus & Cytoplasm SegmentedImages: ANALYSIS\Ilastik_Segmentation\OmeTifName_Field_row_column_NucCytSeg.tif
% Measurments file: ANALYSIS\Analysis_Results\Results_data.mat 
% Step 5 (optional): output montage file in ANALYSIS\MontageforROI_LVxx

%%% INPUT FILE LOCATION AND FORMATTING
%
% The code will run on the Ilastik Probabilities file and FullStacks file 
% The expected file structure is that a master folder will contain all of
% the "ome.tif" data and within that folder, a folder called "ANALYSIS\"
% will contain all the Ilastik Probabilities and FullStacks data 
%
%%% IMPORTANT: Ilastik Output Filename Format for Export Options
% {dataset_dir}/Ilastik_Probabilities/{nickname}_Probabilities.tif

% eg. AllTheRawData\ANALYSIS\FullStacks\Ilastik_Probabilities\SlideName_Field_row_column_Probabilities.tif 
% eg. AllTheRawData\ANALYSIS\FullStacks\SlideName_Field_row_column.tif 


% 1) THE MASTER OME.TIF DATA FOLDER
filename = [];
filename.folders.main = 'Z:\addyourfoldernamehere\';
% Z:\sorger\
% Z:\data\IN_Cell_Analyzer_6000\Danae
% Y:\Danae

% 2) SLIDE SPECIFIC PARAMETERS:
filename.tissues = {'tissue1','tissue2','otherrtissuename'};
                %Name of Ashlared image without '.ome.tif' ending 
                                                 
filename.cycles = 8; % # of cycles to analyse
filename.maxround = filename.cycles*4;
filename.ilastikround = 12; % these are the rounds to be used to segment in Ilastik

% 3) USER DESIRED PARAMETERS 
filename.sizefield = 5000; %Size of field desired for a square 
filename.crops = 1; %# of cropped fields desired per field 
options.DAPIbkgd = 100; % thresh for DAPI used to check if we keep the field
options.DAPIbkgd_crops75prctile = options.DAPIbkgd*5;

% 4) OPTIONS  
% Step 3: SEGMENTATION OPTIONS 
options.nuc = 1;    %Channel nucleus Ilastik probability is in
options.cyt = 2; %Channel cytoplasm Ilastik probability is in
options.backgr = 3; %Channel background Ilastik probability is in 
options.cellsize = 9;
options.bkgdprob_min = 0.75;
options.cytprob_min  = 0.75;
options.nucprob_min  = 0.15;
options.max_prob = 65535; %Maximum Ilastik probability, usually does not change
options.pausetime = 0;

% Step 4: MEASUREMENTS AND FOCI OPTIONS 
options.date = 'XXX';  % add the date here

% Step 5 which slice to use for the montage for ROI
options.pyramidlevel = 3;



% 5) OUTPUT PARAMETERS: DO NOT EDIT 
filename.folders.ashlared = 'DATA\ashlar\';   
filename.ashlarsubfol = 'registration\';
filename.ashlarsuffix = '.ome.tif';
filename.folders.output = 'ANALYSIS\'; 
filename.folders.fullstacks = 'FullStacks\';
filename.folders.coordinates = 'Coordinates_FullStacks\';
filename.folders.ilastikstacks = 'IlastikStacks\';
filename.folders.cropfol = 'CroppedData\';
filename.folders.ilastikprob= 'IlastikStacks\';
filename.folders.ilastikseg = 'Ilastik_Segmentation\';
filename.ilastiksuffix = '_Probabilities.tif';

% filename.folders.ometiff = 'DATA\ashlar\';
filename.folders.montagelowres = 'MontageforROI_Lv3\';
filename.folders.results = 'Analysis_Results\';
filename.folders.ROI = 'ROIs\'; 
filename.suffix = '.tif';
% filename.folders.montage = 'MontageforROI\'; 
filename.dim = ['%02d']; %Delimeter 
filename.overwrite_seg = 0;

%% MAKE SURE TO COMMENT OUT THE OTHER STEPS THAT YOU AREN'T AT YET!!
% Runs Step 1: (Creating FullStacks, CroppedStacks, Montages) 

Step1_fields_preilastik(filename,options) 

% after this step you need to:
% - change the slices to channel eg using the ImageJ macro RUN_Step2_changetifffromZtoXhannel.ijm
% - perform Ilastik training and batch processing to create probability maps

%% clean up the crops to avoid useless ones

Step2b_filtercrops(filename,options)
%% Runs Step 3 and 4: (Segmentation and Measurments)

Step3_segmentfromilastik(filename, options) 
Step4_CycIF_measurements_ilastik(filename, options) 

%% Extract montages for ROI - if applicable
Step5_montageforROI_v2(filename, options)



