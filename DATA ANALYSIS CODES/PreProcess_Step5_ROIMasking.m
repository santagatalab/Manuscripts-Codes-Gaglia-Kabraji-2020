function ROIResults = PreProcess_Step5_ROIMasking(filename, options, ROIInfo, MorpResults)

% Inputs:
% filename.analfolder  = string with path to analysis folder
% filename.resufolder  = string with results subfolder if present
% filename.montfolder  = string with subfolder where montages are found
% filename.folders     = cell array of strings with tissues names
%                        {1 x t} where t is number of tissues
% filename.montsuffix  = string of suffix to montage images ('_montage.tif)
% options.pyramidlevel = from what ome level the masks were made (ie
%                        how much they are scaled. 0 is no scale, 1 is
%                        1/2 each dimension, 2 is 1/4 etc...)
% options.date         = date for saving the results file
% options.FigOpt       = logical to save overlays of masks on montages
% ROIInfo              = structure with groups of ROIs (eg: tumor, blood,
%                        spleen) each group is a cell array containing
%                        masks(either logical or bwlabeled) for each tissue
%                    Ex. ROIInfo.Tumor could be a {1 x t} cell array where
%                        each cell is a [p x q] logical or int. matrix.
%                        If the original tissue image is [m x n] pixels, p 
%                        and q are m/2^(pyramidlevel) and n/2^(pyramidlevel)
%                        respectively
%                    *Note: an empty cell {[]} & a mask of 0s are equivalent
%                    *Note: this structure is made in a step external to 
%                        this Pre-Processing pipeline. 
% MorpResults          = structure with indexes, and centroid (X,Y) vectors

% Outputs:
% ROIResults           = Structure containing label vectors for each ROI
%                        type (taken from ROIInfo). Vectors are [c x 1] and
%                        can be logical (eg spleen or not), or an integer
%                        (eg tumor grades on a scale of 1-4)
% if FigOpt == true,  the code will save images of each tissue where color
%                     ROIs are overlaid on the B&W montage for QC purposes


% initialize all the vectors:
ROITypes = fieldnames(ROIInfo);
if options.FigOpt
    mkdir([filename.analfolder filename.resufolder 'Step5_ROIMasks\' ])
end
ROINum = length(ROITypes);
for t = 1:ROINum
    ROIResults.([ROITypes{t} 'Index']) = zeros(size(MorpResults.X),'uint16');
end

for i = 1:length(filename.folders)
    % load in original image
    mont = imread([filename.analfolder filename.montfolder ...
        filename.folders{i} filename.montsuffix]);
    if options.FigOpt; montage = cat(3,mont, mont, mont); end
    [ysize,xsize]= size(mont);
    % first find all the cells that belong to the ROI
    index = MorpResults.Indexes == i;
    % get the centroids of the cells
    X_tissue = MorpResults.X(index);
    Y_tissue = MorpResults.Y(index);
    % rescale
    X_tissue = round(double(X_tissue)/(2^options.pyramidlevel),0);
    Y_tissue = round(double(Y_tissue)/(2^options.pyramidlevel),0);
    % make sure there are no zeros
    X_tissue(X_tissue==0) = 1;
    Y_tissue(Y_tissue==0) = 1;
    % find indices
    xy_index = sub2ind([ysize,xsize],Y_tissue,X_tissue);
    
    for t = 1:ROINum
        if ~isempty(ROIInfo.(ROITypes{t}){i})
        ROIResults.([ROITypes{t} 'Index'])(index) = ROIInfo.(ROITypes{t}){i}(xy_index);
        % making the check image
        if options.FigOpt
            c_max = 30000/min(max(ROIInfo.(ROITypes{t}){i}(:)),5);
            colors = {[c_max 0 0] [0 c_max 0] [0 0 c_max] [c_max c_max 0]...
                [c_max 0 c_max] [c_max c_max/2 0] [0 c_max/2 c_max/2]...
                [c_max 0 c_max/2] [0 c_max c_max] [c_max c_max c_max]};
                % supports 10 different ROI Types
                % red, green, blue, yellow, purple, orange, teal, magenta, cyan, white 
        montage(:,:,1) = montage(:,:,1)+uint16(ROIInfo.(ROITypes{t}){i}).*colors{t}(1);
        montage(:,:,2) = montage(:,:,2)+uint16(ROIInfo.(ROITypes{t}){i}).*colors{t}(2);
        montage(:,:,3) = montage(:,:,3)+uint16(ROIInfo.(ROITypes{t}){i}).*colors{t}(3);
        end
        end
    end
    if options.FigOpt
        saveas = [filename.analfolder filename.resufolder 'Step5_ROIMasks\' ...
                filename.folders{i} '_check.tif'];
        imwrite(montage,saveas,'tif')
    end
end
end