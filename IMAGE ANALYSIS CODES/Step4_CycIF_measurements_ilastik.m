function Step4_CycIF_measurements_ilastik(filename, options)  



addpath(filename.folders.main)
mkdir([filename.folders.main  filename.folders.output filename.folders.results]); 
tic

for k = 1:length(filename.tissues)
   
    
    fold = filename.tissues{k};
    filename.folders.resultfile = [filename.folders.results fold '_Results_' options.date '.mat'];
    
    addpath(filename.folders.main)
    
    %%%% get coordinates files for centroid rebuilding
    coordmat = [filename.folders.main filename.folders.output filename.folders.coordinates fold '.mat'];
    load(coordmat, 'Coordinates')
    [r, c] = size(Coordinates.Field);
    
    %%%% initialize the variables we are going to measure
    
    if exist([filename.folders.main filename.folders.output filename.folders.resultfile],'file') == 2
        load([filename.folders.main filename.folders.output filename.folders.resultfile])
    else
        Results = [];
        Results.FieldCoord = [];
        Results.FieldFlag = [];
        Results.Area = [];
        Results.Solidity = [];
        Results.CentroidX = [];
        Results.CentroidY = [];
        Results.MedianNucSign = [];
        Results.MedianCytSign = [];
        Results.MeanNucSign = [];
        Results.MeanCytSign = [];
    end
  
    %%%% loop around the split fields
    for i1 = 1:r
        for i2 = 1:c
            
            %%%% load the FullStack and/to check it exists ----------------
            FieldFile = [fold '_Field_' num2str(i1 , filename.dim) '_' num2str(i2 , filename.dim)];
            FileTif = [filename.folders.main filename.folders.output filename.folders.fullstacks FieldFile  filename.suffix];
            
            if exist(FileTif,'file')~=2
                disp(FileTif)
                disp('Error: FullStack not found')
                continue
            end
            disp(FieldFile)
            
            %%%% load nuclear mask, create cytoplasmic mask and save them -
            segfile = [filename.folders.main filename.folders.output filename.folders.ilastikseg FieldFile '_Seg.tif'];
            cyt_segfile = [filename.folders.main filename.folders.output filename.folders.ilastikseg FieldFile '_NucCytSeg.tif'];
            
            if exist(segfile,'file') ~= 2
                disp('Segmentation file not found')
                continue
            end
            
            if size(Results.FieldFlag,1) >= i1 && size(Results.FieldFlag,2) >= i2
                if Results.FieldFlag(i1,i2) == 1
                    disp('Field already analysed')
                    continue
                end
            end
            
            
            NucMask = uint16(imread(segfile));
            if max(NucMask(:)) == 1
                lb_NucMask = uint32(bwlabel(NucMask));
            else
                disp('Seg file already labelled, right?')
                lb_NucMask = NucMask;
            end
            % create cytoplasmic mask
            lb_CytMask = imdilate(lb_NucMask,offsetstrel('ball',ceil(options.cellsize/3),0));
            lb_CytMask(lb_NucMask>0)=0;
            CytMask = uint16(lb_CytMask > 0);
            DAPI_img = uint16(imread(FileTif,'Index',1));
            
            NucCytSeg = cat(3,NucMask,DAPI_img,CytMask);
            imwrite(NucCytSeg,cyt_segfile);

            %%%% measure cell properties and assign them to results -------
            stats_nuc = regionprops(lb_NucMask,'Area','Centroid','Solidity');
            stats_cyt = regionprops(lb_CytMask,'Area');

            Area = {stats_nuc.Area};
            Solidity = {stats_nuc.Solidity};
            Centroid = {stats_nuc.Centroid};
            CentroidVect = cell2mat(Centroid);
            CentroidMat = reshape(CentroidVect,2,length(CentroidVect)/2);
            CytArea = {stats_cyt.Area};
            
            totcells_field = length(Area);
            currentcells = length(Results.Area);
            startcell = currentcells + 1;
            endcell = currentcells + totcells_field;
            
            if totcells_field < 1
                disp(['No cells in Field ' FieldFile])
                continue
            end

            Results.Area(startcell:endcell,1) = cell2mat(Area)';
            Results.Solidity(startcell:endcell,1) = cell2mat(Solidity);
            Results.CentroidX(startcell:endcell,1) = CentroidMat(1,:)+Coordinates.Field{i1,i2}(4);
            Results.CentroidY(startcell:endcell,1) = CentroidMat(2,:)+Coordinates.Field{i1,i2}(2);
            Results.CytArea(startcell:endcell,1) = [cell2mat(CytArea)'; zeros(totcells_field-length(cell2mat(CytArea)),1)];
            Results.FieldCoord(startcell:endcell,1) = zeros(totcells_field,1)+i1;
            Results.FieldCoord(startcell:endcell,2) = zeros(totcells_field,1)+i2;
            
            %%%% cycle around the images and extract measurements ---------
            stackinfo = imfinfo(FileTif);
            maxround = length(stackinfo);
            for i3 = 1:maxround
                % load image
                FluorImage = imread(FileTif,'Index',i3);
                FluorImage_BK = FluorImage - imopen(imclose(FluorImage,strel('disk',ceil(options.cellsize/10))),strel('disk',ceil(options.cellsize*3)));
                % extract pixel values from fluo image using the segmented images    
                Nuclei_stats = regionprops(lb_NucMask,FluorImage_BK,'PixelValues');
                Cytopl_stats = regionprops(lb_CytMask,FluorImage_BK,'PixelValues');
                
                Results.MedianNucSign(startcell:endcell,i3) = cellfun(@median,{Nuclei_stats.PixelValues});
                Results.MedianCytSign(startcell:endcell,i3) = cellfun(@median,{Cytopl_stats.PixelValues});
                Results.MeanNucSign(startcell:endcell,i3) = cellfun(@mean,{Nuclei_stats.PixelValues});
                Results.MeanCytSign(startcell:endcell,i3) = cellfun(@mean,{Cytopl_stats.PixelValues});
     
            end
            Results.FieldFlag(i1,i2)=1; 
            save([filename.folders.main filename.folders.output filename.folders.resultfile],'Results', '-v7.3') %Saving matrix
        end
   end
end