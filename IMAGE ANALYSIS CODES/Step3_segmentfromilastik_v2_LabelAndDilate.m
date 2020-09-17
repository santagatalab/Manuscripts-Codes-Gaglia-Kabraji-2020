function Step3_segmentfromilastik_v2_LabelAndDilate(filename, options) 
% open the tiff file of the full stack images and segments them based on
% Ilastik Probabilities file 
% Inputs FullStacks and Ilastik Probabilities files 

tic
mkdir([filename.folders.main filename.folders.output filename.folders.ilastikseg])

for k = 1:length(filename.tissues) 
    fold = filename.tissues{k}; 
    addpath(filename.folders.main)
    coordmat = [filename.folders.main filename.folders.output filename.folders.coordinates fold '.mat'];
    load(coordmat, 'Coordinates')
    [r, c] = size(Coordinates.Field);
    for i1 = 1:r
        for i2 = 1:c
            core = [fold '_Field_' num2str(i1 , filename.dim) '_' num2str(i2 , filename.dim)];
%             FileTif = [filename.folders.main filename.folders.output  fold filesep filename.folders.fullstacks core  filename.suffix];
            FileTif = [filename.folders.main filename.folders.output  filename.folders.ilastikstacks core  filename.suffix];
            segfile = [filename.folders.main filename.folders.output filename.folders.ilastikseg core '_Seg.tif'];
            IlastikTif = [filename.folders.main filename.folders.output filename.folders.ilastikprob core  filename.ilastiksuffix];
            
            if exist(segfile,'file')==2
                disp(['Image ' core ' alredy segmented'])
                continue
            elseif exist(IlastikTif,'file')~=2
                disp(IlastikTif)
                disp(['No ilastik prob file for ' core])
                continue
            else
                
            end
            
            disp(FileTif)

            try 
                DAPI = imread(FileTif,'Index',1);
                ominfo = imfinfo(FileTif);
                DAPI_last = imread(FileTif,'Index',length(ominfo)-3);
            catch
                disp('Error: DAPI not found')
                continue
            end
            
            
            try
                IlastikRGB = imread(IlastikTif);
            catch
                disp('Error: Ilastik RGB not found')
                disp(IlastikTif)
                continue
            end

            NucProb = IlastikRGB(:,:,options.nuc);
            CytProb = IlastikRGB(:,:,options.cyt);
            BkgProb = IlastikRGB(:,:,options.backgr);
            
            NucProb = imopen(NucProb,strel('octagon',3));
            BkgProb = imopen(BkgProb,strel('octagon',3));
           

            % pre-filter easy ones
            nuc_mask = zeros(size(NucProb));
            temp_mask = zeros(size(NucProb));
            else_mask = zeros(size(NucProb));

            index = BkgProb > options.max_prob*options.bkgdprob_min | CytProb > options.max_prob*options.cytprob_min; 
            else_mask(index) = 1;

            index = NucProb > options.max_prob*options.nucprob_min; 
            temp_mask(index) = 1;

            nuc_mask = temp_mask & ~else_mask;
            nuc_mask = bwareaopen(nuc_mask,options.cellsize*2,4);

            sImage = NucProb;
            ThresImage=nuc_mask;

            seeds = imclose(sImage,strel('disk',2));
            seeds = imgaussfilt(seeds,5,'FilterSize',round(options.cellsize));
            seeds = imregionalmax(seeds);
            seeds=seeds&ThresImage;
            seeds=imdilate(seeds,strel('disk',2));

            % threshold and watershed image
            L=watershed(imimposemin(-sImage,seeds));
            Nuclei= uint32(ThresImage&L);
            
            Nuclei = bwareaopen(Nuclei,options.cellsize*2,4);
            Nuclei = imfill(Nuclei,'holes');
            forcheck = uint16(Nuclei);
            Nuclei = bwlabel(uint32(Nuclei));
            Nuclei = imdilate(Nuclei,strel('disk',2));
            
%            % define the nuclear area less stringently
%            
%            ThreshImg = NucProb > options.max_prob*0.2 & CytProb < options.max_prob*0.7 & BkgProb < options.max_prob*0.50;
%            NucVect = DAPI(ThreshImg);
%             if length(NucVect)>200
%                
%                 ThreshImg = imfill(ThreshImg,'holes');
%            
%                 % watershed image
%                 coImage=imopen(imclose(DAPI,strel('disk',round(options.cellsize/10))),strel('disk',round(options.cellsize*2)));
%                 sImage = DAPI-coImage;
%             
%                 % now define seeds very stringently
%                 seeds_img = NucProb > options.max_prob*0.5 & CytProb < options.max_prob*0.5;
%                 seeds_img = imfill(seeds_img,'holes');
%                 seeds=seeds_img&ThreshImg;
%                 
%                   
%                 L=watershed(imimposemin(-sImage,seeds));
%                 Nuclei=ThreshImg&L;
%                 Nuclei=bwareaopen(Nuclei,options.cellsize);
%                 
%             else
%                 Nuclei = zeros(size(DAPI));
%                 seeds = Nuclei;
%             end
           
           check_img = cat(3,uint16(Nuclei),uint16(DAPI),uint16(forcheck));
           check_file = [filename.folders.main filename.folders.output filename.folders.ilastikseg core  '_checkseg.tif'];
           imwrite(check_img,check_file)
           
           % label the nuclei
           
           saveastiff(uint32(Nuclei),segfile)
           
           if options.pausetime > 0
               disp(['Pausing for ' num2str(options.pausetime) ' seconds'])
               pause(options.pausetime)
           end
        end
    end
end 
