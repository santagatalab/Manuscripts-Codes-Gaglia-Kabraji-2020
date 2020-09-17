function Step1_fields_preilastik(filename, options) 

outputfold = [filename.folders.main filename.folders.output]; 
addpath(filename.folders.main)
mkdir(outputfold) 

%% Splitting into size of field(t) desired 
for i = 1:length(filename.tissues)
    t= filename.sizefield ; %minus 1 to get the field to desired size
    tic
    Coordinates = {};
    fold = filename.tissues{i};
    disp(fold) 
    
    %Making folders
    outputfolder = [outputfold ];
    addpath(outputfold)
    mkdir([outputfolder filename.folders.fullstacks])
    mkdir([outputfold filename.folders.coordinates])  
    mkdir([outputfolder filename.folders.cropfol])
    mkdir([outputfolder filename.folders.ilastikstacks])
    
    filenamesv = [outputfold filename.folders.coordinates fold];
    % this is the folder where the ometiff files are
    omefold = [filename.folders.main filename.folders.ometiff filename.tissues{i} filesep filename.ashlarsubfol];
    omefile = [omefold filename.tissues{i} '.ome.tif'];
    
    % read in the first image of the ometiff file
    DAPICycle0 = imread(omefile, 1);
    [y, x] = size(DAPICycle0); %y= # of pixels in 1 column; x = # of pixels in 1 row
    
    fld_rows = ceil(y/t);
    fld_cols = ceil(x/t);
    
    for r = 1:fld_rows
        for c = 1:fld_cols
            fieldname = [fold '_Field_' num2str(r , filename.dim) '_' num2str(c , filename.dim)];
            outfilename = [outputfolder filename.folders.fullstacks fieldname '.tif'];
            
            pix_rows = 1+t*(r-1):min(y,r*t);
            pix_cols = 1+t*(c-1):min(x,c*t);
            
            field = DAPICycle0(pix_rows, pix_cols);
            if prctile(field(:),97.5) < options.DAPIbkgd 
                kf = 0; %kf= keep field?
            else
                kf = 1;
                imwrite(field, outfilename)
            end
            Coordinates.Field{r,c} = [kf, min(pix_rows), max(pix_rows), min(pix_cols), max(pix_cols)]; %Saving coordinates
            
        end
    end
    

    toc
    save(filenamesv, 'Coordinates', 'x', 'y', 't')
    %% Loop through channels to create stacks
    for ch = 2:filename.maxround
        Cycle = imread(omefile, ch);
        
        if size(Cycle)==size(DAPICycle0)

            for r = 1:fld_rows
                for c = 1:fld_cols
                    coordvect = Coordinates.Field{r,c};
                    if coordvect(1) == 1 % ie keepfield flag is yes
                        r_coord = coordvect(2:3); %Intial, final
                        c_coord = coordvect(4:5); %Intial, final
                        field = Cycle( r_coord(1):r_coord(2), c_coord(1):c_coord(2));
                        fieldname = [fold '_Field_' num2str(r , filename.dim) '_' num2str(c , filename.dim)];

                        outfilename = [outputfolder filename.folders.fullstacks fieldname '.tif'];
                        imwrite(field, outfilename, 'WriteMode','append')  

                        % here we make the stacks for ilastik
                        if ch == filename.ilastikround
                            % we copy the stacks file to the ilastik folder
                            ilafull_filename = [outputfolder filename.folders.ilastikstacks fieldname '.tif'];
                            copyfile(outfilename,ilafull_filename);
                        end
                    end  
                end
            end
        
        fprintf([num2str(ch) ' '])
        else
            break
        end
    end
    fprintf('\n')
    disp(['Tissue ' fold ' done'])
    toc

end