function Step2b_filtercrops_v2(filename, options) 

crop_folder = [filename.folders.main filename.folders.output filename.folders.cropfol];
list_crops = dir(crop_folder);


for i = 3:length(list_crops)
    cropname = [crop_folder list_crops(i).name];
    if strcmp(cropname(end-2:end),'tif') || strcmp(cropname(end-3:end),'tiff')
        DAPIcheck = imread(cropname,'Index',1);
        if prctile(DAPIcheck(:),options.crops.prctile) < options.crops.DAPIbkgd_crops
            delete(cropname)
        end        
    end
end