%%
clear all

codedir = 'Y:\add_foleder_with_functions';
addpath(codedir)

InCell = 'experiment_folder\';
basefolder = [InCell 'experiment_folder\'];
analfolder = [basefolder 'ANALYSIS\'];
resufolder = 'result_folder\';
date = 'add_if_suffix_is_used_for_analysis';

load([ analfolder resufolder 'Results_Aggr_' date '.mat'])
load([ analfolder resufolder 'Results_Morp_' date '.mat'])
load([ analfolder resufolder 'Results_Norm_' date '.mat'])
load([ analfolder resufolder 'Results_Filt_' date '.mat'])
load([ analfolder resufolder 'Results_ROI_' date '.mat'])
load([ analfolder resufolder 'Results_Settings_' date '.mat'])

filename.basefolder = basefolder;
filename.analfolder = analfolder;
filename.resufolder = resufolder;
options.date = date;
save([ analfolder resufolder 'Results_Settings_' date '.mat'], 'filename','options','-append')


%% Define your cell types!
CellType = [];
CellType.NameList = {{'Immune','Epithelial','Stroma','Other'}};

% worth explaining what these columns mean?
CellType.Classes =  { ...
              {'Immune',           1, 1         ,{2}} ... 
             ,{'Epithelial',       1, 2         ,{18}}  ...
             ,{'Stroma',           1, 3         ,{[20 -2],[48 -2 -18]}} ... ,{'Immune Other',        100   ,{7,10,11,12,19,24}} ...
             ,{'Other',            1, 4         ,{[-2 -18 -20 -48]}} ...
             }; 

% ABlist is an ordered list of all the channels (antibodies) used for cell typing
CellType.ABlist = cellfun(@(x) x{4},CellType.Classes,'un',0);
CellType.ABlist = unique(cell2mat(cellfun(@(x) abs([x{:}]),CellType.ABlist,'un',0)));
% use this list to specify which markers are Nuclear vs Cytoplasmic
    % 1 for nuclear (default), 0 for cytoplasmic
CellType.NvC = ones(length(CellType.ABlist),'logical');
CellType.NvC = zeros(length(CellType.ABlist),'logical'); 
    % for manual override of all markers
% CellType.NvC(1,3,4) = 0; 
    % for manual override of a few markers

CellType.index = [];
CellType.codes  = [];
CellType.layer = [];
CellType.names  = {};
CellType.layerjump = 10;  % either 10 or 100 usually

for i = 1:length(CellType.Classes)
    CellType.index = [CellType.index; i];
    CellType.codes  = [CellType.codes;  CellType.Classes{i}{3}];
    CellType.layer = [CellType.layer; CellType.Classes{i}{2}];
    CellType.names  = [CellType.names;  CellType.Classes{i}{1}];
    
end

%%% Classify your cells!

totcells = length(NormResults.MedianNucMat(:,1));

%%% 20200617 addition: make a data vec out of only 
TypeData = zeros(size(NormResults.MedianNucMat),'int16');
TypeData(:,CellType.ABlist) = NormResults.MedianNucMat(:,CellType.ABlist).*CellType.NvC + ...
    NormResults.MedianCytMat(:,CellType.ABlist).*~CellType.NvC;


for layer = 1:max(CellType.layer)
    layer
    % define how many types of cells this layer has
    types = CellType.index(CellType.layer==layer);
    Matrix{layer} = zeros(totcells,length(types));
    ProbType{layer} = zeros(totcells,length(types));
    
    for type = 1:length(types)
        % count the number of criteria used
        num_criteria = size(CellType.Classes{types(type)}{4},2);
        
        ind_mat = zeros(totcells,num_criteria);
        ind_prob= zeros(totcells,num_criteria);
        for i3 = 1:num_criteria
            crit = CellType.Classes{types(type)}{4}{i3};
           
            temp_ind = ones(totcells,1);
            temp_prob= zeros(totcells,1)+NaN;
            for i4 = 1:length(crit)
                crit_ind = sign(crit(i4)) * NormResults.MedianNucMat(:,abs(crit(i4))) > 0 ;%.* uint16(Filter.all(:,abs(crit(i4)))) 
                temp_ind = temp_ind & crit_ind;
                if sign(crit(i4)) > 0
                    test = crit_ind.*double(NormResults.MedianNucMat(:,crit(i4)));
                    test(test==0) = NaN;
                    temp_prob = min([temp_prob, test],[],2);
                end 
            end
            ind_mat(:,i3)  = temp_ind;
            ind_prob(:,i3) = temp_prob;
        end
        ind_vect = max(ind_mat,[],2);
        Matrix{layer}(ind_vect>0,type) = CellType.codes(types(type)); 
        ProbType{layer}(ind_vect>0,type) = max(ind_prob(ind_vect>0,:),[],2);
    end   
end
clear temp_ind

%%%% check and resolve conflicts

% - starting from the first layer we check for conflicts and resolve each
% layer into one column instead of one matrix
% - then we check in the next matrix and first take out all the calls that
% should not have been made in the first place, because the second layer
% does not match the next one
% loop through the two-way comparison for each layer
thresh = 10;
test = [];
CleanMatrix = Matrix;
CleanProbType = ProbType;
CellType.Matrix = [];

for i1 = 1:length(Matrix)
    i1
    % from the second layer onwards check that the subcalls were made from
    % the right branch of the layer above, if not set them to zero
    if i1 > 1
        for j = 1:size(CleanMatrix{i1},2)
            wrong_branch_index = CleanMatrix{i1}(:,j) > 0 & CellType.Matrix(:,i1-1) ~= floor(CleanMatrix{i1}(:,j)/CellType.layerjump);
            CleanMatrix{i1}(wrong_branch_index,j) = 0;
            CleanProbType{i1}(wrong_branch_index,j) = 0;
        end
    end
    for j1 = 1:size(CleanMatrix{i1},2)
        for j2 = j1+1:size(CleanMatrix{i1},2)
            mat  = [CleanMatrix{i1}(:,j1) CleanMatrix{i1}(:,j2)];
            prob = [CleanProbType{i1}(:,j1) CleanProbType{i1}(:,j2)];
            
            % find conflicts
            index_conflict = mat(:,1) > 0 & mat(:,2) > 0;
            mat_conflict = mat(index_conflict,:);
            prob_conflict = prob(index_conflict,:);
            
            for doubt = 1:size(mat_conflict,1)
                [~,ind_worse] = min(prob_conflict(doubt,:));
                mat_conflict(doubt,ind_worse) = 0;
            end

            % if the difference in probability is less than say 0.1 put
            % both the celltypecodes to 0
            diff_prob = abs(prob_conflict(:,1)-prob_conflict(:,2));
            mat_conflict(diff_prob < thresh,:) = 0;

            mat(index_conflict,:) = mat_conflict;

            CleanMatrix{i1}(:,j1) = mat(:,1);
            CleanMatrix{i1}(:,j2) = mat(:,2);
           
        end
    end
    CellType.Matrix(:,i1) = max(CleanMatrix{i1},[],2);
end

% save some space
CellType.Matrix = uint16(CellType.Matrix);

% check for clean up efficiency
for i = 1:length(Matrix)
    matcheck = CleanMatrix{i};
    for j1 = 1:size(matcheck,2)
        for j2 = j1+1:size(matcheck,2)
            disp([num2str(i) ' ' num2str(j1) ' ' num2str(j2) ' ' num2str(sum(matcheck(:,j1)>0 & matcheck(:,j2)>0))])
        end
    end
end


count = 0;
for i = 1:size(CellType.Matrix,2)
    types = unique(CellType.Matrix(:,i));
    for t = 1:length(types)
        count = count + 1;
        rep(count,1) = i;
        rep(count,2) = types(t);
        rep(count,3) = sum(CellType.Matrix(:,i)==types(t));
    end
end
CellType.Count = rep;


save([filename.analfolder filename.resufolder 'Results_CellType_' options.date '.mat'],'CellType')

disp('DONE!')
%% PRINT PLOTS FOR EACH CELL TYPE

numlayers = size(CellType.Matrix,2);
ind_nospleen = ones(size(CellType.Matrix,2),1);
outplotdir= [filename.analfolder filename.resufolder 'Step4_CellType_Plots\'];
mkdir(outplotdir)

for chan = 1:options.maxround
    for layer = 1:numlayers
        
        % find the cell types of the layer and label them
        types = sort(unique(CellType.Matrix(:,layer)),'ascend');
        labels = cell(1,length(CellType.NameList{1,layer})+1);

        if types(1) == 0
            labels{1} = 'n/a';
            labels(2:length(CellType.NameList{1,layer})+1) = [CellType.NameList{1,layer}];
        else
            labels{1:length(CellType.NameList{1,layer})} = CellType.NameList{1,layer};
        end

        nuc = [];
        h_n = [];
        cyt = [];
        h_c = [];
        ind_cellnum = [];
        for type = 1:length(types)
            index = CellType.Matrix(:,layer) == types(type) & ind_nospleen & Filter.all(:,chan) == 1; 
            [nuc(:,type),h_n(:,type)] = hist(NormResults.MedianNucMat(index,chan),linspace(-4000,4000,200));
            [cyt(:,type),h_c(:,type)] = hist(NormResults.MedianCytMat(index,chan),linspace(-4000,4000,200));
            nuc(:,type) = nuc(:,type)/sum(nuc(:,type));
            cyt(:,type) = cyt(:,type)/sum(cyt(:,type));
            ind_cellnum(type) = sum(index);
            
        end
        figure(chan)
        subplot(3,numlayers,layer)
        plot(h_n,nuc)
        title(['Nuclear ' options.Markers{chan}])
        legend(labels,'Location','northeast')
        subplot(3,numlayers,layer+numlayers)
        plot(h_c,cyt)
        title(['Cytoplasmic ' options.Markers{chan}])
        legend(labels,'Location','northeast')
        subplot(3,numlayers,layer+2*numlayers)
        bar(ind_cellnum)
        set(gca,'xticklabel',labels)
        ylabel('cell number')
        xtickangle(90)
        ylim([0 max(ind_cellnum(2:end))*1.05])       
    end
    set(gcf,'Units','Normalized','OuterPosition',[0.04 0.06 0.88 0.90])
    saveas(figure(chan),[outplotdir, [num2str(chan) '_' options.Markers{chan}]],'jpeg');
end
if options.figOpt == 0
    close all
end
