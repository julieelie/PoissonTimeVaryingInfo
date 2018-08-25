function [Info]=info_model_Calculus_wrapper(Trials, FirstTimePoint, LastTimePoint, Bootstrapsum, CatList)

if nargin<4
    Bootstrapsum = 0; %0 to take all the samples equally into account,
    % set to 1 to randomly choose (resampling) the same number of samples for each grouping factor within each grouping factor,
    % set to 2 to randomly choose (resampling) the same number of samples for each
    % grouping factor accross all the input data
end

% Change the grouping parameter to a vector if not a vector
if nargin<5
    NbStim = length(Trials);
    CatList_in = 1:NbStim;
    NbCat = NbStim;
elseif iscell(CatList) && ischar(CatList{1})
    NbStim = length(Trials);
    CatList_in = nan(1,NbStim);
    % Number of categories of stims
    IdCats = unique(CatList);
    NbCat = length(IdCats);
    for cc =1:NbCat
        CatList_in(find(strcmp(CatList, IdCats(cc))))=cc;
    end
elseif iscell(CatList) && isnumeric(CatList{1})
    CatList_in = cell2mat(CatList);
     % Number of categories of stims
    IdCats = unique(CatList_in);
    NbCat = length(IdCats);
elseif isnumeric(CatList)
    CatList_in = CatList;
     % Number of categories of stims
    IdCats = unique(CatList_in);
    NbCat = length(IdCats);
end

% Set the maximum value of Y investigated for the calculation of
% information
MaxY = 2*(LastTimePoint - FirstTimePoint +1);

% Format the input
Info.Inputdata = nan(1,NbCat);
NbEltperCat = nan(1,NbCat);
Info.Inputdata_samples = cell(1,NbCat);
    
for cc = 1:NbCat
    CatIndices = find(CatList_in == cc);
    NbEltperCat(cc) = length(CatIndices);
    Local_d = cell2mat(Trials(CatIndices));
    if ~Bootstrapsum
        Info.Inputdata(cc) = mean(sum(Local_d(:,FirstTimePoint:LastTimePoint),2));
        Info.Inputdata_samples{cc} = sum(Local_d(:,FirstTimePoint:LastTimePoint),2);
    elseif Bootstrapsum == 1
        NSamples = size(Local_d,1);
        Samples = nan(1,NSamples);
        for sp=1:NSamples
            Samples(sp) = randperm(NSamples,1);
        end
        Info.Inputdata(cc) = mean(sum(Local_d(Samples,FirstTimePoint:LastTimePoint),2));
        Info.Inputdata_samples{cc} = sum(Local_d(Samples,FirstTimePoint:LastTimePoint),2);
    elseif Bootstrapsum == 2
        NSamples = size(Local_d,1);
        Samples = nan(1,NSamples);
        Local_d = cell2mat(Trials);
        for sp=1:NSamples
            Samples(sp) = randperm(size(Local_d,1),1);
        end
        Info.Inputdata(cc) = mean(sum(Local_d(Samples,FirstTimePoint:LastTimePoint),2));
        Info.Inputdata_samples{cc} = sum(Local_d(Samples,FirstTimePoint:LastTimePoint),2);
    end
end

% Entropy of the dataset
PEltperCat = NbEltperCat./sum(NbEltperCat);
Info.data_entropy = -sum(PEltperCat.*log2(PEltperCat));
    
% Calculate information    
[Info.value,Info.P_YgivenS,~] = info_model_Calculus(Info.Inputdata, MaxY);
end