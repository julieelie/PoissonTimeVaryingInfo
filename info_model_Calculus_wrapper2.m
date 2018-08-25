function [Info]=info_model_Calculus_wrapper2(PSTH, FirstTimePoint, LastTimePoint, CatList, CatListRand, Response_samprate)


NbStim = size(PSTH,1);
% Change the grouping parameter to a vector if not a vector
if nargin<4
    fprintf(1,'No input for categories\nNo information about category is calculated\n');
    CatList_in = [];
    NbCat = 1;
elseif iscell(CatList) && ischar(CatList{1})
    CatList_in = nan(1,NbStim);
    % Number of categories of stims
    IdCats = unique(CatList);
    NbCat = length(IdCats);
    for ss =1:NbCat
        CatList_in(find(strcmp(CatList, IdCats(ss))))=ss;
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

if nargin<5
    fprintf(1,'No input for random categories\nNo information about random category is calculated\n');
    CatListRand_in = [];
    NbCatRand = 1;
elseif isempty(CatListRand)
    fprintf(1,'No input for random categories\nNo information about random category is calculated\n');
    CatListRand_in = [];
    NbCatRand = 1;
elseif iscell(CatListRand) && ischar(CatListRand{1})
    CatListRand_in = nan(1,NbStim);
    % Number of categories of stims
    IdCatsRand = unique(CatListRand);
    NbCatRand = length(IdCatsRand);
    for ss =1:NbCat
        CatListRand_in(find(strcmp(CatListRand, IdCatsRand(ss))))=ss;
    end
elseif iscell(CatListRand) && isnumeric(CatListRand{1})
    CatListRand_in = cell2mat(CatListRand);
     % Number of categories of stims
    IdCatsRand = unique(CatListRand_in);
    NbCatRand = length(IdCatsRand);
elseif isnumeric(CatListRand)
    CatListRand_in = CatListRand;
     % Number of categories of stims
    IdCatsRand = unique(CatListRand_in);
    NbCatRand = length(IdCatsRand);
end

% Set the maximum value of Y investigated for the calculation of
% information
MaxY = 2*(LastTimePoint - FirstTimePoint +1)*1000/Response_samprate; % response sampling  rate should be in hertz


% Format the input
Info.InputdataStim = nan(1,NbStim);
if iscell(PSTH)
    PSTH_Local = cell2mat(PSTH);
else
    PSTH_Local = PSTH;
end

Info.InputdataStim = sum(PSTH_Local(:,FirstTimePoint:LastTimePoint),2);

% Entropy of the stimulus dataset
Info.stim_entropy = log2(NbStim);
    
% Calculate information about stimuli   
[Info.stim_value,Info.P_YgivenS,~,~] = info_poisson_model_calculus(Info.InputdataStim, 'Ymax',MaxY);

if ~isempty(CatList_in)
    % Derive P_YgivenC and calculate information about categories
    Info.P_YgivenC = nan(MaxY+1,NbCat);
    P_YgivenS_scaled = Info.P_YgivenS ./ repmat(sum(Info.P_YgivenS,1),size(Info.P_YgivenS,1),1); % Make sure that all distributions of probabilities sum to 1 for each stimulus
    for Cat=1:NbCat
        Info.P_YgivenC(:,Cat) = mean(P_YgivenS_scaled(:,find(CatList_in == Cat)),2);
    end
    % Conditional entropy
    H_ycat = sum(sum(-Info.P_YgivenC.*log2(Info.P_YgivenC+(Info.P_YgivenC==0))))/NbCat;
    % entropy of the response
    P_y_i = mean(Info.P_YgivenC,2);
    % Calculate log probability
    H_y = sum(-P_y_i.*log2(P_y_i + (P_y_i==0)));% entropy of the response

    %% Calculate the mutual information about categories
    Info.cat_value = (H_y - H_ycat);
    % Entropy of the categories in the dataset
    Info.cat_entropy = log2(NbCat);
end

if ~isempty(CatListRand_in)
    % Derive P_YgivenCRand and calculate information about categories
    Info.P_YgivenCRand = nan(MaxY+1,NbCatRand);
    P_YgivenS_scaled = Info.P_YgivenS ./ repmat(sum(Info.P_YgivenS,1),size(Info.P_YgivenS,1),1); % Make sure that all distributions of probabilities sum to 1 for each stimulus
    for Cat=1:NbCatRand
        Info.P_YgivenCRand(:,Cat) = mean(P_YgivenS_scaled(:,find(CatListRand_in == Cat)),2);
    end
    % Conditional entropy
    H_ycat = sum(sum(-Info.P_YgivenCRand.*log2(Info.P_YgivenCRand+(Info.P_YgivenCRand==0))))/NbCatRand;
    % entropy of the response
    P_y_i = mean(Info.P_YgivenCRand,2);
    % Calculate log probability
    H_y = sum(-P_y_i.*log2(P_y_i + (P_y_i==0)));% entropy of the response

    %% Calculate the mutual information about categories
    Info.catRand_value = (H_y - H_ycat);
    % Entropy of the categories in the dataset
    Info.catRand_entropy = log2(NbCatRand);
end
end