function [ParamModel, Data, InputData, Wins]=info_callsemantic(PSTH,JackKnifeTrials,VocType, ParamModel,  Calfilename)
FIG=0;

%% Deals with input parameters
if nargin<4
    ParamModel = struct();
end
if  ~isfield(ParamModel,'MinWin') || isempty(ParamModel.MinWin)
    ParamModel.MinWin = 10; % end point of the first analysis window (spectrogram and neural response)
end
if ~isfield(ParamModel,'MaxWin') || isempty(ParamModel.MaxWin)
    ParamModel.MaxWin = 1000; %end point of the last anaysis window for...
    ... neural response and end point of the largest analysis window for...
        ... spectrogram
end

if ~isfield(ParamModel,'Increment') || isempty(ParamModel.Increment)
    ParamModel.Increment = 10; %increase the size of the spectro window with a Xms pace
end
if ~isfield(ParamModel,'NeuroBin') || isempty(ParamModel.NeuroBin)
    ParamModel.NeuroBin = 10; % size of the window (ms) within which the neural response is analyzed
                               % The end of the window of analysis is
                               % determined by the Increment and ResDelay (see below).
end
if ~isfield(ParamModel,'ResDelay') || isempty(ParamModel.ResDelay)
    ParamModel.ResDelay = 0; % Delay in ms between the end of the...
    ... spectrogram window and the end of the neural response window
end

% Number of bootstraps
if ~isfield(ParamModel, 'NbBoot_Info') || isempty(ParamModel.NbBoot_Info)
    ParamModel.NbBoot_Info = 10;
end

% now use the predetermined sets of JK indices
Avail_NJK = length(ParamModel.SetIndices_JK);


if ParamModel.NbBoot_Info > Avail_NJK
    fprintf('WARNING: Only %d possible calculations of sets of JK points while you are asking for %d in the calculation of information\n', Avail_NJK, ParamModel.NbBoot_Info);
    ParamModel.NbBoot_Info = Avail_NJK;
else
    fprintf('WARNING: %d possible calculations of sets of JK points. you are asking for %d in the calculation of information\n', Avail_NJK, ParamModel.NbBoot_Info);
end

%% Set up parameters of the function
% define the list of end points of time bins at which information is
% calculated
Wins = ParamModel.MinWin:ParamModel.Increment:ParamModel.MaxWin;

% # of bins in the data
WinNum = length(Wins);

% Number of stims in the data set
NbStims = length(PSTH);

% Random grouping of vocalizations along categories
VocTypeRand = VocType(randperm(length(VocType)));

%% Configure Parallel computing
if ~isempty(strfind(getenv('HOSTNAME'),'.savio')) || ~isempty(strfind(getenv('HOSTNAME'),'.brc'))
    MyParPool = parpool(min(ParamModel.NbBoot_Info,str2num(getenv('SLURM_CPUS_ON_NODE'))),'IdleTimeout', Inf);
    system('mkdir -p /global/scratch/$USER/$SLURM_JOB_ID')
    [~,JobID] = system('echo $SLURM_JOB_ID');
    parcluster.JobStorageLocation = ['/global/scratch/jelie/' JobID];    
end

%% Initialize output variables for the calculation of instantaneous information
Rate4InfoStim = nan(NbStims,WinNum);
Rate4InfoStim_Boot = cell(1,ParamModel.NbBoot_Info);
CsteRate4InfoStim_Boot = cell(1,ParamModel.NbBoot_Info);
Stim_entropy = nan(1,WinNum);
Category_entropy = nan(1,WinNum);
Stim_info = nan(1,WinNum);
Category_info = nan(1,WinNum);
CategoryRand_info = nan(1,WinNum);

P_YgivenS = cell(1,WinNum);
P_YgivenC = cell(1,WinNum);
P_YgivenCRand = cell(1,WinNum);
P_YgivenS_Bootstrap = cell(1,ParamModel.NbBoot_CumInfo);
P_YgivenC_Bootstrap = cell(1,ParamModel.NbBoot_CumInfo);
P_YgivenCRand_Bootstrap = cell(1,ParamModel.NbBoot_CumInfo);
P_YgivenS_Bootstrap_csteRate = cell(1,ParamModel.NbBoot_CumInfo);
P_YgivenC_Bootstrap_csteRate = cell(1,ParamModel.NbBoot_CumInfo);
stim_info_JKBoot_infT_mean = nan(ParamModel.NbBoot_Info,WinNum);
category_info_JKBoot_infT_mean = nan(ParamModel.NbBoot_Info,WinNum);
categoryRand_info_JKBoot_infT_mean = nan(ParamModel.NbBoot_Info,WinNum);

stim_info_JKBoot_infT = cell(ParamModel.NbBoot_Info,1);
category_info_JKBoot_infT = cell(ParamModel.NbBoot_Info,1);
categoryRand_info_JKBoot_infT = cell(ParamModel.NbBoot_Info,1);
Stim_info_JKBoot = cell(ParamModel.NbBoot_Info,1);
Category_info_JKBoot = cell(ParamModel.NbBoot_Info,1);
CategoryRand_info_JKBoot = cell(ParamModel.NbBoot_Info,1);
stim_info_JKBoot_infT_var = nan(ParamModel.NbBoot_Info,WinNum);
category_info_JKBoot_infT_var = nan(ParamModel.NbBoot_Info,WinNum);
categoryRand_info_JKBoot_infT_var = nan(ParamModel.NbBoot_Info,WinNum);
stim_info_JKBoot_infT_csteRate_var = nan(ParamModel.NbBoot_Info,1);
category_info_JKBoot_infT_csteRate_var = nan(ParamModel.NbBoot_Info,1);
stim_info_JKBoot_infT_csteRate_mean = nan(ParamModel.NbBoot_Info,1);
category_info_JKBoot_infT_csteRate_mean=nan(ParamModel.NbBoot_Info,1);
        


%% Now loop through bins: bin spike patterns, calculate instantaneous information and output the distribution probabilities of neural responses given the stimulus or category
%parfor
parfor ww = 1:WinNum
    Tstart=tic;
    fprintf(1,'Instantaneous info exact spike patterns %d/%d window\n', ww, WinNum);
    Win = Wins(ww);
    FirstTimePoint = (Win - ParamModel.NeuroBin+ ParamModel.ResDelay)*ParamModel.Response_samprate/1000 +1;
    LastTimePoint = (Win + ParamModel.ResDelay)*ParamModel.Response_samprate/1000;
     
    % Calculating info about the stims and the categories for actual values of spike rates
    [Local_Output] = info_model_Calculus_wrapper2(PSTH, FirstTimePoint, LastTimePoint,VocType,VocTypeRand,ParamModel.Response_samprate);
    Rate4InfoStim(:,ww) = Local_Output.InputdataStim;
    Stim_info(ww) = Local_Output.stim_value;
    Category_info(ww) = Local_Output.cat_value;
    CategoryRand_info(ww) = Local_Output.catRand_value;
    P_YgivenS{ww} = Local_Output.P_YgivenS;
    P_YgivenC{ww} = Local_Output.P_YgivenC;
    P_YgivenCRand{ww} = Local_Output.P_YgivenCRand;
    Stim_entropy(ww) = Local_Output.stim_entropy;
    Category_entropy(ww) = Local_Output.cat_entropy; 
   fprintf('Instantaneous Info: Done bin %d/%d after %f sec\n', ww, WinNum, toc(Tstart));
end

% Calculate the info if the neural coding was using an average rate for
% each stimulus
Tstart2=tic;
if size(PSTH,1)==NbStims
    PSTH_Local = cell2mat(PSTH);
else
    PSTH_Local = cell2mat(PSTH');
end
CsteRate = mean(PSTH_Local(:,1:ParamModel.MaxWin),2);
PSTH_cste = repmat(CsteRate,1,ParamModel.NeuroBin);
% Calculating info about the stims and the categories for time averaged value of spike rates
[Local_Output] = info_model_Calculus_wrapper2(PSTH_cste, 1, ParamModel.NeuroBin,VocType,[],ParamModel.Response_samprate);
CsteRate4InfoStim = Local_Output.InputdataStim;
Stim_info_csteRate = Local_Output.stim_value;
Category_info_csteRate = Local_Output.cat_value;
P_YgivenS_csteRate = Local_Output.P_YgivenS;
P_YgivenC_csteRate = Local_Output.P_YgivenC;
fprintf('Instantaneous Info: Done constant Rate after %f sec\n', toc(Tstart2));

    
% Bootstrapping the calculation of information with Jackknife estimations of spike rate
%parfor
parfor bb=1:ParamModel.NbBoot_Info
    fprintf(1,'%d/%d bootstrap instantaneous info with Jackknife estimates of spike rates\n', bb, ParamModel.NbBoot_Info);
    
    % Initialize some variables for that bootstrap
    NJKsets = size(ParamModel.SetIndices_JK{bb},1);
    Rate4InfoStim_Boot{bb} = cell(NJKsets,1);
    Stim_info_JKSets = nan(NJKsets,WinNum);
    Category_info_JKSets = nan(NJKsets,WinNum);
    CategoryRand_info_JKSets = nan(NJKsets,WinNum);
    P_YgivenS_JKSets = cell(NJKsets, WinNum);
    P_YgivenC_JKSets = cell(NJKsets, WinNum);
    P_YgivenCRand_JKSets = cell(NJKsets, WinNum);
    CsteRate4InfoStim_Boot{bb} = cell(NJKsets,1);
    Stim_info_csteRate_JKSets = nan(NJKsets,1);
    Category_info_csteRate_JKSets = nan(NJKsets,1);
    P_YgivenS_csteRate_JKSets = cell(NJKsets, 1);
    P_YgivenC_csteRate_JKSets = cell(NJKsets, 1);

    % Loop through JK sets for that bootstrap
    for jk = 1:NJKsets
        % systematically choose without replacement a different set of JK trials for the stims for each bootstrap
        JackKnifePSTH = cell(NbStims,1);
        for st = 1:NbStims
            PSTH_Local = JackKnifeTrials{st};
            JackKnifePSTH{st} = PSTH_Local(ParamModel.SetIndices_JK{bb}(jk,st),:);
        end

        % Then run the calculation of information on all windows
        Rate4InfoStim_Boot{bb}{jk} = nan(NbStims,WinNum);
        for ww_in = 1:WinNum
            fprintf(1,'JK instantaneous info Boostrap %d/%d set %d/%d window %d/%d\n',bb,ParamModel.NbBoot_Info, jk,NJKsets, ww_in, WinNum);
            Win = Wins(ww_in);
            FirstTimePoint = (Win - ParamModel.NeuroBin+ ParamModel.ResDelay)*ParamModel.Response_samprate/1000 +1;
            LastTimePoint = (Win + ParamModel.ResDelay)*ParamModel.Response_samprate/1000;        
        
            [Local_Output] = info_model_Calculus_wrapper2(JackKnifePSTH, FirstTimePoint, LastTimePoint, VocType, VocTypeRand,ParamModel.Response_samprate);
        
            Rate4InfoStim_Boot{bb}{jk}(:,ww_in) = Local_Output.InputdataStim;
            Stim_info_JKSets(jk,ww_in) = Local_Output.stim_value;
            Category_info_JKSets(jk,ww_in) = Local_Output.cat_value;
            CategoryRand_info_JKSets(jk,ww_in) = Local_Output.catRand_value;
            if bb <= ParamModel.NbBoot_CumInfo
                P_YgivenS_JKSets{jk,ww_in} = Local_Output.P_YgivenS;
                P_YgivenC_JKSets{jk,ww_in} = Local_Output.P_YgivenC;
                P_YgivenCRand_JKSets{jk,ww_in} = Local_Output.P_YgivenCRand;
            end
        end
        
        % Calculate the info if the neural coding was using an average rate for
        % each stimulus
        PSTH_Local = cell2mat(JackKnifePSTH);
        CsteRate = mean(PSTH_Local(:,1:ParamModel.MaxWin),2);
        PSTH_cste = repmat(CsteRate,1,ParamModel.NeuroBin);
        % Calculating info about the stims and the categories for actual values of spike rates
        [Local_Output] = info_model_Calculus_wrapper2(PSTH_cste, 1, ParamModel.NeuroBin,VocType, [], ParamModel.Response_samprate);
        CsteRate4InfoStim_Boot{bb}{jk} = Local_Output.InputdataStim;
        Stim_info_csteRate_JKSets(jk) = Local_Output.stim_value;
        Category_info_csteRate_JKSets(jk) = Local_Output.cat_value;
        P_YgivenS_csteRate_JKSets{jk} = Local_Output.P_YgivenS;
        P_YgivenC_csteRate_JKSets{jk} = Local_Output.P_YgivenC;
        fprintf('Instantaneous Info: Done constant Rate Boostrap %d/%d set %d/%d\n', bb,ParamModel.NbBoot_Info,jk,NJKsets);
    end
    % Estimate information for infinite number of trials
    % ordonnée à l'origine: b = (y1*x2 - y2*x1)/(x2-x1) = I1 - N2 * (I2 -
    % I1) / (N1 - N2) = I1 * N1 - (N1-1)*I2 (when N2=N1-1)
    stim_info_JKBoot_infT_mean(bb,:) = ParamModel.Mean_Ntrials_perstim .* Stim_info - (ParamModel.Mean_Ntrials_perstim - 1) .* mean(Stim_info_JKSets,1);
    category_info_JKBoot_infT_mean(bb,:) = ParamModel.Mean_Ntrials_perstim .* Category_info - (ParamModel.Mean_Ntrials_perstim - 1) .* mean(Category_info_JKSets,1);
    categoryRand_info_JKBoot_infT_mean(bb,:) = ParamModel.Mean_Ntrials_perstim .* CategoryRand_info - (ParamModel.Mean_Ntrials_perstim - 1) .* mean(CategoryRand_info_JKSets,1);
    Stim_info_JKBoot{bb} = Stim_info_JKSets;
    Category_info_JKBoot{bb} = Category_info_JKSets;
    Category_info_JKBoot{bb} = CategoryRand_info_JKSets;
    

    stim_info_JKBoot_infT{bb} = ParamModel.Mean_Ntrials_perstim .* repmat(Stim_info,NJKsets,1) - (ParamModel.Mean_Ntrials_perstim - 1) .* Stim_info_JKSets;
    category_info_JKBoot_infT{bb} = ParamModel.Mean_Ntrials_perstim .* repmat(Category_info,NJKsets,1) - (ParamModel.Mean_Ntrials_perstim - 1) .* Category_info_JKSets;
    categoryRand_info_JKBoot_infT{bb} = ParamModel.Mean_Ntrials_perstim .* repmat(CategoryRand_info,NJKsets,1) - (ParamModel.Mean_Ntrials_perstim - 1) .* CategoryRand_info_JKSets;
    stim_info_JKBoot_infT_var(bb,:) = var(stim_info_JKBoot_infT{bb},0,1);
    category_info_JKBoot_infT_var(bb,:) = var(category_info_JKBoot_infT{bb}, 0,1);
    categoryRand_info_JKBoot_infT_var(bb,:) = var(categoryRand_info_JKBoot_infT{bb}, 0,1);
    
    stim_info_JKBoot_infT_csteRate = ParamModel.Mean_Ntrials_perstim .* repmat(Stim_info_csteRate,NJKsets,1) - (ParamModel.Mean_Ntrials_perstim - 1) .* Stim_info_csteRate_JKSets;
    category_info_JKBoot_infT_csteRate = ParamModel.Mean_Ntrials_perstim .* repmat(Category_info_csteRate,NJKsets,1) - (ParamModel.Mean_Ntrials_perstim - 1) .* Category_info_csteRate_JKSets;
    stim_info_JKBoot_infT_csteRate_var(bb) = var(stim_info_JKBoot_infT_csteRate,0,1);
    category_info_JKBoot_infT_csteRate_var(bb) = var(category_info_JKBoot_infT_csteRate, 0,1);
    stim_info_JKBoot_infT_csteRate_mean(bb) = mean(stim_info_JKBoot_infT_csteRate,1);
    category_info_JKBoot_infT_csteRate_mean(bb) = mean(category_info_JKBoot_infT_csteRate,1);
    
    if bb <= ParamModel.NbBoot_CumInfo
        P_YgivenS_Bootstrap{bb} = P_YgivenS_JKSets;
        P_YgivenC_Bootstrap{bb} = P_YgivenC_JKSets;
        P_YgivenCRand_Bootstrap{bb} = P_YgivenCRand_JKSets;
        P_YgivenS_Bootstrap_csteRate{bb} = P_YgivenS_csteRate_JKSets;
        P_YgivenC_Bootstrap_csteRate{bb} = P_YgivenC_csteRate_JKSets;
    end
end
stim_info_bcorr = mean(stim_info_JKBoot_infT_mean,1);
category_info_bcorr = mean(category_info_JKBoot_infT_mean,1);
categoryRand_info_bcorr = mean(categoryRand_info_JKBoot_infT_mean,1);
stim_info_err = (mean(stim_info_JKBoot_infT_var,1)).^0.5;
category_info_err = (mean(category_info_JKBoot_infT_var,1)).^0.5;
categoryRand_info_err = (mean(categoryRand_info_JKBoot_infT_var,1)).^0.5;

stim_info_bcorr_csteRate = mean(stim_info_JKBoot_infT_csteRate_mean);
category_info_bcorr_csteRate = mean(category_info_JKBoot_infT_csteRate_mean);
stim_info_err_csteRate = (mean(stim_info_JKBoot_infT_csteRate_var,1)).^0.5;
category_info_err_csteRate = (mean(category_info_JKBoot_infT_csteRate_var,1)).^0.5;



% Stuff in results in structure
InputData.VocType = VocType;
InputData.VocTypeRand = VocTypeRand;
InputData.Rate4InfoStim_Boot = Rate4InfoStim_Boot;
InputData.Rate4InfoStim = Rate4InfoStim;
InputData.CsteRate4InfoStim_Boot = CsteRate4InfoStim_Boot;
InputData.CsteRate4InfoStim = CsteRate4InfoStim;

Data.P_YgivenS_Bootstrap = P_YgivenS_Bootstrap;
Data.P_YgivenC_Bootstrap = P_YgivenC_Bootstrap;
Data.P_YgivenCRand_Bootstrap = P_YgivenCRand_Bootstrap;

Data.P_YgivenS_Bootstrap_csteRate = P_YgivenS_Bootstrap_csteRate;
Data.P_YgivenC_Bootstrap_csteRate = P_YgivenC_Bootstrap_csteRate;

Data.P_YgivenS = P_YgivenS;
Data.P_YgivenC = P_YgivenC;
Data.P_YgivenCRand = P_YgivenCRand;
Data.stim_info_bcorr = stim_info_bcorr;
Data.category_info_bcorr = category_info_bcorr;
Data.categoryRand_info_bcorr = categoryRand_info_bcorr;
Data.stim_info_err = stim_info_err;
Data.category_info_err = category_info_err;
Data.categoryRand_info_err = categoryRand_info_err;

Data.P_YgivenS_csteRate = P_YgivenS_csteRate;
Data.P_YgivenC_csteRate = P_YgivenC_csteRate;
Data.stim_info_bcorr_csteRate = stim_info_bcorr_csteRate;
Data.category_info_bcorr_csteRate = category_info_bcorr_csteRate;
Data.stim_info_err_csteRate = stim_info_err_csteRate;
Data.category_info_err_csteRate = category_info_err_csteRate;


Data.stim_entropy = Stim_entropy;
Data.category_entropy = Category_entropy;
Data.stim_info = Stim_info;
Data.stim_info_JKBoot = Stim_info_JKBoot;
Data.category_info = Category_info;
Data.category_info_JKBoot = Category_info_JKBoot;
Data.categoryRand_info = CategoryRand_info;
Data.categoryRand_info_JKBoot = CategoryRand_info_JKBoot;

Data.stim_info_bcorr_Boot = stim_info_JKBoot_infT_mean;
Data.category_info_bcorr_Boot = category_info_JKBoot_infT_mean;
Data.categoryRand_info_bcorr_Boot = categoryRand_info_JKBoot_infT_mean;

Data.stim_info_JKBoot_infT = stim_info_JKBoot_infT;
Data.category_info_JKBoot_infT = category_info_JKBoot_infT;
Data.categoryRand_info_JKBoot_infT = categoryRand_info_JKBoot_infT;
Data.stim_info_Boot_var = stim_info_JKBoot_infT_var;
Data.category_info_Boot_var = category_info_JKBoot_infT_var;
Data.categoryRand_info_Boot_var = categoryRand_info_JKBoot_infT_var;

%% Save what we have for now
if exist(Calfilename, 'file')==2
    save(Calfilename,'Data','ParamModel','Wins','InputData','-append');
else
    save(Calfilename,'Data','ParamModel','Wins','InputData');
end

 %% Plot the results if requested
 if FIG
     ColorCode = get(groot,'DefaultAxesColorOrder');
     ColorCode = [ColorCode ; 0.85 0.6940 0.556; 0.301 0.078 0.556; 0.929 0.184 0.188; 0.494 0.6740 0.933;0.466 0.745 0.184;0.635 0.694 0.7410];
     figure()
     for ss=1:NbStims
         plot(InputData.Rate4InfoStim(ss,:)./ParamModel.NeuroBin,'LineWidth',2, 'Color','g')
         hold on
         for bb=1:ParamModel.NbBoot_CumInfo
             NJK = length(InputData.Rate4InfoStim_Boot{bb});
             Local_bootrate = nan(NJK, size(InputData.Rate4InfoStim_Boot{1}{1},2));
             for jk=1:NJK
                plot(InputData.Rate4InfoStim_Boot{bb}{jk}(ss,:)./ParamModel.NeuroBin, 'Color','k')
                hold on
                if bb==1 && jk==1
                    plot(mean(Local_bootrate,1), 'LineWidth',2,'Color','r')%this is just for the legend
                    hold on
                    legend('Actual spike rate','individual bootstrap', 'Average bootstrapped spike rate')
                end
                Local_bootrate(jk,:) = InputData.Rate4InfoStim_Boot{bb}{jk}(ss,:);
             end
             plot(mean(Local_bootrate,1)./ParamModel.NeuroBin, 'LineWidth',2,'Color','r')
         end
         hold off
         Xtickposition=get(gca,'XTick');
         set(gca,'XTickLabel', Xtickposition*ParamModel.Increment)
         xlabel('Time ms')
         ylabel('Spike rate spike/ms')
         title(sprintf('Stim %d/%d',ss, NbStims))
         pause(1)
     end
     
     figure()
     plot(var(InputData.Rate4InfoStim),'LineWidth',2)
     hold on
     for bb=1:ParamModel.NbBoot_CumInfo
         NJK = length(InputData.Rate4InfoStim_Boot{bb});
         for jk=1:NJK
            plot(var(InputData.Rate4InfoStim_Boot{bb}{jk}))
            if bb==1 && jk==1
                legend('Actual', 'Bootstrapped')
            end
            hold on
         end
     end
     hold off
     Xtickposition=get(gca,'XTick');
     set(gca,'XTickLabel', Xtickposition*ParamModel.Increment)
     xlabel('Time ms')
     ylabel('Stimulus spike rate variance')
     pause()
     
     figure()
     subplot(3,1,1)
     
     % Plotting the instantaneous information data
     plot(1:WinNum,Data.stim_info,'LineWidth',2, 'Color',ColorCode(5,:))
     hold on
     plot(1:WinNum,Data.stim_info_bcorr,'LineWidth',2, 'Color',ColorCode(5,:), 'LineStyle', '--')
     hold on
     plot(Data.stim_entropy, 'LineStyle','-.','Color','r')
     legend('Information', 'Information biais corrected','Information upper-bound given dataset size', 'Location','NorthEast');
     hold on
     line([0 WinNum], [0 0], 'LineStyle','-.','Color','k')
     hold on
     ylim([-0.5 log2(NbStims)+1])
     Xtickposition=get(gca,'XTick');
     set(gca,'XTickLabel', Xtickposition*ParamModel.NeuroBin)
     xlabel('Time ms')
     ylabel('Stimulus Information in bits')
     shadedErrorBar([],Data.stim_info_bcorr, Data.stim_info_err,{'Color',ColorCode(5,:), 'LineStyle','-', 'LineWidth',1},1)
     title('Instantaneous information about stimulus')
     hold off
     
     subplot(3,1,2)
     
     % Plotting the instantaneous information data
     plot(1:WinNum,Data.category_info,'LineWidth',2, 'Color',ColorCode(7,:))
     hold on
     plot(1:WinNum,Data.category_info_bcorr,'LineWidth',2, 'Color',ColorCode(7,:), 'LineStyle', '--')
     hold on
     plot(Data.category_entropy, 'LineStyle','-.','Color','r')
     legend('Information', 'Information biais corrected','Information upper-bound given dataset size', 'Location','NorthEast');
     hold on
     line([0 WinNum], [0 0], 'LineStyle','-.','Color','k')
     hold on
     ylim([-0.5 log2(NbStims)+1])
     Xtickposition=get(gca,'XTick');
     set(gca,'XTickLabel', Xtickposition*ParamModel.NeuroBin)
     xlabel('Time ms')
     ylabel('Semantic Information in bits')
     shadedErrorBar([],Data.category_info_bcorr, Data.category_info_err,{'Color',ColorCode(7,:), 'LineStyle','-', 'LineWidth',1},1)
     title('Instantaneous information about semantic categories')
     hold off
     
     subplot(3,1,3)
     % Plotting the instantaneous information data
     plot(1:WinNum,Data.category_info*100 ./ Data.stim_info,'LineWidth',2, 'Color',ColorCode(1,:))
     hold on
     plot(1:WinNum,Data.category_info_bcorr*100 ./ Data.stim_info_bcorr,'LineWidth',2, 'Color',ColorCode(1,:), 'LineStyle', '--')
     hold on
     plot(Data.category_entropy*100 ./ Data.stim_entropy, 'LineStyle','-.','Color','r')
     legend('% Semantic Information', ' % Semantic  Information biais corrected',' % Semantic Information upper-bound given dataset size', 'Location','NorthEast');
     hold on
     line([0 WinNum], [0 0], 'LineStyle','-.','Color','k')
     hold on
     %ylim([-0.5 log2(NbStims)+1])
     Xtickposition=get(gca,'XTick');
     set(gca,'XTickLabel', Xtickposition*ParamModel.NeuroBin)
     xlabel('Time ms')
     ylabel('% Semantic Information')
     title('Percentage of Instantaneous information about semantic categories')
     hold off
     
 end
 
%% get rid of temporary files for parallel computing
if ~isempty(strfind(getenv('HOSTNAME'),'.savio')) || ~isempty(strfind(getenv('HOSTNAME'),'.brc'))
    delete(MyParPool);
    system(['rm -r ' parcluster.JobStorageLocation])
end

end