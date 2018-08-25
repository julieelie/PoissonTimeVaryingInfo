function [ParamModel, Data, InputData, Wins]=info_cuminfo_callsemantic(PSTH,JackKnifeTrials,VocType, ParamModel,  Calfilename)
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

if ~isfield(ParamModel,'MaxWin_cumInfo') || isempty(ParamModel.MaxWin_cumInfo)
    ParamModel.MaxWin_cumInfo = 600; %end point of the last anaysis window for...
    ... the calculation of cumulative information
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

if ~isfield(ParamModel, 'NbBoot_CumInfo') || isempty(ParamModel.NbBoot_CumInfo)
    ParamModel.NbBoot_CumInfo = 10;
end

% Checking the number of possible jackknike given the dataset size
% NbStim = length(JackKnifeTrials);
% NJK = nan(NbStim,1);
% for st = 1:NbStim
%     PSTH_Local = JackKnifeTrials{st};
%     NJK(st) = size(PSTH_Local,1);
% end
% Avail_NJK = min(NJK);

% now use the predetermined sets of JK indices
Avail_NJK = length(ParamModel.SetIndices_JK);


if ParamModel.NbBoot_Info > Avail_NJK
    fprintf('WARNING: Only %d possible calculations of sets of JK points while you are asking for %d in the calculation of information\n', Avail_NJK, ParamModel.NbBoot_Info);
    ParamModel.NbBoot_Info = Avail_NJK;
else
    fprintf('WARNING: %d possible calculations of sets of JK points. you are asking for %d in the calculation of information\n', Avail_NJK, ParamModel.NbBoot_Info);
end
if ParamModel.NbBoot_CumInfo > Avail_NJK
    fprintf('WARNING: Only %d possible calculations of sets of JK points while you are asking for %d in the calculation of cumulative information\n', Avail_NJK, ParamModel.NbBoot_CumInfo);
    ParamModel.NbBoot_CumInfo = Avail_NJK;
else
    fprintf('WARNING: %d possible calculations of sets of JK points. you are asking for %d in the calculation of cumulative information\n', Avail_NJK, ParamModel.NbBoot_CumInfo);
end



% Set parameters for the number of samples that should be tested in the
% MonteCarlo estimation of the cumulative information with fix number of
% samples
if ~isfield(ParamModel, 'NumSamples_MC_Cum_Info')
    ParamModel.NumSamples_MC_Cum_Info = 10^6; %Set the number of samples for the Monte Carlo approximation of the cumulative information 10^7 takes too much memory prefers lower numbers
elseif isempty(ParamModel.NumSamples_MC_Cum_Info)
    fprintf(1,'No Monte Carlo approximation of the cumulative information with a fix number of samples will be calculated\n');
end

% Set parameters for the max number of samples that should be tested in the
% optimal MonteCarlo estimation of the cumulative information
if ~isfield(ParamModel, 'MaxNumSamples_MCopt_Cum_Info')
    ParamModel.MaxNumSamples_MCopt_Cum_Info = 5*10^6; %Set the number of samples for the Monte Carlo approximation of the cumulative information 10^7 takes too much memory prefers lower numbers
elseif isempty(ParamModel.MaxNumSamples_MCopt_Cum_Info)
    fprintf(1,'No Monte Carlo approximation of the cumulative information with an optimized number of samples will be calculated\n');
end

% Set the Parameters of the Markov approximation of the cumulative
% information
if ~isfield(ParamModel, 'MarkovParameters_Cum_Info')
    ParamModel.MarkovParameters_Cum_Info = [2 3 4;1 1 1];
elseif isempty(ParamModel.MarkovParameters_Cum_Info)
    fprintf(1,'No Markov approximation of the cumulative information will be calculated\n');
end

% Set the fix time history to calculate the exact cumulative information
% (above 4 will certainly bug the computer asking for too big matrices)
if ~isfield(ParamModel, 'ExactHist')
    ParamModel.ExactHist = 4;
elseif isempty(ParamModel.ExactHist)
    fprintf(1,'There will be no exact calculation of cumulative information\n');
end

if nargin<5
    saveonline = 0;
else
    saveonline = 1;
end

%% Set up parameters of the function
% define the list of end points of time bins at which information is
% calculated
Wins = ParamModel.MinWin:ParamModel.Increment:ParamModel.MaxWin;
Wins_cumInfo = ParamModel.MinWin:ParamModel.Increment:ParamModel.MaxWin_cumInfo;

% # of bins in the data
WinNum = length(Wins);
WinNum_cumInfo = length(Wins_cumInfo);

% Number of stims in the data set
NbStims = length(PSTH);

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

P_YgivenS = cell(1,WinNum);
P_YgivenC = cell(1,WinNum);
P_YgivenS_Bootstrap = cell(1,ParamModel.NbBoot_CumInfo);
P_YgivenC_Bootstrap = cell(1,ParamModel.NbBoot_CumInfo);
P_YgivenS_Bootstrap_csteRate = cell(1,ParamModel.NbBoot_CumInfo);
P_YgivenC_Bootstrap_csteRate = cell(1,ParamModel.NbBoot_CumInfo);
stim_info_JKBoot_infT_mean = nan(ParamModel.NbBoot_Info,WinNum);
category_info_JKBoot_infT_mean = nan(ParamModel.NbBoot_Info,WinNum);

stim_info_JKBoot_infT = cell(ParamModel.NbBoot_Info,1);
category_info_JKBoot_infT = cell(ParamModel.NbBoot_Info,1);
Stim_info_JKBoot = cell(ParamModel.NbBoot_Info,1);
Category_info_JKBoot = cell(ParamModel.NbBoot_Info,1);
stim_info_JKBoot_infT_var = nan(ParamModel.NbBoot_Info,WinNum);
category_info_JKBoot_infT_var = nan(ParamModel.NbBoot_Info,WinNum);
stim_info_JKBoot_infT_csteRate_var = nan(ParamModel.NbBoot_Info,1);
category_info_JKBoot_infT_csteRate_var = nan(ParamModel.NbBoot_Info,1);
stim_info_JKBoot_infT_csteRate_mean = nan(ParamModel.NbBoot_Info,1);
category_info_JKBoot_infT_csteRate_mean=nan(ParamModel.NbBoot_Info,1);
    


%% Now loop through bins: bin spike patterns, calculate instantaneous information and output the distribution probabilities of neural responses givn the stimulus or category
%parfor
parfor ww = 1:WinNum
    Tstart=tic;
    fprintf(1,'Instantaneous info exact spike patterns %d/%d window\n', ww, WinNum);
    Win = Wins(ww);
    FirstTimePoint = (Win - ParamModel.NeuroBin+ ParamModel.ResDelay)*ParamModel.Response_samprate/1000 +1;
    LastTimePoint = (Win + ParamModel.ResDelay)*ParamModel.Response_samprate/1000;
     
    % Calculating info about the stims and the categories for actual values of spike rates
    [Local_Output] = info_model_Calculus_wrapper2(PSTH, FirstTimePoint, LastTimePoint,VocType,ParamModel.Response_samprate);
    Rate4InfoStim(:,ww) = Local_Output.InputdataStim;
    Stim_info(ww) = Local_Output.stim_value;
    Category_info(ww) = Local_Output.cat_value;
    P_YgivenS{ww} = Local_Output.P_YgivenS;
    P_YgivenC{ww} = Local_Output.P_YgivenC;
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
% Calculating info about the stims and the categories for actual values of spike rates
[Local_Output] = info_model_Calculus_wrapper2(PSTH_cste, 1, ParamModel.NeuroBin,VocType,ParamModel.Response_samprate);
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
    
    % Choosing a different set of JK trials for the stims for each bootstrap
%     NbStim = length(JackKnifeTrials);
%     JackKnifePSTH = cell(1,NbStim);
%     for st = 1:NbStim
%         PSTH_Local = JackKnifeTrials{st};
%         NJK = size(PSTH_Local,1);
%         JackKnifePSTH{st} = PSTH_Local(randperm(NJK,1),:);
%     end
    
    % Initialize some variables for that bootstrap
    NJKsets = size(ParamModel.SetIndices_JK{bb},1);
    Rate4InfoStim_Boot{bb} = cell(NJKsets,1);
    Stim_info_JKSets = nan(NJKsets,WinNum);
    Category_info_JKSets = nan(NJKsets,WinNum);
    P_YgivenS_JKSets = cell(NJKsets, WinNum);
    P_YgivenC_JKSets = cell(NJKsets, WinNum);
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
        
            [Local_Output] = info_model_Calculus_wrapper2(JackKnifePSTH, FirstTimePoint, LastTimePoint, VocType, ParamModel.Response_samprate);
        
            Rate4InfoStim_Boot{bb}{jk}(:,ww_in) = Local_Output.InputdataStim;
            Stim_info_JKSets(jk,ww_in) = Local_Output.stim_value;
            Category_info_JKSets(jk,ww_in) = Local_Output.cat_value;
            if bb <= ParamModel.NbBoot_CumInfo
                P_YgivenS_JKSets{jk,ww_in} = Local_Output.P_YgivenS;
                P_YgivenC_JKSets{jk,ww_in} = Local_Output.P_YgivenC;
            end
        end
        
        % Calculate the info if the neural coding was using an average rate for
        % each stimulus
        PSTH_Local = cell2mat(JackKnifePSTH);
        CsteRate = mean(PSTH_Local(:,1:ParamModel.MaxWin),2);
        PSTH_cste = repmat(CsteRate,1,ParamModel.NeuroBin);
        % Calculating info about the stims and the categories for actual values of spike rates
        [Local_Output] = info_model_Calculus_wrapper2(PSTH_cste, 1, ParamModel.NeuroBin,VocType,ParamModel.Response_samprate);
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
    Stim_info_JKBoot{bb} = Stim_info_JKSets;
    Category_info_JKBoot{bb} = Category_info_JKSets;
    

    stim_info_JKBoot_infT{bb} = ParamModel.Mean_Ntrials_perstim .* repmat(Stim_info,NJKsets,1) - (ParamModel.Mean_Ntrials_perstim - 1) .* Stim_info_JKSets;
    category_info_JKBoot_infT{bb} = ParamModel.Mean_Ntrials_perstim .* repmat(Category_info,NJKsets,1) - (ParamModel.Mean_Ntrials_perstim - 1) .* Category_info_JKSets;
    stim_info_JKBoot_infT_var(bb,:) = var(stim_info_JKBoot_infT{bb},0,1);
    category_info_JKBoot_infT_var(bb,:) = var(category_info_JKBoot_infT{bb}, 0,1);
    
    stim_info_JKBoot_infT_csteRate = ParamModel.Mean_Ntrials_perstim .* repmat(Stim_info_csteRate,NJKsets,1) - (ParamModel.Mean_Ntrials_perstim - 1) .* Stim_info_csteRate_JKSets;
    category_info_JKBoot_infT_csteRate = ParamModel.Mean_Ntrials_perstim .* repmat(Category_info_csteRate,NJKsets,1) - (ParamModel.Mean_Ntrials_perstim - 1) .* Category_info_csteRate_JKSets;
    stim_info_JKBoot_infT_csteRate_var(bb) = var(stim_info_JKBoot_infT_csteRate,0,1);
    category_info_JKBoot_infT_csteRate_var(bb) = var(category_info_JKBoot_infT_csteRate, 0,1);
    stim_info_JKBoot_infT_csteRate_mean(bb) = mean(stim_info_JKBoot_infT_csteRate,1);
    category_info_JKBoot_infT_csteRate_mean(bb) = mean(category_info_JKBoot_infT_csteRate,1);
    
    if bb <= ParamModel.NbBoot_CumInfo
        P_YgivenS_Bootstrap{bb} = P_YgivenS_JKSets;
        P_YgivenC_Bootstrap{bb} = P_YgivenC_JKSets;
        P_YgivenS_Bootstrap_csteRate{bb} = P_YgivenS_csteRate_JKSets;
        P_YgivenC_Bootstrap_csteRate{bb} = P_YgivenC_csteRate_JKSets;
    end
end
stim_info_bcorr = mean(stim_info_JKBoot_infT_mean,1);
category_info_bcorr = mean(category_info_JKBoot_infT_mean,1);
stim_info_err = (mean(stim_info_JKBoot_infT_var,1)).^0.5;
category_info_err = (mean(category_info_JKBoot_infT_var,1)).^0.5;

stim_info_bcorr_csteRate = mean(stim_info_JKBoot_infT_csteRate_mean);
category_info_bcorr_csteRate = mean(category_info_JKBoot_infT_csteRate_mean);
stim_info_err_csteRate = (mean(stim_info_JKBoot_infT_csteRate_var,1)).^0.5;
category_info_err_csteRate = (mean(category_info_JKBoot_infT_csteRate_var,1)).^0.5;



% Stuff in results in structure
InputData.VocType = VocType;
InputData.Rate4InfoStim_Boot = Rate4InfoStim_Boot;
InputData.Rate4InfoStim = Rate4InfoStim;
InputData.CsteRate4InfoStim_Boot = CsteRate4InfoStim_Boot;
InputData.CsteRate4InfoStim = CsteRate4InfoStim;

Data.P_YgivenS_Bootstrap = P_YgivenS_Bootstrap;
Data.P_YgivenC_Bootstrap = P_YgivenC_Bootstrap;

Data.P_YgivenS_Bootstrap_csteRate = P_YgivenS_Bootstrap_csteRate;
Data.P_YgivenC_Bootstrap_csteRate = P_YgivenC_Bootstrap_csteRate;

Data.P_YgivenS = P_YgivenS;
Data.P_YgivenC = P_YgivenC;
Data.stim_info_bcorr = stim_info_bcorr;
Data.category_info_bcorr = category_info_bcorr;
Data.stim_info_err = stim_info_err;
Data.category_info_err = category_info_err;

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

Data.stim_info_bcorr_Boot = stim_info_JKBoot_infT_mean;
Data.category_info_bcorr_Boot = category_info_JKBoot_infT_mean;

Data.stim_info_JKBoot_infT = stim_info_JKBoot_infT;
Data.category_info_JKBoot_infT = category_info_JKBoot_infT;
Data.stim_info_Boot_var = stim_info_JKBoot_infT_var;
Data.category_info_Boot_var = category_info_JKBoot_infT_var;

 %% Save what we have for now
 if saveonline
     if exist(Calfilename, 'file')==2
        save(Calfilename,'Data','VocType','ParamModel','Wins','InputData','-append');
     else
         save(Calfilename,'Data','VocType','ParamModel','Wins','InputData');
     end
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
 
%% Now calculating cumulative information if requested by the presence of input parameters
if ~isempty(ParamModel.ExactHist) || ~isempty(ParamModel.MarkovParameters_Cum_Info) || ~isempty(ParamModel.NumSamples_MC_Cum_Info) || ~isempty(ParamModel.MaxNumSamples_MCopt_Cum_Info)
    fprintf(1,'Starting calculation of cumulative information\n')
    
    % Initializing output variables
    Data.cum_info_stim = struct();
    Data.cum_info_cat = struct();
    
    %% Running optimal Monte Carlo sampling cumulative information
    if ~isempty(ParamModel.MaxNumSamples_MCopt_Cum_Info)
        % Cumulative information for stimuli
        [Data.cum_info_stim.MonteCarloOpt_raw, Data.cum_info_stim.MonteCarloOpt_bcorr, Data.cum_info_stim.MonteCarloOpt_err, Data.cum_info_stim.MonteCarloOpt_Samples] = cumulative_info_poisson_model_MCJK_wrapper(Data.P_YgivenS, Data.P_YgivenS_Bootstrap, ParamModel.Mean_Ntrials_perstim, WinNum_cumInfo, ParamModel.NbBoot_CumInfo, ParamModel.MaxNumSamples_MCopt_Cum_Info);
        
        % Filling in the first value of cumulative info with information
        % value at bin 1
        Data.cum_info_stim.MonteCarloOpt_raw(1) = Data.stim_info(1);
        Data.cum_info_stim.MonteCarloOpt_bcorr(1) = Data.stim_info_bcorr(1);
        Data.cum_info_stim.MonteCarloOpt_err(1) = Data.stim_info_err(1);
        
        % Cumulative information for categories
        [Data.cum_info_cat.MonteCarloOpt_raw, Data.cum_info_cat.MonteCarloOpt_bcorr, Data.cum_info_cat.MonteCarloOpt_err, Data.cum_info_cat.MonteCarloOpt_Samples] = cumulative_info_poisson_model_MCJK_wrapper(Data.P_YgivenC, Data.P_YgivenC_Bootstrap, ParamModel.Mean_Ntrials_perstim,WinNum_cumInfo, ParamModel.NbBoot_CumInfo, ParamModel.MaxNumSamples_MCopt_Cum_Info);
        
        % Filling in the first value of cumulative info with information
        % value at bin 1
        Data.cum_info_cat.MonteCarloOpt_raw(1) = Data.category_info(1);
        Data.cum_info_cat.MonteCarloOpt_bcorr(1) = Data.category_info_bcorr(1);
        Data.cum_info_cat.MonteCarloOpt_err(1) = Data.category_info_err(1);
        
        %% Save what we have for now
        if saveonline
            if exist(Calfilename, 'file')==2
                save(Calfilename,'Data','ParamModel','-append');
            else
                save(Calfilename,'Data','ParamModel');
            end
        end
        %% Cumulative information if the neuron was using a constant rate coding
        P_YgivenS_local = cell(WinNum_cumInfo,1);
        P_YgivenS_Bootstrap_local = cell(1,ParamModel.NbBoot_CumInfo);
        P_YgivenC_local = cell(WinNum_cumInfo,1);
        P_YgivenC_Bootstrap_local = cell(1,ParamModel.NbBoot_CumInfo);
        for ww=1:WinNum_cumInfo
            P_YgivenS_local{ww} = Data.P_YgivenS_csteRate;
            P_YgivenC_local{ww} = Data.P_YgivenC_csteRate;
        end
        for bb=1:ParamModel.NbBoot_CumInfo
            NJKsets = size(ParamModel.SetIndices_JK{bb},1);
            P_YgivenS_Bootstrap_local{bb} = cell(NJKsets, WinNum_cumInfo);
            P_YgivenC_Bootstrap_local{bb} = cell(NJKsets, WinNum_cumInfo);
            for jk=1:NJKsets
                for ww=1:WinNum_cumInfo
                    P_YgivenS_Bootstrap_local{bb}{jk,ww} = Data.P_YgivenS_Bootstrap_csteRate{bb}{jk};
                    P_YgivenC_Bootstrap_local{bb}{jk,ww} = Data.P_YgivenC_Bootstrap_csteRate{bb}{jk};
                end
            end
        end
        
        % Cumulative information for stimuli
        [~, Data.cum_info_stim.MonteCarloOpt_bcorr_csteRate, Data.cum_info_stim.MonteCarloOpt_err_csteRate, ~] = cumulative_info_poisson_model_MCJK_wrapper(P_YgivenS_local, P_YgivenS_Bootstrap_local, ParamModel.Mean_Ntrials_perstim, WinNum_cumInfo, ParamModel.NbBoot_CumInfo, ParamModel.MaxNumSamples_MCopt_Cum_Info);
        
        % Filling in the first value of cumulative info with information
        % value at bin 1
        Data.cum_info_stim.MonteCarloOpt_bcorr_csteRate(1) = Data.stim_info_bcorr_csteRate(1);
        Data.cum_info_stim.MonteCarloOpt_err_csteRate(1) = Data.stim_info_err_csteRate(1);
        
        % Cumulative information for categories
        [~, Data.cum_info_cat.MonteCarloOpt_bcorr_csteRate, Data.cum_info_cat.MonteCarloOpt_err_csteRate, ~] = cumulative_info_poisson_model_MCJK_wrapper(P_YgivenC_local, P_YgivenC_Bootstrap_local, ParamModel.Mean_Ntrials_perstim,WinNum_cumInfo, ParamModel.NbBoot_CumInfo, ParamModel.MaxNumSamples_MCopt_Cum_Info);
        
        % Filling in the first value of cumulative info with information
        % value at bin 1
        Data.cum_info_cat.MonteCarloOpt_bcorr_csteRate(1) = Data.category_info_bcorr_csteRate(1);
        Data.cum_info_cat.MonteCarloOpt_err_csteRate(1) = Data.category_info_err_csteRate(1);
        
        %% Save what we have for now
        if saveonline
            if exist(Calfilename, 'file')==2
                save(Calfilename,'Data','ParamModel','-append');
            else
                save(Calfilename,'Data','ParamModel');
            end
        end
    end
    %% Running through time bins and calculate cumulative information with other methods
    if ~isempty(ParamModel.ExactHist) || ~isempty(ParamModel.MarkovParameters_Cum_Info) || ~isempty(ParamModel.NumSamples_MC_Cum_Info)
        % Initialize output variables
        % ... for Markov approximation
        if ~isempty(ParamModel.MarkovParameters_Cum_Info)
            HY_Markov_stim = struct();
            HY_Markov_cat = struct();
            for ss=1:size(ParamModel.MarkovParameters_Cum_Info,2)
                ModelType = sprintf('MarkovEst%d',ParamModel.MarkovParameters_Cum_Info(1,ss));
                Data.cum_info_stim.(sprintf('%s',ModelType))= nan(1,WinNum_cumInfo);
                Data.cum_info_stim.(sprintf('%s',ModelType))(1,1)= Data.stim_info(1);
                Data.cum_info_cat.(sprintf('%s',ModelType))= nan(1,WinNum_cumInfo);
                Data.cum_info_cat.(sprintf('%s',ModelType))(1,1)= Data.category_info(1);
                HY_Markov_stim.(sprintf('%s',ModelType)) = nan(1,WinNum_cumInfo);
                HY_Markov_cat.(sprintf('%s',ModelType)) = nan(1,WinNum_cumInfo);
            end
        end
        % ... for Monte Carlo with fix number of samples
        if ~isempty(ParamModel.NumSamples_MC_Cum_Info)
            for ss=1:length(ParamModel.NumSamples_MC_Cum_Info)
                ModelType=sprintf('MonteCarlo%d',log10(ParamModel.NumSamples_MC_Cum_Info(ss)));
                Data.cum_info_stim.(sprintf('%s',ModelType))= nan(2,WinNum_cumInfo);
                Data.cum_info_stim.(sprintf('%s',ModelType))(1,1)= Data.stim_info(1);
                Data.cum_info_stim.(sprintf('%s_Samples',ModelType)) = cell(1, WinNum_cumInfo);
                Data.cum_info_cat.(sprintf('%s',ModelType))= nan(2,WinNum_cumInfo);
                Data.cum_info_cat.(sprintf('%s',ModelType))(1,1)= Data.category_info(1);
                Data.cum_info_cat.(sprintf('%s_Samples',ModelType)) = cell(1, WinNum_cumInfo);
            end
        end
        % ... for the exact calculation with short time history
        if ~isempty(ParamModel.ExactHist)
            ModelType = sprintf('ExactMem0_%d',ParamModel.ExactHist);
            Data.cum_info_stim.(sprintf('%s',ModelType)) = nan(1,WinNum_cumInfo);
            Data.cum_info_cat.(sprintf('%s',ModelType)) = nan(1,WinNum_cumInfo);
            Data.cum_info_stim.(sprintf('%s',ModelType))(1) = Data.stim_info(1);
            Data.cum_info_cat.(sprintf('%s',ModelType))(1) = Data.category_info(1);
        end
        
        % Loop through time bins
        for ww =2:WinNum_cumInfo
            Tstart2=tic;
            
            % Monte Carlo with fix number of samples
            if ~isempty(ParamModel.NumSamples_MC_Cum_Info)
                for ss=1:length(ParamModel.NumSamples_MC_Cum_Info)
                    ModelType=sprintf('MonteCarlo%d',log10(ParamModel.NumSamples_MC_Cum_Info(ss)));
                    % Monte Carlo estimation with full memory cumulative information stimuli
                    [Data.cum_info_stim.(sprintf('%s',ModelType))(:,ww),~,~,Data.cum_info_stim.(sprintf('%s_Samples',ModelType)){ww}]=info_cumulative_model_Calculus(Data.P_YgivenS(1:ww),'Model#',1,'CalMode','MonteCarlo', 'MCParameter',ParamModel.NumSamples_MC_Cum_Info(ss));
                    % Monte Carlo estimation with full memory cumulative information
                    % categories
                    [Data.cum_info_cat.(sprintf('%s',ModelType))(:,ww),~,~,Data.cum_info_cat.(sprintf('%s_Samples',ModelType)){ww}]=info_cumulative_model_Calculus(Data.P_YgivenC(1:ww),'Model#',2,'CalMode','MonteCarlo', 'MCParameter',ParamModel.NumSamples_MC_Cum_Info(ss));
                end
            end
            
            % Exact calculation with short time history
            if ~isempty(ParamModel.ExactHist)
                ModelType = sprintf('ExactMem0_%d',ParamModel.ExactHist);
                % Exact calculation cumulative information on stims with ParamModel.ExactHist*10 ms memory
                [Data.cum_info_stim.(sprintf('%s',ModelType))(ww),~]=info_cumulative_model_Calculus(Data.P_YgivenS(1:ww),'Model#',1,'CalMode','Exact_Mem', 'Exact_history',ParamModel.ExactHist);
                % Exact calculation cumulative information on categories with ParamModel.ExactHist*10 ms memory
                [Data.cum_info_cat.(sprintf('%s',ModelType))(ww),~]=info_cumulative_model_Calculus(Data.P_YgivenC(1:ww),'Model#',2,'CalMode','Exact_Mem', 'Exact_history',ParamModel.ExactHist);
            end
            % Markov approximation
            if ~isempty(ParamModel.MarkovParameters_Cum_Info)
                for ss=1:size(ParamModel.MarkovParameters_Cum_Info,2)
                    ModelType = sprintf('MarkovEst%d',ParamModel.MarkovParameters_Cum_Info(1,ss));
                    if ww>2
                        % Markov chain estimation of cumulative information on stims
                        [Data.cum_info_stim.(sprintf('%s',ModelType))(ww), HY_Markov_stim.(sprintf('%s',ModelType))(ww), ~]=info_cumulative_model_Calculus(Data.P_YgivenS(1:ww),'Model#',1,'CalMode','MarkovChain', 'MarkovParameters',ParamModel.MarkovParameters_Cum_Info(:,ss), 'HY_old', HY_Markov_stim.(sprintf('%s',ModelType))(ww-1));
                        % Markov chain estimation of cumulative information on
                        % categories
                        [Data.cum_info_cat.(sprintf('%s',ModelType))(ww), HY_Markov_cat.(sprintf('%s',ModelType))(ww),~]=info_cumulative_model_Calculus(Data.P_YgivenC(1:ww),'Model#',1,'CalMode','MarkovChain', 'MarkovParameters',ParamModel.MarkovParameters_Cum_Info(:,ss),'HY_old', HY_Markov_cat.(sprintf('%s',ModelType))(ww-1));
                    else
                        % Markov chain estimation of cumulative information on stims
                        [Data.cum_info_stim.(sprintf('%s',ModelType))(ww), HY_Markov_stim.(sprintf('%s',ModelType))(ww), ~]=info_cumulative_model_Calculus(Data.P_YgivenS(1:ww),'Model#',1,'CalMode','MarkovChain', 'MarkovParameters',ParamModel.MarkovParameters_Cum_Info(:,ss));
                        % Markov chain estimation of cumulative information on
                        % categories
                        [Data.cum_info_cat.(sprintf('%s',ModelType))(ww), HY_Markov_cat.(sprintf('%s',ModelType))(ww),~]=info_cumulative_model_Calculus(Data.P_YgivenC(1:ww),'Model#',1,'CalMode','MarkovChain', 'MarkovParameters',ParamModel.MarkovParameters_Cum_Info(:,ss));
                    end
                end
            end
            fprintf('Done cumulative information bin %d/%d after %f sec\n', ww, WinNum_cumInfo, toc(Tstart2));
        end

        %% Save what we have for now
        if saveonline
            if exist(Calfilename, 'file')==2
                save(Calfilename,'Data','ParamModel','-append');
            else
                save(Calfilename,'Data','ParamModel');
            end
        end
    end
    %% Plot the cumulative information
    if FIG
        figure()
        subplot(3,1,1)
        Biais_MC = Data.cum_info_stim.MonteCarloOpt_raw - Data.cum_info_stim.MonteCarloOpt_bcorr;
        plot(1:WinNum,Data.stim_info_bcorr,'LineWidth',2, 'Color',ColorCode(5,:))
        hold on
        plot(1:WinNum_cumInfo, Data.cum_info_stim.MonteCarloOpt_raw, 'LineWidth',2, 'Color',ColorCode(6,:))
        hold on
        plot(1:WinNum_cumInfo,Data.cum_info_stim.MonteCarloOpt_bcorr, 'LineWidth',2, 'Color',ColorCode(1,:))
        hold on
        plot(1:WinNum_cumInfo, Biais_MC, 'LineWidth',2, 'Color', [ColorCode(6,:) 0.5])
        hold on
        line(1:WinNum_cumInfo, Data.stim_entropy(1:WinNum_cumInfo), 'LineStyle','-.','Color','r')
        legend('Information','Cumulative Information','Biais corrected Cumulative Information', 'Biais on cumulative information', 'Information upper-bound given dataset size', 'Location','NorthEast');
        hold on
        line([0 WinNum_cumInfo], [0 0], 'LineStyle','-.','Color','k')
        hold on
        ylim([-0.5 log2(NbStims)+2])
        Xtickposition=get(gca,'XTick');
        set(gca,'XTickLabel', Xtickposition*ParamModel.NeuroBin)
        xlabel('Time (ms)')
        ylabel('Information (bits)')
        shadedErrorBar([],Data.stim_info_bcorr, Data.stim_info_err,{'Color',ColorCode(5,:), 'LineStyle','-', 'LineWidth',1},1)
        hold on
        shadedErrorBar([],Data.cum_info_stim.MonteCarloOpt_raw, Data.cum_info_stim.MonteCarloOpt_err,{'Color',ColorCode(6,:), 'LineStyle','-', 'LineWidth',1},1)
        hold on
        shadedErrorBar([],Data.cum_info_stim.MonteCarloOpt_bcorr, Data.cum_info_stim.MonteCarloOpt_err,{'Color',ColorCode(1,:), 'LineStyle','-', 'LineWidth',1},1)
        hold on
        plot(1:WinNum_cumInfo, Data.stim_entropy(1:WinNum_cumInfo), 'LineWidth',1, 'LineStyle','-.','Color', [0.8 0.5 0.3])
        hold off
        title('Stimulus cumulative information')
        
        subplot(3,1,2)
        Biais_MC = Data.cum_info_cat.MonteCarloOpt_raw - Data.cum_info_cat.MonteCarloOpt_bcorr;
        plot(1:WinNum,Data.category_info_bcorr,'LineWidth',2, 'Color',ColorCode(5,:))
        hold on
        plot(1:WinNum_cumInfo, Data.cum_info_cat.MonteCarloOpt_raw, 'LineWidth',2, 'Color',ColorCode(6,:))
        hold on
        plot(1:WinNum_cumInfo,Data.cum_info_cat.MonteCarloOpt_bcorr, 'LineWidth',2, 'Color',ColorCode(1,:))
        hold on
        plot(1:WinNum_cumInfo, Biais_MC, 'LineWidth',2, 'Color', [ColorCode(6,:) 0.5])
        hold on
        line(1:WinNum_cumInfo, Data.category_entropy(1:WinNum_cumInfo), 'LineStyle','-.','Color','r')
        legend('Information','Cumulative Information','Biais corrected Cumulative Information', 'Biais on cumulative information', 'Information upper-bound given dataset size', 'Location','NorthEast');
        hold on
        line([0 WinNum_cumInfo], [0 0], 'LineStyle','-.','Color','k')
        hold on
        ylim([-0.5 log2(NbStims)+2])
        Xtickposition=get(gca,'XTick');
        set(gca,'XTickLabel', Xtickposition*ParamModel.NeuroBin)
        xlabel('Time (ms)')
        ylabel('Information (bits)')
        shadedErrorBar([],Data.category_info_bcorr, Data.category_info_err,{'Color',ColorCode(5,:), 'LineStyle','-', 'LineWidth',1},1)
        hold on
        shadedErrorBar([],Data.cum_info_cat.MonteCarloOpt_raw, Data.cum_info_cat.MonteCarloOpt_err,{'Color',ColorCode(6,:), 'LineStyle','-', 'LineWidth',1},1)
        hold on
        shadedErrorBar([],Data.cum_info_cat.MonteCarloOpt_bcorr, Data.cum_info_cat.MonteCarloOpt_err,{'Color',ColorCode(1,:), 'LineStyle','-', 'LineWidth',1},1)
        hold on
        plot(1:WinNum_cumInfo, Data.category_entropy(1:WinNum_cumInfo), 'LineWidth',1, 'LineStyle','-.','Color', [0.8 0.5 0.3])
        hold off
        title('Semantic cumulative information')
        
        subplot(3,1,3)
        plot(1:WinNum_cumInfo,Data.cum_info_cat.MonteCarloOpt_raw*100 ./ Data.cum_info_stim.MonteCarloOpt_raw,'LineWidth',2, 'Color',ColorCode(1,:))
        hold on
        plot(1:WinNum_cumInfo,Data.cum_info_cat.MonteCarloOpt_bcorr*100 ./ Data.cum_info_stim.MonteCarloOpt_bcorr,'LineWidth',2, 'Color',ColorCode(1,:), 'LineStyle', '--')
        hold on
        plot(Data.category_entropy*100 ./ Data.stim_entropy, 'LineStyle','-.','Color','r')
        legend('% Semantic Information', ' % Semantic  Information biais corrected',' % Semantic Information upper-bound given dataset size', 'Location','NorthEast');
        hold on
        line([0 WinNum], [0 0], 'LineStyle','-.','Color','k')
        hold on
        ylim([-0.5 log2(NbStims)+1])
        Xtickposition=get(gca,'XTick');
        set(gca,'XTickLabel', Xtickposition*ParamModel.NeuroBin)
        xlabel('Time ms')
        ylabel('% Semantic Cumulative Information')
        title('Percentage of cumulative information about semantic categories')
        hold off
    end
    
    
    
    %% Calculating Jackknife values of cumulative information for all
    % methods except Monte Carlo Optimal which is done properly in the previous
    % section for the optimal # of samples only
    
    if ~isempty(ParamModel.ExactHist) || ~isempty(ParamModel.MarkovParameters_Cum_Info)
        fprintf(1,'Bootstraping calculation of cumulative information\n');
    
        if ~isempty(ParamModel.ExactHist)
            Cum_info_stim_ExactMem0_4_Bootstrap = nan(ParamModel.NbBoot_CumInfo,WinNum_cumInfo);
            Cum_info_cat_ExactMem0_4_Bootstrap = nan(ParamModel.NbBoot_CumInfo,WinNum_cumInfo);
        end

        if ~isempty(ParamModel.MarkovParameters_Cum_Info)
            Cum_info_stim_LastMarkov_Bootstrap = nan(ParamModel.NbBoot_CumInfo,WinNum_cumInfo);
            Cum_info_cat_LastMarkov_Bootstrap = nan(ParamModel.NbBoot_CumInfo,WinNum_cumInfo);
        end
    
    
        %parfor
        parfor bb=1:ParamModel.NbBoot_CumInfo
            Tstart3=tic;
            fprintf('Bootstrap CumInfo %d/%d\n', bb, ParamModel.NbBoot_CumInfo);
            if ~isempty(ParamModel.MarkovParameters_Cum_Info)
                HY_Markov4_stim_bbstim = nan(1,WinNum_cumInfo);
                HY_Markov4_cat_bbstim = nan(1,WinNum_cumInfo);
            end
            for ww=2:WinNum_cumInfo
                fprintf('Bootstrap CumInfo %d/%d Time point %d/%d\n',bb, ParamModel.NbBoot_CumInfo, ww, WinNum_cumInfo);

                % First for the cumulative information about stimuli
                P_YgivenS_local = P_YgivenS_Bootstrap(bb,1:ww);

                if ~isempty(ParamModel.ExactHist)
                    % Exact calculation with 50 ms memory
                    [Cum_info_stim_ExactMem0_4_Bootstrap(bb,ww),~]=info_cumulative_model_Calculus(P_YgivenS_local,'Model#',bb,'CalMode','Exact_Mem', 'Exact_history',4);
                end

                if ~isempty(ParamModel.MarkovParameters_Cum_Info)
                    if ww==2
                        % Markov chain estimation 50 ms
                        [Cum_info_stim_LastMarkov_Bootstrap(bb,ww),HY_Markov4_stim_bbstim(ww),~]=info_cumulative_model_Calculus(P_YgivenS_local,'CalMode','MarkovChain', 'MarkovParameters',ParamModel.MarkovParameters_Cum_Info(:,end));
                    else
                        % Markov chain estimation 50 ms
                        [Cum_info_stim_LastMarkov_Bootstrap(bb,ww),HY_Markov4_stim_bbstim(ww),~]=info_cumulative_model_Calculus(P_YgivenS_local,'CalMode','MarkovChain', 'MarkovParameters',ParamModel.MarkovParameters_Cum_Info(:,end), 'HY_old', HY_Markov4_stim_bbstim(ww-1));
                    end
                end

                % Then same thing for the cumulative information about categories
                P_YgivenC_local = P_YgivenC_Bootstrap(bb,1:ww);

                if ~isempty(ParamModel.ExactHist)
                    % Exact calculation with 40 ms memory
                    [Cum_info_cat_ExactMem0_4_Bootstrap(bb,ww),~]=info_cumulative_model_Calculus(P_YgivenC_local,'Model#',bb,'CalMode','Exact_Mem', 'Exact_history',4);
                end

                if ~isempty(ParamModel.MarkovParameters_Cum_Info)
                    if ww==2
                        % Markov chain estimation 50 ms
                        [Cum_info_cat_LastMarkov_Bootstrap(bb,ww),HY_Markov4_cat_bbstim(ww),~]=info_cumulative_model_Calculus(P_YgivenC_local,'CalMode','MarkovChain', 'MarkovParameters',ParamModel.MarkovParameters_Cum_Info(:,end));
                    else
                        % Markov chain estimation 50 ms
                        [Cum_info_cat_LastMarkov_Bootstrap(bb,ww),HY_Markov4_cat_bbstim(ww),~]=info_cumulative_model_Calculus(P_YgivenC_local,'CalMode','MarkovChain', 'MarkovParameters',ParamModel.MarkovParameters_Cum_Info(:,end), 'HY_old', HY_Markov4_cat_bbstim(ww-1));
                    end
                end
            end
            fprintf('Done bootstrap %d on cumulative information after %f sec\n', bb, toc(Tstart3));
        end
    
        if ~isempty(ParamModel.ExactHist)
            Data.cum_info_stim.ExactMem0_4_BootJK = Cum_info_stim_ExactMem0_4_Bootstrap;
            Data.cum_info_cat.ExactMem0_4_BootJK = Cum_info_cat_ExactMem0_4_Bootstrap;
        end
    
        if ~isempty(ParamModel.MarkovParameters_Cum_Info)
            ModelTypeMarJK = sprintf('MarkovEst%d_BootJK',ParamModel.MarkovParameters_Cum_Info(1,end));
            Data.cum_info_stim.(sprintf('%s',ModelTypeMarJK)) = Cum_info_stim_LastMarkov_Bootstrap;
            Data.cum_info_cat.(sprintf('%s',ModelTypeMarJK)) = Cum_info_cat_LastMarkov_Bootstrap;

        end
    
        %% Save what we have for now
        if saveonline
            if exist(Calfilename, 'file')==2
                save(Calfilename,'Data','ParamModel','-append');
            else
                save(Calfilename,'Data','ParamModel');
            end
        end
    end
end

%% get rid of temporary files for parallel computing
if ~isempty(strfind(getenv('HOSTNAME'),'.savio')) || ~isempty(strfind(getenv('HOSTNAME'),'.brc'))
    delete(MyParPool);
    system(['rm -r ' parcluster.JobStorageLocation])
end

end