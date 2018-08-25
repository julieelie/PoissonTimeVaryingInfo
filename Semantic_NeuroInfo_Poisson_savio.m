function [PG_Index,FanoFactor_Index, Wins, Ymean, Yvar] = Semantic_NeuroInfo_Poisson_savio(MatfilePath, SWITCH, ParamModel,Cellname)
%[Calfilename_local] = Semantic_NeuroInfo_Poisson_savio(MatfilePath, SWITCH, ParamModel,Cellname)
%[Calfilename_local] = Semantic_NeuroInfo_Poisson_savio(MatfilePath,CIType, SWITCH, ParamModel,Cellname)
% [OptimalFreqCutOff] = Semantic_NeuroInfo_Poisson_savio(MatfilePath, SWITCH, ParamModel,Cellname)
% [PG_Index,FanoFactor_Index, Wins] = Semantic_NeuroInfo_Poisson_savio(MatfilePath, SWITCH, ParamModel,Cellname)
min_nargin = 2;
%% Get the environment to figure out on which machine/cluster we are
fprintf(1,'The environment is: %s\n',getenv('HOSTNAME'))


if contains(getenv('HOSTNAME'),'ln') || contains(getenv('HOSTNAME'),'.savio') || contains(getenv('HOSTNAME'),'.brc')%savio Cluster
    Savio=1;
    fprintf(1, 'We are on savio!\n')
    addpath(genpath('/global/home/users/jelie/CODE/SingleUnitModels'));
    addpath(genpath('/global/home/users/jelie/CODE/GeneralCode'));
    addpath(genpath('/global/home/users/jelie/CODE/tlab/src/slurmbot/matlab'));
elseif ismac()
    Savio=0;
    Me = 1;
    fprintf(1, 'We are on my Mac Book Pro!\n')
    addpath(genpath('/Users/elie/Documents/CODE/SingleUnitModels'));
    addpath(genpath('/Users/elie/Documents/CODE/GeneralCode'));
    addpath(genpath('/Users/elie/Documents/CODE/STRFLab/trunk'));
else %we are on strfinator or a cluster machine
    Savio = 0;
    Me = 0;
    fprintf(1, 'Hum.. We must be on one of the lab machine\n')
    addpath(genpath('/auto/fhome/julie/Code/SingleUnitModels'));
    addpath(genpath('/auto/fhome/julie/Code/GeneralCode'));
    addpath(genpath('/auto/fhome/julie/Code/strflab'));
end
%% Start a timer for the function
TimerVal=tic;

if nargin<1
    MatfilePath = '/auto/tdrive/julie/k6/julie/matfile/FirstVoc1s_Site3_L1250R1650_e13_s0_ss1.mat';
end

%% Deal with input parameters
if nargin<(min_nargin)
    SWITCH.InfoCal=1;
elseif nargin>min_nargin && exist('CIType','var')
    if strcmp(CIType, 'Info')
        SWITCH.InfoCal=1;
        SWITCH.CumInfoCal=0;
    end
    if contains(CIType, 'CI')
        SWITCH.InfoCal=0;
        SWITCH.CumInfoCal=1;
    end
end
if exist('CIType','var')
    fprintf(1,'You are running Semantic_NeuroInfo_Poisson_savio with the type of information set up as %s\n', CIType);
end

if ~exist('SWITCH', 'var')
    SWITCH = struct();
end
if ~isfield(SWITCH,'FanoFactor') || isempty(SWITCH.FanoFactor)
    SWITCH.FanoFactor=0;
end
if ~isfield(SWITCH,'BestBin') || isempty(SWITCH.BestBin)
    SWITCH.BestBin=0;
end
if ~isfield(SWITCH,'InfoCal') || isempty(SWITCH.InfoCal)
    SWITCH.InfoCal=0;%Set to 1 if you want to calculate information on spike trains and change the name of the output file so they indicate "Info"
end
if ~isfield(SWITCH,'CumInfoCal') || isempty(SWITCH.CumInfoCal)
    SWITCH.CumInfoCal=1;%Set to 1 if you want to calculate information on spike trains and change the name of the output file so they indicate "Info"
end


if nargin<(min_nargin+2)
    ParamModel = struct();
end
if ~isfield(ParamModel,'LINK') || isempty(ParamModel.LINK)
    ParamModel.LINK='log'; %'identity'
end
if ~isfield(ParamModel,'DISTR') || isempty(ParamModel.DISTR)
    ParamModel.DISTR='poisson';%'normal'
end
% if ~isfield(ParamModel,'NeuroRes') || isempty(ParamModel.NeuroRes)
%     ParamModel.NeuroRes = 'count_gaussfiltered';
% end
if  ~isfield(ParamModel,'MinWin') || isempty(ParamModel.MinWin)
    ParamModel.MinWin = 10; % end point of the first analysis window (spectrogram and neural response)
end
if ~isfield(ParamModel,'MaxWin') || isempty(ParamModel.MaxWin)
    ParamModel.MaxWin = 600; %end point of the last anaysis window for...
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
    ParamModel.NbBoot_Info = 20; %20
end

if ~isfield(ParamModel, 'NbBoot_CumInfo') || isempty(ParamModel.NbBoot_CumInfo)
    ParamModel.NbBoot_CumInfo = 20; %20
end


if nargin<(min_nargin+3)
    [path,Cellname,ext]=fileparts(MatfilePath);
end

%% Load the unit matfile
%Res=load('/Users/elie/Documents/MATLAB/data/matfile/GreBlu9508M/ZS_Site2_L1100R1450_21.mat')
%Res=load('/Users/elie/Documents/MATLAB/data/matfile/WholeVocMat/WholeVoc_Site2_L1100R1450_e21_s0_ss1.mat')
%Res=load('/Users/elie/Documents/MATLAB/data/matfile/WholeVocMat/WholeVoc_Site2_L2000R1600_e27_s1_ss1.mat')
if Savio %savio Cluster
    Dir_local='/global/scratch/jelie/MatFiles/FirstVoc1sMat/';
    [~, File, ext] = fileparts(MatfilePath);
    Res = load(fullfile(Dir_local, [File ext]));
    %Res=loadfromTdrive_savio(MatfilePath, Dir_local);
elseif Me
    Dir_local='/Users/elie/Documents/CODE/data/matfile/FirstVoc1sMat/';
    if ~exist('ext','var')
        [~,Cellname,ext]=fileparts(MatfilePath);
    end
    Res = load([Dir_local Cellname ext]);
else
    Res = load(MatfilePath);
end



%% Get ready saving files and directories
%OutputDir_final=fullfile('/auto','tdrive','julie','k6','julie','matfile','ModMatSavio');

if Savio
    OutputDir_local='/global/scratch/jelie/MatFiles/ModMatInfo';
    OutputDirEx_local='/global/home/users/jelie/JobExecutableFiles';
    OutputDir_final=fullfile('/auto','tdrive','julie','k6','julie','matfile','ModMatInfo');
elseif Me
    OutputDir_local='/users/elie/Documents/CODE/data/matfile/ModMatInfo';
    OutputDir_final=OutputDir_local;
else
    OutputDir_final=fullfile('/auto','tdrive','julie','k6','julie','matfile','ModMatInfo');
    OutputDir_local=OutputDir_final;
end
Calfilename_local=fullfile(OutputDir_local,['InfoPoissonKDEF_' Res.Site '.mat']);
calfilename_final=fullfile(OutputDir_final,['InfoPoissonKDEF_' Res.Site '.mat']);

outfilename_local=fullfile(OutputDir_local,['slurm_out*' Res.Site '*.txt']);
outfilename_final=fullfile(OutputDir_final,['slurm_out*' Res.Site '*.txt']);


if Savio
%    This code was checking the existance of the
%        file on tdrive I'm now only running on savio and updating data on
%        tdrive
%     try
%        DoneCalc=loadfromTdrive_savio(calfilename_final,
%        OutputDir_local,1); 
%         fprintf('Found some data for this unit\n')
%         PrevData = 1;
%     catch err
%         fprintf('No Previous Data available or complete, working from the first window\nThe error is %s\n',err.identifier, err.message);
%         PrevData = 0;
%     end
    if exist(Calfilename_local, 'file') == 2
        fprintf('Found some data for this unit\n')
        PrevData = 1;
    else
        fprintf('No Previous Data available New file will be created\n');
        PrevData = 0;
    end
else
    if exist(calfilename_final, 'file') == 2
        fprintf('Found some data for this unit\n')
        PrevData = 1;
    else
        fprintf('No Previous Data available New file will be created\n');
        PrevData = 0;
    end
end



%% Get the data ready
% Need to get rid of mlnoise sections and whine sections when they
% exist. I construct a vector of indices of the right sections
DataSel=zeros(1,length(Res.VocType));
nvoc=0;
voctype=Res.VocType;
for ii=1:length(Res.VocType)
    if strcmp(voctype{ii}, 'Ag')
        nvoc=nvoc+1;
        DataSel(nvoc)=ii;
    elseif strcmp(voctype{ii}, 'Be')
        nvoc=nvoc+1;
        DataSel(nvoc)=ii;
    elseif strcmp(voctype{ii}, 'DC')
        nvoc=nvoc+1;
        DataSel(nvoc)=ii;
    elseif strcmp(voctype{ii}, 'Di')
        nvoc=nvoc+1;
        DataSel(nvoc)=ii;
    elseif strcmp(voctype{ii}, 'LT')
        nvoc=nvoc+1;
        DataSel(nvoc)=ii;
    elseif strcmp(voctype{ii}, 'Ne')
        nvoc=nvoc+1;
        DataSel(nvoc)=ii;
    elseif strcmp(voctype{ii}, 'Te')
        nvoc=nvoc+1;
        DataSel(nvoc)=ii;
    elseif strcmp(voctype{ii}, 'Th')
        nvoc=nvoc+1;
        DataSel(nvoc)=ii;
    elseif strcmp(voctype{ii}, 'song')
        nvoc=nvoc+1;
        DataSel(nvoc)=ii;
        %     elseif strcmp(voctype{dd}, 'Wh')
        %         nvoc=nvoc+1;
        %         DataSel(nvoc)=dd;
    end
end
DataSel=DataSel(1:nvoc);

%% Compute coherence and or the spectrum of the KDE of PSTH to determine the optimal window size, the scale at which neural response information is maximized
if SWITCH.BestBin
    ParamModel.Response_samprate = Res.Response_samprate;
    
    % COHERENCE on RAW spike trains
    %Data processing group spike trains in two sets containing all stims
    %but half of the trials each
    [HalfTrain1, HalfTrain2, ResponseDuration]=organiz_SpikeArray4coherence(Res.Spike_array(DataSel),ParamModel);
    
    % Compute coherence
    [CoherenceStruct]=compute_coherence_mean(HalfTrain1, HalfTrain2,Res.Response_samprate);
    CoherenceStruct.ResponseDuration = ResponseDuration;
    OptimalFreqCutOff.CoherenceStruct = CoherenceStruct;
    
    
    % Power spectrum density of KDE of spike rates
    % First retrieve the KDE of the spike patterns
    SignalTot_local = nan(nvoc,ParamModel.MaxWin);
    for vv=1:nvoc
        SignalTot_local(vv,:) = Res.PSTH_KDE_Filtered{DataSel(vv)}(1:ParamModel.MaxWin) .* hann(ParamModel.MaxWin)';
    end
    
    % Calculate the cumulative power spectrum of the signal
    Window = 0.2*Res.Response_samprate;
    SignalTot_1dim=reshape(SignalTot_local',[size(SignalTot_local,1)*size(SignalTot_local,2),1]);
    [Pxx,F] = pwelch(SignalTot_1dim, Window, [],[],Res.Response_samprate);
    if sum(isnan(Pxx))
        fprintf('Power spectrum cannot be calculated for this cell');
        OptimalFreqCutOff.PowerSpectrum.Pxx_Perc = [];
        OptimalFreqCutOff.PowerSpectrum.F = [];
        OptimalFreqCutOff.PowerSpectrum.Thresh = [];
        OptimalFreqCutOff.PowerSpectrum.LowerBandPassOpt = [];
        OptimalFreqCutOff.PowerSpectrum.LowerBandPassOpt_as4slope = [];
    else
        OptimalFreqCutOff.PowerSpectrum.Pxx_Perc = 100*cumsum(Pxx / sum(Pxx));
        OptimalFreqCutOff.PowerSpectrum.F = F;
        
        % Identify the frequency cut-off corresponding to the list of
        % thresholds
        OptimalFreqCutOff.PowerSpectrum.Thresh = 80:99;
        OptimalFreqCutOff.PowerSpectrum.LowerBandPassOpt = nan(1,length(80:99));
        for tt = 1:length(OptimalFreqCutOff.PowerSpectrum.Thresh)
            Thresh = OptimalFreqCutOff.PowerSpectrum.Thresh(tt);
            IndMax=find(OptimalFreqCutOff.PowerSpectrum.Pxx_Perc > Thresh,1);
            if numel(IndMax)~=0
                OptimalFreqCutOff.PowerSpectrum.LowerBandPassOpt(tt) = F(IndMax);
            end
        end
        
        % Identify optimal frequency cut-off (point where an increase of 10Hz per 1%
        % increase of cumsum power)
        dPower = OptimalFreqCutOff.PowerSpectrum.Pxx_Perc(2:end) - OptimalFreqCutOff.PowerSpectrum.Pxx_Perc(1:(end-1));
        dF = F(2:end)-F(1:(end-1));
        CutOffIndOpt = find(dF./dPower>10,1);
        OptimalFreqCutOff.PowerSpectrum.LowerBandPassOpt_as4slope = F(CutOffIndOpt-1);
    end
    
    % save data for each semantic cell in its own file
    if PrevData
        save(Calfilename_local,'MatfilePath', 'OptimalFreqCutOff', '-append');
    else
        save(Calfilename_local,'MatfilePath', 'OptimalFreqCutOff');
    end
    
    %     OLD CODE
    %     ParamModel.NeuroRes = 'count_gaussfiltered';
    %Data processing
    %     [HalfTrain1, HalfTrain2, NumTrials]=organiz_data4coherence(Res.Trials_GaussFiltered(DataSel),Res.PSTH_GaussFiltered(DataSel),ParamModel);
    %     [CoherenceStruct]=compute_coherence_mean(HalfTrain1, HalfTrain2,Res.Response_samprate);
    % Compute coherence on gaussian filtered spike trains
    %     OptimalFreqCutOff.CoherenceGaussFilt = CoherenceStruct.freqCutoff;
    
    % Calculate the frequency of the gaussian filtered PSTH below which 99%
    % of the spectrum power density is contained
    %    OptimalFreqCutOff.Thresh = 80:99;
    %     for kk=1:5
    %         % First retrieve the PSTH calculated with the same # of nearest
    %         % neighbour Ntrial/d where d=1:Ntrials
    %         SignalTot_local = nan(nvoc,size(Res.PSTH_GaussFiltered{1},2));
    %         if kk<5 % Treating Neigh = NT/2, NT/3, NT/4, NT/5
    %             for vv=1:nvoc
    %                 NT = size(Res.PSTH_GaussFiltered{DataSel(vv)},1);
    %                 SignalTot_local(vv,:) = Res.PSTH_GaussFiltered{DataSel(vv)}(NT-kk,:);
    %             end
    %         else % Treating Neigh =1 = NT/NT
    %             for vv=1:nvoc
    %                 SignalTot_local(vv,:) = Res.PSTH_GaussFiltered{DataSel(vv)}(1,:);
    %             end
    %         end
    %         SignalTot_1dim=reshape(SignalTot_local',[size(SignalTot_local,1)*size(SignalTot_local,2),1]);
    %
    %         % Calculate the power spectrum of the signal
    %         Window = 0.2*Res.Response_samprate;
    %         [Pxx,F] = pwelch(SignalTot_1dim, Window, [],[],10000);
    %         Pxx_Perc = 100*cumsum(Pxx / sum(Pxx));
    %
    %         % Identify the frequency cut-off corresponding to the list of
    %         % thresholds
    %         OptimalFreqCutOff.(sprintf('PowerSpectrumDensityKth%d',kk)) = nan(1,length(80:99));
    %         for tt = 1:length(OptimalFreqCutOff.Thresh)
    %             Thresh = OptimalFreqCutOff.Thresh(tt);
    %             IndMax=find(Pxx_Perc > Thresh);
    %             OptimalFreqCutOff.(sprintf('PowerSpectrumDensityKth%d',kk))(tt) = F(IndMax(1));
    %         end
    %     end
    
    % According to this code 10ms is the best size for 97% of cells see
    % fig BestPSTHBin.fig
    % At that window size, the values of the FanoFactor over cells is very
    % close to 1. see fig PoissonFanoFactor.fig
end    

%% Estimate Poisson assumption for data at the choosen bining
if SWITCH.FanoFactor
    [PG_Index,FanoFactor_Index, Wins,Ymean,Yvar] = PoissonGaussianNeuralResponses(Res.Trials(DataSel),ParamModel,SWITCH,Cellname);
   % [PG_Index,FanoFactor_Index, Wins,~,~] = PoissonGaussianNeuralResponses(Res.Trials(DataSel),ParamModel,SWITCH,Cellname);
    if PrevData
        save(Calfilename_local,'PG_Index', 'FanoFactor_Index','Wins','-append');
    else
        save(Calfilename_local,'PG_Index', 'FanoFactor_Index','Wins');
    end
    
end

%% Calculate information about stimuli along time
if SWITCH.InfoCal
    ParamModel.Response_samprate = Res.Response_samprate;
    
    
    % Find out the number of trials per stimulus and feed-in PSTH and
    % PSTH_JackKnife
    Ntrials_perstim = nan(length(DataSel),1);
    for vv=1:nvoc
            Ntrials_perstim(vv) = length(Res.Trials{DataSel(vv)});
    end
    ParamModel.Mean_Ntrials_perstim = mean(Ntrials_perstim);
    if mean(Ntrials_perstim - 1) ~= (mean(Ntrials_perstim) - 1)
        fprintf('check number of trials something weirds going on...\n')
        return
    end
        
    % add as input to the parameters the set of indices for predetermined
    % JK sets
    ParamModel.SetIndices_JK = Res.SetIndices_JK;
    
    % Calculate information
    %Calfilename_localKth = sprintf('%s_Kth%d_%s',calfilename_local(1:(end-4)),Kth_i,calfilename_local((end-4):end));
    [ParamModel, Data, InputData, Wins]=info_callsemantic(Res.PSTH_KDE_Filtered(DataSel),Res.JackKnife_KDE_Filtered(DataSel),Res.VocType(DataSel), ParamModel,Calfilename_local);
    
    % save the wave sections
    Data.SectionWave = Res.SectionWave(DataSel);
    
    if exist(Calfilename_local, 'file') == 2
        fprintf(1,'appending to the file\n');
        save(Calfilename_local,'Data', 'InputData','Wins','ParamModel','-append');
    else
        save(Calfilename_local,'Data', 'InputData','Wins','ParamModel');
    end
    
    ElapsedTime = toc(TimerVal);
    Days = floor(ElapsedTime/(60*60*24));
    ETRem = ElapsedTime - Days*60*60*24;
    Hours = floor(ETRem/(60*60));
    ETRem = ETRem - Hours*60*60;
    Minutes = floor(ETRem/60);
    ETRem = ETRem-60*Minutes;
    fprintf(1,'Code run for %d days %dh %dmin and %dsec\n',Days,Hours,Minutes,ETRem);
    
end

%% Calculate cumulative information about stimuli along time
if SWITCH.CumInfoCal
    if exist(Calfilename_local, 'file') == 2
        fprintf(1,'loading the data\n');
        load(Calfilename_local,'Data', 'InputData');
    else
        fprintf(1,'No Data available for this cell, run the calculation of information first!\n');
    end
    if nargin<2
        fprintf('No input for the type of cumulative information!!!\n');
        return
    else
        ParamModel.CIType = CIType;
    end
    ParamModel.Response_samprate = Res.Response_samprate;
    ParamModel.MaxNumSamples_MCopt_Cum_Info = 5*10^6;
    
    % Find out the number of trials per stimulus
    Ntrials_perstim = nan(length(DataSel),1);
    % Find out # of trials per stim
    for vv=1:nvoc
            Ntrials_perstim(vv) = length(Res.Trials{DataSel(vv)});
    end
    ParamModel.Mean_Ntrials_perstim = mean(Ntrials_perstim);
    if mean(Ntrials_perstim - 1) ~= (mean(Ntrials_perstim) - 1)
        fprintf('check number of trials something weirds going on...\n')
        return
    end
        
    % add as input to the parameters the set of indices for predetermined
    % JK sets
    ParamModel.SetIndices_JK = Res.SetIndices_JK;
    
    % Calculate cumulative information
    [Calfilename_local_final]=cuminfo_callsemantic(Data, ParamModel,  Calfilename_local);
    fprintf('Data saved to %s\n', Calfilename_local_final)
end
   

%%
%  if Savio
%      KeepLocalFile = 1;
%      [Status1]=transfertoTdrive_savio(calfilename_local,calfilename_final,KeepLocalFile);
%      [Status2]=transfertoTdrive_savio(outfilename_local,[OutputDir_final '/'],KeepLocalFile);
%      if ~(Status1 || Status2)
%          system(['mv ' OutputDirEx_local '/JobToDoSavio/Ex*' Res.Site '*.txt ' OutputDirEx_local '/JobDoneSavio/'])
%      end
%      fprintf(1,'Ready to quit');
%      quit
%  end

end

