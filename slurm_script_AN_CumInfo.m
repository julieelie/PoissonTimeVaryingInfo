function slurm_script_AN_CumInfo(Cell)
addpath(genpath('/global/home/users/jelie/CODE/SingleUnitModels'));
addpath(genpath('/global/home/users/jelie/CODE/GeneralCode'));
addpath(genpath('/global/home/users/jelie/CODE/tlab/src'));
rmpath(genpath('/global/home/users/jelie/CODE/tlab/src/hedi'));
Path2Data='/auto/tdrive/julie/NeuralData/SemanticInfoPoisson/';
Storage_path = '/global/scratch/jelie/MatFiles/';

if ~isempty(strfind(getenv('HOSTNAME'),'.savio')) || ~isempty(strfind(getenv('HOSTNAME'),'.brc'))
        MyParPool = parpool(str2num(getenv('SLURM_CPUS_ON_NODE')),'IdleTimeout', Inf);
        system('mkdir -p /global/scratch/$USER/$SLURM_JOB_ID')
        [~,JobID] = system('echo $SLURM_JOB_ID');
        parcluster.JobStorageLocation = ['/global/scratch/jelie/' JobID];
end

Cum_boot=10;
fprintf('------------------------------------------------------\n')
fprintf('---------- %s cell Dataset-----------\n', Cell)
fprintf('------------------------------------------------------\n')
fprintf('**************** Calculating cumulative information *****************\n')
From_file = sprintf('%sInfoCumInfoSpikeCount_AN_JK_KDE_%s.mat',Path2Data,Cell);
KeepF=1;

[Data]=loadfromTdrive_savio(From_file,Storage_path, KeepF);
P_YgivenS = Data.P_YgivenS;
P_YgivenS_BootJK = Data.P_YgivenS_BootJK;
Bin_Trials = Data.Bin_Trials;
clear Data

%% Calculate the Monte Carlo estimation with optimum number of samples
MaxMCParameter = 5*10^6;
ConvThresh = 0.2;
IncrMCParameter = 10^5;
NTrials = size(Bin_Trials,2);
clear Bin_Trials
Nb_Win = length(P_YgivenS);
Icum_EstMonteCarloOpt = nan(1,Nb_Win);
Icum_EstMonteCarloOpt_err = nan(1,Nb_Win);
Icum_EstMonteCarloOpt_bcorr = nan(1,Nb_Win);
MC_Samp = nan(1,Nb_Win);
tstart2 = tic;
fprintf('**** Monte Carlo with optimal # samples and Jackknife *****\n');
Error_local = 0;
tt=1;
while ((Error_local<=(3*ConvThresh)) && (tt<Nb_Win))
    tt= tt+1;
    tstart = tic;
    fprintf('Time point %d/%d\n', tt, Nb_Win);
    P_YgivenS_local = P_YgivenS(1:tt);
    P_YgivenS_JK_local = cell(1:Cum_boot);
    for bb=1:Cum_boot
        P_YgivenS_JK_local{bb} = P_YgivenS_BootJK{bb}(:,1:tt);
    end
    
    % Monte Carlo estimation with full memory
    [Icum_EstMonteCarloOpt(tt),~, ~, Icum_EstMonteCarloOpt_bcorr(tt),Icum_EstMonteCarloOpt_err(tt),MC_Samp(tt)]=cumulative_info_poisson_model_calculus_MCJK(P_YgivenS_local, P_YgivenS_JK_local,NTrials,'ConvThresh',ConvThresh, 'MaxMCParameter',MaxMCParameter, 'IncrMCParameter', IncrMCParameter);
    fprintf('# MC samples %d Error:%.2f\n', MC_Samp(tt),Icum_EstMonteCarloOpt_err(tt));
    telapsed = toc(tstart);
    fprintf('Elapsed time: %d s\n', telapsed)
    Error_local=Icum_EstMonteCarloOpt_err(tt);
end
fprintf('Calculations stop at %d with an error of: %.2f\n', tt, Icum_EstMonteCarloOpt_err(tt))
telapsed2 = toc(tstart2);
fprintf('MC Opt total elapsed time: %d s\n', telapsed2)
load(sprintf('%sInfoCumInfoSpikeCount_AN_JK_KDE_%s.mat',Storage_path,Cell),'Info_bcorr');
Icum_EstMonteCarloOpt(1) = Info_bcorr(1); % Initializing the first value of cumulative info
Icum_EstMonteCarloOpt_bcorr(1) = Info_bcorr(1); % Initializing the first value of cumulative info
save(sprintf('%sInfoCumInfoSpikeCount_AN_JK_KDE_%s.mat',Storage_path,Cell),'Nb_Win','Cum_boot','Icum_EstMonteCarloOpt','Icum_EstMonteCarloOpt_bcorr','Icum_EstMonteCarloOpt_err','MC_Samp','-append');
clear P_Y* Icum*

%% Calculate the cumulative information with Markov and exact calculation with 5 bin memory
fprintf('**** Exact calculation and Markov with 5 bins memory *****\n');
load(sprintf('%sInfoCumInfoSpikeCount_AN_JK_KDE_%s.mat',Storage_path,Cell),'P_YgivenS');
Nb_Win = length(P_YgivenS);
Icum_ExactMem0_5 = nan(1,Nb_Win);
Icum_EstMarkov5 = Icum_ExactMem0_5;
HY_Markov5 = nan(1,Nb_Win);
for tt=2:Nb_Win
    tstart = tic;
    fprintf('Time point %d/%d\n', tt, Nb_Win);
    P_YgivenS_local = P_YgivenS(1:tt);
    
    % Exact calculation with 5 bins memory in the past (5 bins total)
    [Icum_ExactMem0_5(tt),~]=info_cumulative_model_Calculus(P_YgivenS_local,'CalMode','Exact_Mem', 'Exact_history',5);
    if tt==2
        % Markov chain estimation 5 bins memory
        [Icum_EstMarkov5(tt),HY_Markov5(tt),~]=info_cumulative_model_Calculus(P_YgivenS_local,'CalMode','MarkovChain', 'MarkovParameters',[5,1]);
    else
        % Markov chain estimation 5 bins memory
        [Icum_EstMarkov5(tt),HY_Markov5(tt),~]=info_cumulative_model_Calculus(P_YgivenS_local,'CalMode','MarkovChain', 'MarkovParameters',[5,1], 'HY_old', HY_Markov5(tt-1));
    end
    telapsed = toc(tstart);
    fprintf('Markov + Exact calculation: total elapsed time: %d s\n', telapsed)
end
load(sprintf('%sInfoCumInfoSpikeCount_AN_JK_KDE_%s.mat',Storage_path,Cell),'Info_bcorr');
Icum_ExactMem0_5(1) =  Info_bcorr(1);
Icum_EstMarkov5(1) = Info_bcorr(1);
save(sprintf('%sInfoCumInfoSpikeCount_AN_JK_KDE_%s.mat',Storage_path,Cell),'Icum_ExactMem0_5','Icum_EstMarkov5','-append');
clear P_Y* Icum* HY*

% jackknife the calculations
fprintf('**** Jackknife Exact calculation and Markov with 5 bins memory *****\n');
load(sprintf('%sInfoCumInfoSpikeCount_AN_JK_KDE_%s.mat',Storage_path,Cell),'P_YgivenS_BootJK', 'Nb_Win','Cum_boot','NTrials','Icum_ExactMem0_5','Icum_EstMarkov5');
Icum_ExactMem0_5_JK_mean = nan(Cum_boot,Nb_Win);
Icum_ExactMem0_5_JK_var = nan(Cum_boot,Nb_Win);
Icum_EstMarkov5_JK_mean = nan(Cum_boot,Nb_Win);
Icum_EstMarkov5_JK_var = nan(Cum_boot,Nb_Win);

parfor bb=1:Cum_boot
    Nb_Sets = size(P_YgivenS_BootJK{bb},1);
    Icum_ExactMem0_5_JK_Set = nan(Nb_Sets, Nb_Win);
    Icum_EstMarkov5_JK_Set = nan(Nb_Sets, Nb_Win);
    for jkk=1:Nb_Sets
        P_YgivenS_JK_local = P_YgivenS_BootJK{bb}(jkk,:);
        HY_Markov5 = nan(1,Nb_Win);
        for tt=2:Nb_Win
            tstart = tic;
            fprintf('Bootstrap %d/%d Time point %d/%d\n',bb,Cum_boot, tt, Nb_Win);
        
            P_YgivenS_local = P_YgivenS_JK_local(1:tt);
    
            % Exact calculation with 5 bins memory
            [Icum_ExactMem0_5_JK_Set(jkk,tt),~]=info_cumulative_model_Calculus(P_YgivenS_local,'CalMode','Exact_Mem', 'Exact_history',5);
            if tt==2
                % Markov chain estimation 5 bins memory
                [Icum_EstMarkov5_JK_Set(jkk,tt),HY_Markov5(tt),~]=info_cumulative_model_Calculus(P_YgivenS_local,'CalMode','MarkovChain', 'MarkovParameters',[5,1]);
            else
                % Markov chain estimation 5 bins memory
                [Icum_EstMarkov5_JK_Set(jkk,tt),HY_Markov5(tt),~]=info_cumulative_model_Calculus(P_YgivenS_local,'CalMode','MarkovChain', 'MarkovParameters',[5,1], 'HY_old', HY_Markov5(tt-1));
            end
            telapsed = toc(tstart);
            fprintf('Markov + Exact calculation: total elapsed time: %d s\n(Bootstrap %d/%d Time point %d/%d)\n',telapsed,bb,Cum_boot, tt, Nb_Win)
        end
    end
    % Calculate the Jack-knife biais corrected values of information and their errors
    Icum_ExactMem0_5_JK_Setcorrected = NTrials .* repmat(Icum_ExactMem0_5,Nb_Sets,1) - (NTrials-1) .* Icum_ExactMem0_5_JK_Set;
    Icum_ExactMem0_5_JK_mean(bb,:) = mean(Icum_ExactMem0_5_JK_Setcorrected,1);
    Icum_ExactMem0_5_JK_var(bb,:) = var(Icum_ExactMem0_5_JK_Setcorrected,0,1);
    Icum_EstMarkov5_JK_Setcorrected = NTrials .* repmat(Icum_EstMarkov5,Nb_Sets,1) - (NTrials-1) .* Icum_EstMarkov5_JK_Set;
    Icum_EstMarkov5_JK_mean(bb,:) = mean(Icum_EstMarkov5_JK_Setcorrected,1);
    Icum_EstMarkov5_JK_var(bb,:) = var(Icum_EstMarkov5_JK_Setcorrected,0,1);
end

load(sprintf('%sInfoCumInfoSpikeCount_AN_JK_KDE_%s.mat',Storage_path,Cell),'Info_bcorr');
Icum_ExactMem0_5_JK_mean(:,1) = repmat(Info_bcorr(1),Cum_boot,1);
Icum_ExactMem0_5_bcorr = mean(Icum_ExactMem0_5_JK_mean,1);
Icum_ExactMem0_5_err = (mean(Icum_ExactMem0_5_JK_var,1)).^0.5;

Icum_EstMarkov5_JK_mean(:,1) = repmat(Info_bcorr(1),Cum_boot,1);
Icum_EstMarkov5_bcorr = mean(Icum_EstMarkov5_JK_mean,1);
Icum_EstMarkov5_err = (mean(Icum_EstMarkov5_JK_var,1)).^0.5;

save(sprintf('%sInfoCumInfoSpikeCount_AN_JK_KDE_%s.mat',Storage_path,Cell),'Icum_ExactMem0_5_JK_mean','Icum_EstMarkov5_JK_mean','Icum_ExactMem0_5_JK_var','Icum_EstMarkov5_JK_var','Icum_ExactMem0_5_bcorr','Icum_EstMarkov5_bcorr','Icum_ExactMem0_5_err','Icum_EstMarkov5_err','-append');
clear P_Y* Icum* HY*

%% Calculate theoretical values if the spike rate was exactly known
fprintf('**** Theoretical values: Monte Carlo 10^6, Exact calculation and Markov with 5 bins memory *****\n');
load(sprintf('%sInfoCumInfoSpikeCount_AN_JK_KDE_%s.mat',Storage_path,Cell),'P_YgivenS_Theo', 'Nb_Win','Cum_boot');
Icum_EstMonteCarlo6_Theo = nan(1,Nb_Win);
Icum_ExactMem0_5_Theo = nan(1,Nb_Win);
Icum_EstMarkov5_Theo = nan(1,Nb_Win);
HY_Markov5 = nan(1,Nb_Win);
for tt=2:Nb_Win
    tstart = tic;
    fprintf('Time point %d/%d\n', tt, Nb_Win);
    P_YgivenS_local = P_YgivenS_Theo(1:tt);
    % Monte Carlo estimation with full memory
    [Icum_EstMonteCarlo6_local,~]=info_cumulative_model_Calculus(P_YgivenS_local,'CalMode','MonteCarlo', 'MCParameter',10^6);
    Icum_EstMonteCarlo6_Theo(tt) = Icum_EstMonteCarlo6_local(1);
    
    % Exact calculation with 5 bins memory
    [Icum_ExactMem0_5_Theo(tt),~]=info_cumulative_model_Calculus(P_YgivenS_local,'CalMode','Exact_Mem', 'Exact_history',5);
    if tt==2
        % Markov chain estimation 5 bins memory
        [Icum_EstMarkov5_Theo(tt),HY_Markov5(tt),~]=info_cumulative_model_Calculus(P_YgivenS_local,'CalMode','MarkovChain', 'MarkovParameters',[5,1]);
    else
        % Markov chain estimation 5 bins memory
        [Icum_EstMarkov5_Theo(tt),HY_Markov5(tt),~]=info_cumulative_model_Calculus(P_YgivenS_local,'CalMode','MarkovChain', 'MarkovParameters',[5,1], 'HY_old', HY_Markov5(tt-1));
    end
    telapsed = toc(tstart);
    fprintf('Markov + Exact calculation: total elapsed time: %d s\n', telapsed)
end
load(sprintf('%sInfoCumInfoSpikeCount_AN_JK_KDE_%s.mat',Storage_path,Cell),'Info_bcorr');
Icum_EstMonteCarlo6_Theo(1) = Info_bcorr(1); % Initializing the first value of cumulative info
Icum_ExactMem0_5_Theo(1) =  Info_bcorr(1);
Icum_EstMarkov5_Theo(1) = Info_bcorr(1);
save(sprintf('%sInfoCumInfoSpikeCount_AN_JK_KDE_%s.mat',Storage_path,Cell),'Icum_EstMonteCarlo6_Theo','Icum_ExactMem0_5_Theo','Icum_EstMarkov5_Theo', '-append');
clear P_Y* Icum* HY*
end