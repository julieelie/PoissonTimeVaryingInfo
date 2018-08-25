function slurm_script_AN_CumInfo_patchcorrection(Cell)
addpath(genpath('/global/home/users/jelie/CODE/SingleUnitModels'));
addpath(genpath('/global/home/users/jelie/CODE/GeneralCode'));
addpath(genpath('/global/home/users/jelie/CODE/tlab/src'));
rmpath(genpath('/global/home/users/jelie/CODE/tlab/src/hedi'));
Storage_path = '/global/scratch/jelie/MatFiles/';

fprintf('------------------------------------------------------\n')
fprintf('---------- %s cell Dataset-----------\n', Cell)
fprintf('------------------------------------------------------\n')
fprintf('**************** Calculating cumulative information *****************\n')
Cum_boot=10;

%% jackknife the calculations
fprintf('**** Jackknife Exact calculation and Markov with 5 bins memory *****\n');
load(sprintf('%sInfoCumInfoSpikeCount_AN_JK_KDE_%s.mat',Storage_path,Cell), 'Nb_Win','Cum_boot','NTrials');
InputData = load(sprintf('%sInfoCumInfoSpikeCount_AN_JK_KDE_%s.mat',Storage_path,Cell),'P_YgivenS','P_YgivenS_BootJK');
Nb_Win = length(InputData.P_YgivenS);
Icum_ExactMem0_5_JK_mean = nan(Cum_boot,Nb_Win);
Icum_ExactMem0_5_JK_var = nan(Cum_boot,Nb_Win);
Icum_EstMarkov5_JK_mean = nan(Cum_boot,Nb_Win);
Icum_EstMarkov5_JK_var = nan(Cum_boot,Nb_Win);
NTrials=10;
parfor bb=1:Cum_boot
    Nb_Sets = size(InputData.P_YgivenS_BootJK{bb},1);
    Icum_ExactMem0_5_JK_Set = nan(Nb_Sets, Nb_Win);
    Icum_EstMarkov5_JK_Set = nan(Nb_Sets, Nb_Win);
    for jkk=1:Nb_Sets
        P_YgivenS_JK_local = InputData.P_YgivenS_BootJK{bb}(jkk,:);
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
    LocalInput = load(sprintf('%sInfoCumInfoSpikeCount_AN_JK_KDE_%s.mat',Storage_path,Cell),'Icum_ExactMem0_5','Icum_EstMarkov5');
    % Calculate the Jack-knife biais corrected values of information and their errors
    Icum_ExactMem0_5_JK_Setcorrected = NTrials .* repmat(LocalInput.Icum_ExactMem0_5,Nb_Sets,1) - (NTrials-1) .* Icum_ExactMem0_5_JK_Set;
    Icum_ExactMem0_5_JK_mean(bb,:) = mean(Icum_ExactMem0_5_JK_Setcorrected,1);
    Icum_ExactMem0_5_JK_var(bb,:) = var(Icum_ExactMem0_5_JK_Setcorrected,0,1);
    Icum_EstMarkov5_JK_Setcorrected = NTrials .* repmat(LocalInput.Icum_EstMarkov5,Nb_Sets,1) - (NTrials-1) .* Icum_EstMarkov5_JK_Set;
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
Nb_Win = length(P_YgivenS_Theo);
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