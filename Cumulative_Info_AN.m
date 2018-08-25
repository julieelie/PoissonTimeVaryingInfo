function Cumulative_Info_AN(InputFile,kk)
% Cell name
Cell = 'TI';
% on the cluster
addpath(genpath('/global/home/users/jelie/CODE/SingleUnitModels'));
addpath(genpath('/global/home/usres/jelie/CODE/GeneralCode'));
addpath(genpath('/global/home/users/jelie/CODE/tlab/src'));
rmpath(genpath('/global/home/usres/jelie/CODE/tlab/src/hedi'));

%% Configure Parallel computing
if ~isempty(strfind(getenv('HOSTNAME'),'.savio')) || ~isempty(strfind(getenv('HOSTNAME'),'.brc'))
    MyParPool = parpool(str2num(getenv('SLURM_CPUS_ON_NODE')),'IdleTimeout', Inf);
    system('mkdir -p /global/scratch/$USER/$SLURM_JOB_ID')
    [~,JobID] = system('echo $SLURM_JOB_ID');
    parcluster.JobStorageLocation = ['/global/scratch/jelie/' JobID];    
end

%% load data and start calculations
Dir_local='/global/scratch/jelie/MatFiles/';
DataCell=loadfromTdrive_savio(InputFile, Dir_local,1);

Cum_boot=20;
fprintf('------------------------------------------------------\n')
fprintf('---------- Temporal code Dataset-----------\n')
fprintf('------------------------------------------------------\n')
fprintf('**************** Calculating cumulative information *****************\n')
Nb_Kth = min([length(DataCell.Kth_Neigh_JK) length(DataCell.Kth_Neigh)]);

P_YgivenS = DataCell.(sprintf('Kth%d',kk)).P_YgivenS;
P_YgivenS_BootJK = DataCell.(sprintf('Kth%d',kk)).P_YgivenS_BootJK;
Bin_Trials = DataCell.(sprintf('Kth%d',kk)).Bin_Trials;

%% Calculate the Monte Carlo estimation with optimum number of samples
MaxMCParameter = 5*10^6;
ConvThresh = 0.2;
IncrMCParameter = 10^5;
NTrials = size(Bin_Trials,2);
Nb_Win = length(P_YgivenS);
Icum_EstMonteCarloOpt = nan(1,Nb_Win);
Icum_EstMonteCarloOpt_err = nan(1,Nb_Win);
Icum_EstMonteCarloOpt_bcorr = nan(1,Nb_Win);
Icum_EstMonteCarloOpt(1) = DataCell.(sprintf('Kth%d',kk)).Info(1); % Initializing the first value of cumulative info
Icum_EstMonteCarloOpt_bcorr(1) = DataCell.(sprintf('Kth%d',kk)).Info(1); % Initializing the first value of cumulative info

save(sprintf('%sInfoCumInfoSpikeCount_AN_JK_KNeigh%dms%d%s.mat',Dir_local,DataCell.StepInfo,kk,Cell),'-struct','DataCell');

MC_Samp = nan(1,Nb_Win);
tstart2 = tic;
fprintf('**** Monte Carlo with optimal # samples and Jackknife kk=%d/%d *****\n', kk, Nb_Kth);
for tt=2:Nb_Win
    tstart = tic;
    fprintf('Time point %d/%d\n', tt, Nb_Win);
    P_YgivenS_local = P_YgivenS(1:tt);
    P_YgivenS_JK_local = P_YgivenS_BootJK(1:Cum_boot,1:tt);
    
    % Monte Carlo estimation with full memory
    [Icum_EstMonteCarloOpt(tt),~, ~, Icum_EstMonteCarloOpt_bcorr(tt),Icum_EstMonteCarloOpt_err(tt),MC_Samp(tt)]=cumulative_info_poisson_model_calculus_MCJK(P_YgivenS_local, P_YgivenS_JK_local,NTrials,'ConvThresh',ConvThresh, 'MaxMCParameter',MaxMCParameter, 'IncrMCParameter', IncrMCParameter);
    fprintf('# MC samples %d Error:%.2f\n', MC_Samp(tt),Icum_EstMonteCarloOpt_err(tt));
    telapsed = toc(tstart);
    fprintf('Elapsed time: %d s\n', telapsed)
    if Icum_EstMonteCarloOpt_err(tt)>(3*ConvThresh)
        fprintf('Error is too high at this point stop here: %.2f\n', Icum_EstMonteCarloOpt_err(tt))
        return
    end
    DataCell.(sprintf('Kth%d',kk)).Nb_Win = Nb_Win;
    DataCell.(sprintf('Kth%d',kk)).Cum_boot = Cum_boot;
    DataCell.(sprintf('Kth%d',kk)).Icum_EstMonteCarloOpt = Icum_EstMonteCarloOpt;
    DataCell.(sprintf('Kth%d',kk)).Icum_EstMonteCarloOpt_bcorr = Icum_EstMonteCarloOpt_bcorr;
    DataCell.(sprintf('Kth%d',kk)).Icum_EstMonteCarloOpt_err = Icum_EstMonteCarloOpt_err;
    DataCell.(sprintf('Kth%d',kk)).MC_Samp =MC_Samp;
    %if tt/10 == fix(tt/10)
    save(sprintf('%sInfoCumInfoSpikeCount_AN_JK_KNeigh%dms%d%s.mat',Dir_local,DataCell.StepInfo,kk,Cell),'-struct','DataCell','-append');
%     save(sprintf('%sInfoCumInfoSpikeCount_AN_JK_KNeigh%d%s.mat',Dir_local,kk,Cell),'-struct','DataCell','-append');
    %end
    telapsed = toc(tstart);
    fprintf('Elapsed time after saving line: %d s\n', telapsed)
end
telapsed2 = toc(tstart2);
fprintf('MC Opt total elapsed time: %d s\n', telapsed2)
quit

end