function [Icum_EstMonteCarloOpt, Icum_EstMonteCarloOpt_bcorr, Icum_EstMonteCarloOpt_err, MC_Samp] = cumulative_info_poisson_model_MCJK_wrapper(P_YgivenS, P_YgivenS_BootJK, NTrials,Nb_Win, Cum_boot, MaxMCParameter, FirstStep, IncrMCParameter, ConvThresh, Verbose)
%% Treat input arguments
if nargin<10
    Verbose = 0;
end
if nargin<9
    ConvThresh = 0.2;
end
if nargin<8
    IncrMCParameter = 10^5;
end
if nargin<7
    FirstStep=1;
end
if nargin<6
    MaxMCParameter = 5*10^6;
end
if nargin<5
    Cum_boot = length(P_YgivenS_BootJK);
end
if nargin<4
    Nb_Win = length(P_YgivenS);
end

%% Configure Parallel computing
if ~isempty(strfind(getenv('HOSTNAME'),'.savio')) || ~isempty(strfind(getenv('HOSTNAME'),'.brc'))
    delete(gcp('nocreate'))
    MyParPool = parpool(str2num(getenv('SLURM_CPUS_ON_NODE')),'IdleTimeout', Inf);
    system('mkdir -p /global/scratch/$USER/$SLURM_JOB_ID')
    [~,JobID] = system('echo $SLURM_JOB_ID');
    parcluster.JobStorageLocation = ['/global/scratch/jelie/' JobID];    
end

%% initialize output variables
Icum_EstMonteCarloOpt = nan(1,Nb_Win);
Icum_EstMonteCarloOpt_err = nan(1,Nb_Win);
Icum_EstMonteCarloOpt_bcorr = nan(1,Nb_Win);
MC_Samp = nan(1,Nb_Win);

%% Loop through time bins and calculate cumulative information as long as...
... the error does not hit the threshold = 3 times the convergence error value ...
    ...(3*ConvThresh) that cumulative_info_poisson_model_calculus_MCJK is trying to match...
    ... with as many samples as necessary, up to MaxMCParameter samples
tstart2 = tic;
fprintf('**** Monte Carlo with optimal # samples and Jackknife from time step %d *****\n', FirstStep);
Error_local = 0;
tt=FirstStep;
while ((Error_local<=(3*ConvThresh)) && (tt<Nb_Win))
    tt= tt+1;
    tstart = tic;
    fprintf('Time point %d/%d\n', tt, Nb_Win);
    P_YgivenS_local = P_YgivenS(1:tt);
    P_YgivenS_JK_local = cell(1,Cum_boot);
    for bb=1:Cum_boot
        P_YgivenS_JK_local{bb} = P_YgivenS_BootJK{bb}(:,1:tt);
    end
    
    % Monte Carlo estimation with full memory
    [Icum_EstMonteCarloOpt(tt),~, ~, Icum_EstMonteCarloOpt_bcorr(tt),Icum_EstMonteCarloOpt_err(tt),MC_Samp(tt)]=cumulative_info_poisson_model_calculus_MCJK(P_YgivenS_local, P_YgivenS_JK_local,NTrials,'ConvThresh',ConvThresh, 'MaxMCParameter',MaxMCParameter, 'IncrMCParameter', IncrMCParameter, 'Verbose',Verbose);
    fprintf('# MC value: %.2f MC samples %d Error:%.2f\n', Icum_EstMonteCarloOpt(tt), MC_Samp(tt),Icum_EstMonteCarloOpt_err(tt));
    telapsed = toc(tstart);
    fprintf('Elapsed time: %d s\n', telapsed)
    Error_local=Icum_EstMonteCarloOpt_err(tt);
end
fprintf('Calculations stop at %d with an error of: %.2f\n', tt, Icum_EstMonteCarloOpt_err(tt))
telapsed2 = toc(tstart2);
fprintf('MC Opt total elapsed time: %d s\n', telapsed2)

%% get rid of temporary files for parallel computing
if ~isempty(strfind(getenv('HOSTNAME'),'.savio')) || ~isempty(strfind(getenv('HOSTNAME'),'.brc'))
    delete(MyParPool);
    system(['rm -r ' parcluster.JobStorageLocation])
end
end