function [Model] = info_cumulative_wrapper(ParamModel,SWITCH,Model,mm,x_stim_indices_wholeset)
fprintf('Pre-process data for cumulative information calculations\n')
X_stim_indices_wholeset = cell(10,1);

mm_local=0;
%Pre-process data for a parfor loop
P_YgivenS_allModel = cell(10,1);
if ParamModel.ModelChoice(1) && ~SWITCH.AllAlpha
    % ACoustic Model
    mm_local=mm_local+1;
    P_YgivenS_allModel{mm_local}=Model.Acoustic.P_YgivenS_all1(1:mm,1);
    X_stim_indices_wholeset{mm_local}=x_stim_indices_wholeset(1:mm);
end
if ParamModel.ModelChoice(2) && ~SWITCH.AllAlpha
    % Semantic Model
    mm_local=mm_local+1;
    P_YgivenS_allModel{mm_local}=Model.Semantic.P_YgivenS_all1(1:mm,1);
    X_stim_indices_wholeset{mm_local}=x_stim_indices_wholeset(1:mm);
    % Semantic Model on dataset of the last window only
%     mm_local=mm_local+1;
%     P_YgivenS_allModel{mm_local}=Model.Semantic.P_YgivenS_all1_enddataset(1:mm,1);
%     X_stim_indices_wholeset{mm_local}=cell(mm,1);
%     X_stim_indices_wholeset{mm_local}(:)=x_stim_indices_wholeset(end);
end
if ParamModel.ModelChoice(4) && ~SWITCH.AllAlpha
    % AcSemAc
    mm_local=mm_local+1;
    P_YgivenS_allModel{mm_local}=Model.AcSemAc.P_YgivenS_all1(1:mm,1);
    X_stim_indices_wholeset{mm_local}=x_stim_indices_wholeset(1:mm);
end
if ParamModel.ModelChoice(5) && ~SWITCH.AllAlpha
    % AcSemSem
    mm_local=mm_local+1;
    P_YgivenS_allModel{mm_local}=Model.AcSemSem.P_YgivenS_all1(1:mm,1);
    X_stim_indices_wholeset{mm_local}=x_stim_indices_wholeset(1:mm);
end
% Floor
% mm_local=mm_local+1;
% P_YgivenS_allModel{mm_local}=Model.Floor.P_YgivenS_all1(1:mm,1);
%Ceiling
mm_local=mm_local+1;
P_YgivenS_allModel{mm_local}=Model.Ceiling.P_YgivenS_all1(1:mm,1);
X_stim_indices_wholeset{mm_local}=x_stim_indices_wholeset(1:mm);
%Ceiling on dataset of the last windo only
% mm_local=mm_local+1;
% P_YgivenS_allModel{mm_local}=Model.Ceiling.P_YgivenS_all1_enddataset(1:mm,1);
% X_stim_indices_wholeset{mm_local}=cell(mm,1);
% X_stim_indices_wholeset{mm_local}(:)=x_stim_indices_wholeset(end);

% AR
if SWITCH.AR
    mm_local=mm_local+1;
    P_YgivenS_allModel{mm_local}=Model.AR.P_YgivenS_all1(1:mm,1);
    X_stim_indices_wholeset{mm_local}=x_stim_indices_wholeset(1:mm);
end

Icum_ExactMem0_5 = nan(mm_local,1);
Icum_EstMonteCarlo = nan(2,mm_local);
Icum_EstMonteCarlo2 = nan(2,mm_local);
Icum_EstMonteCarlo3 = nan(2,mm_local);
Icum_EstMonteCarlo4 = nan(2,mm_local);
Icum_EstMarkov2 = nan(mm_local,1);
Icum_EstMarkov3 = nan(mm_local,1);
Icum_EstMarkov4 = nan(mm_local,1);
Icum_EstMarkov5 = nan(mm_local,1);

fprintf('Calculate Cumulative information for all models from win %d\n',1)
for modelrun=1:mm_local
    fprintf('Cumulative info %d/%d\n', modelrun, mm_local);
    % Monte Carlo estimation with full memory
    [Icum_EstMonteCarlo(:,modelrun),~]=info_cumulative_model_Calculus(P_YgivenS_allModel{modelrun},'StimIndicesAll',X_stim_indices_wholeset{modelrun},'StimIndicesLast',X_stim_indices_wholeset{modelrun}{mm},'Model#',modelrun,'CalMode','MonteCarlo', 'MCParameter',ParamModel.NumSamples_MC_Cum_Info(1));
    % Monte Carlo estimation with full memory
    [Icum_EstMonteCarlo4(:,modelrun),~]=info_cumulative_model_Calculus(P_YgivenS_allModel{modelrun},'StimIndicesAll',X_stim_indices_wholeset{modelrun},'StimIndicesLast',X_stim_indices_wholeset{modelrun}{mm},'Model#',modelrun,'CalMode','MonteCarlo', 'MCParameter',ParamModel.NumSamples_MC_Cum_Info(2));
    % Monte Carlo estimation with full memory
    [Icum_EstMonteCarlo3(:,modelrun),~]=info_cumulative_model_Calculus(P_YgivenS_allModel{modelrun},'StimIndicesAll',X_stim_indices_wholeset{modelrun},'StimIndicesLast',X_stim_indices_wholeset{modelrun}{mm},'Model#',modelrun,'CalMode','MonteCarlo', 'MCParameter',ParamModel.NumSamples_MC_Cum_Info(3));
    % Monte Carlo estimation with full memory
    [Icum_EstMonteCarlo2(:,modelrun),~]=info_cumulative_model_Calculus(P_YgivenS_allModel{modelrun},'StimIndicesAll',X_stim_indices_wholeset{modelrun},'StimIndicesLast',X_stim_indices_wholeset{modelrun}{mm},'Model#',modelrun,'CalMode','MonteCarlo', 'MCParameter',ParamModel.NumSamples_MC_Cum_Info(4));
    % Exact calculation with 50 ms memory
    [Icum_ExactMem0_5(modelrun),~]=info_cumulative_model_Calculus(P_YgivenS_allModel{modelrun},'StimIndicesAll',X_stim_indices_wholeset{modelrun},'StimIndicesLast',X_stim_indices_wholeset{modelrun}{mm},'Model#',modelrun,'CalMode','Exact_Mem', 'Exact_history',5);
    % Markov chain estimation 20 ms
    [Icum_EstMarkov2(modelrun),~]=info_cumulative_model_Calculus(P_YgivenS_allModel{modelrun},'StimIndicesAll',X_stim_indices_wholeset{modelrun},'StimIndicesLast',X_stim_indices_wholeset{modelrun}{mm},'Model#',modelrun,'CalMode','MarkovChain', 'MarkovParameters',[2,1]);
    % Markov chain estimation 30 ms
    [Icum_EstMarkov3(modelrun),~]=info_cumulative_model_Calculus(P_YgivenS_allModel{modelrun},'StimIndicesAll',X_stim_indices_wholeset{modelrun},'StimIndicesLast',X_stim_indices_wholeset{modelrun}{mm},'Model#',modelrun,'CalMode','MarkovChain', 'MarkovParameters',[3,1]);
    % Markov chain estimation 40 ms
    [Icum_EstMarkov4(modelrun),~]=info_cumulative_model_Calculus(P_YgivenS_allModel{modelrun},'StimIndicesAll',X_stim_indices_wholeset{modelrun},'StimIndicesLast',X_stim_indices_wholeset{modelrun}{mm},'Model#',modelrun,'CalMode','MarkovChain', 'MarkovParameters',[4,1]);
    % Markov chain estimation 50 ms
    [Icum_EstMarkov5(modelrun),~]=info_cumulative_model_Calculus(P_YgivenS_allModel{modelrun},'StimIndicesAll',X_stim_indices_wholeset{modelrun},'StimIndicesLast',X_stim_indices_wholeset{modelrun}{mm},'Model#',modelrun,'CalMode','MarkovChain', 'MarkovParameters',[5,1]);
end

mm_local=0;
%Post-process data for a parfor loop
if ParamModel.ModelChoice(1) && ~SWITCH.AllAlpha
    fprintf('**CumInfo on Acoustic**\n')
    % ACoustic Model
    mm_local=mm_local+1;
    Model.Acoustic.cum_info.ExactHist5(mm,1)=Icum_ExactMem0_5(mm_local);
    Model.Acoustic.cum_info.EstMonteCarlo(mm,1:2)=Icum_EstMonteCarlo(:,mm_local);
    Model.Acoustic.cum_info.EstMarkov2(mm,1)=Icum_EstMarkov2(mm_local);
    Model.Acoustic.cum_info.EstMarkov3(mm,1)=Icum_EstMarkov3(mm_local);
    Model.Acoustic.cum_info.EstMarkov4(mm,1)=Icum_EstMarkov4(mm_local);
    Model.Acoustic.cum_info.EstMarkov5(mm,1)=Icum_EstMarkov5(mm_local);
end
if ParamModel.ModelChoice(2) && ~SWITCH.AllAlpha
    fprintf('**CumInfo on Semantic**\n')
    % Semantic Model
    mm_local=mm_local+1;
    Model.Semantic.cum_info.current_dataset.ExactHist5(mm,1)=Icum_ExactMem0_5(mm_local);
    Model.Semantic.cum_info.current_dataset.EstMonteCarlo(mm,1:2)=Icum_EstMonteCarlo(:,mm_local);
    Model.Semantic.cum_info.current_dataset.EstMarkov2(mm,1)=Icum_EstMarkov2(mm_local);
    Model.Semantic.cum_info.current_dataset.EstMarkov3(mm,1)=Icum_EstMarkov3(mm_local);
    Model.Semantic.cum_info.current_dataset.EstMarkov4(mm,1)=Icum_EstMarkov4(mm_local);
    Model.Semantic.cum_info.current_dataset.EstMarkov5(mm,1)=Icum_EstMarkov5(mm_local);
%     mm_local=mm_local+1;
%     Model.Semantic.cum_info.end_dataset.ExactHist5(mm,1)=Icum_ExactMem0_5(mm_local);
%     Model.Semantic.cum_info.end_dataset.EstMonteCarlo(mm,1:2)=Icum_EstMonteCarlo(:,mm_local);
%     Model.Semantic.cum_info.end_dataset.EstMarkov2(mm,1)=Icum_EstMarkov2(mm_local);
%     Model.Semantic.cum_info.end_dataset.EstMarkov3(mm,1)=Icum_EstMarkov3(mm_local);
%     Model.Semantic.cum_info.end_dataset.EstMarkov4(mm,1)=Icum_EstMarkov4(mm_local);
%     Model.Semantic.cum_info.end_dataset.EstMarkov5(mm,1)=Icum_EstMarkov5(mm_local);
end
if ParamModel.ModelChoice(4) && ~SWITCH.AllAlpha
    fprintf('**CumInfo on AcSemAc**\n')
    % AcSemAc
    mm_local=mm_local+1;
    Model.AcSemAc.cum_info.ExactHist5(mm,1)=Icum_ExactMem0_5(mm_local);
    Model.AcSemAc.cum_info.EstMonteCarlo(mm,1:2)=Icum_EstMonteCarlo(:,mm_local);
    Model.AcSemAc.cum_info.EstMarkov2(mm,1)=Icum_EstMarkov2(mm_local);
    Model.AcSemAc.cum_info.EstMarkov3(mm,1)=Icum_EstMarkov3(mm_local);
    Model.AcSemAc.cum_info.EstMarkov4(mm,1)=Icum_EstMarkov4(mm_local);
    Model.AcSemAc.cum_info.EstMarkov5(mm,1)=Icum_EstMarkov5(mm_local);
end
if ParamModel.ModelChoice(5) && ~SWITCH.AllAlpha
    fprintf('**CumInfo on AcSemSem**\n')
    % AcSemSem
    mm_local=mm_local+1;
    Model.AcSemSem.cum_info.ExactHist5(mm,1)=Icum_ExactMem0_5(mm_local);
    Model.AcSemSem.cum_info.EstMonteCarlo(mm,1:2)=Icum_EstMonteCarlo(:,mm_local);
    Model.AcSemSem.cum_info.EstMarkov2(mm,1)=Icum_EstMarkov2(mm_local);
    Model.AcSemSem.cum_info.EstMarkov3(mm,1)=Icum_EstMarkov3(mm_local);
    Model.AcSemSem.cum_info.EstMarkov4(mm,1)=Icum_EstMarkov4(mm_local);
    Model.AcSemSem.cum_info.EstMarkov5(mm,1)=Icum_EstMarkov5(mm_local);
end

% Floor
% fprintf('**CumInfo on Floor**\n')
% mm_local=mm_local+1;
% Model.Floor.cum_info.ExactHist5(mm,1)=Icum_ExactMem0_5(mm_local);
% Model.Floor.cum_info.EstMonteCarlo(mm,1:3)=Icum_EstMonteCarlo10_7(:,mm_local);
% Model.Floor.cum_info.EstMarkov2(mm,1)=Icum_EstMarkov2(mm_local);
% Model.Floor.cum_info.EstMarkov3(mm,1)=Icum_EstMarkov3(mm_local);
% Model.Floor.cum_info.EstMarkov4(mm,1)=Icum_EstMarkov4(mm_local);
% Model.Floor.cum_info.EstMarkov5(mm,1)=Icum_EstMarkov5(mm_local);


%Ceiling
fprintf('**CumInfo on Ceiling**\n')
mm_local=mm_local+1;
Model.Ceiling.cum_info.current_dataset.ExactHist5(mm,1)=Icum_ExactMem0_5(mm_local);
Model.Ceiling.cum_info.current_dataset.EstMonteCarlo(mm,1:2)=Icum_EstMonteCarlo(:,mm_local);
Model.Ceiling.cum_info.current_dataset.EstMarkov2(mm,1)=Icum_EstMarkov2(mm_local);
Model.Ceiling.cum_info.current_dataset.EstMarkov3(mm,1)=Icum_EstMarkov3(mm_local);
Model.Ceiling.cum_info.current_dataset.EstMarkov4(mm,1)=Icum_EstMarkov4(mm_local);
Model.Ceiling.cum_info.current_dataset.EstMarkov5(mm,1)=Icum_EstMarkov5(mm_local);

% mm_local=mm_local+1;
% Model.Ceiling.cum_info.end_dataset.ExactHist5(mm,1)=Icum_ExactMem0_5(mm_local);
% Model.Ceiling.cum_info.end_dataset.EstMonteCarlo(mm,1:2)=Icum_EstMonteCarlo(:,mm_local);
% Model.Ceiling.cum_info.end_dataset.EstMarkov2(mm,1)=Icum_EstMarkov2(mm_local);
% Model.Ceiling.cum_info.end_dataset.EstMarkov3(mm,1)=Icum_EstMarkov3(mm_local);
% Model.Ceiling.cum_info.end_dataset.EstMarkov4(mm,1)=Icum_EstMarkov4(mm_local);
% Model.Ceiling.cum_info.end_dataset.EstMarkov5(mm,1)=Icum_EstMarkov5(mm_local);


% AR
if SWITCH.AR
    fprintf('**CumInfo on AR**\n')
    mm_local=mm_local+1;
    Model.AR.cum_info.ExactHist5(mm,1)=Icum_ExactMem0_5(mm_local);
    Model.AR.cum_info.EstMonteCarlo(mm,1:2)=Icum_EstMonteCarlo(:,mm_local);
    Model.AR.cum_info.EstMarkov2(mm,1)=Icum_EstMarkov2(mm_local);
    Model.AR.cum_info.EstMarkov3(mm,1)=Icum_EstMarkov3(mm_local);
    Model.AR.cum_info.EstMarkov4(mm,1)=Icum_EstMarkov4(mm_local);
    Model.AR.cum_info.EstMarkov5(mm,1)=Icum_EstMarkov5(mm_local);
end