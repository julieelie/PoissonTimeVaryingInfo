function [Icum_EstMonteCarlo6,Icum_EstMonteCarlo4,Icum_EstMonteCarlo3,Icum_EstMonteCarlo2,Icum_EstMonteCarlo6_JK, Icum_EstMonteCarlo4_JK, Icum_EstMonteCarlo3_JK, Icum_EstMonteCarlo2_JK] = cumInfo_MC_4AN(P_YgivenS, P_YgivenS_JK)

Nb_Win = length(P_YgivenS);
Icum_EstMonteCarlo6 = nan(1,Nb_Win);
Icum_EstMonteCarlo4 = Icum_EstMonteCarlo6;
Icum_EstMonteCarlo3 = Icum_EstMonteCarlo6;
Icum_EstMonteCarlo2 = Icum_EstMonteCarlo6;

Icum_EstMonteCarlo6_JK = nan(size(P_YgivenS_JK,1),Nb_Win);
Icum_EstMonteCarlo4_JK = Icum_EstMonteCarlo6_JK;
Icum_EstMonteCarlo3_JK = Icum_EstMonteCarlo6_JK;
Icum_EstMonteCarlo2_JK = Icum_EstMonteCarlo6_JK;


for tt=2:Nb_Win
    tstart = tic;
    fprintf('Time point %d/%d\n', tt, Nb_Win);
    P_YgivenS_local = P_YgivenS(1:tt);
    P_YgivenS_JK_local = P_YgivenS_JK(:,1:tt);
    % Monte Carlo estimation with full memory
    [Icum_EstMonteCarlo6(tt),~, ~, Icum_EstMonteCarlo6_JK(:,tt),~,~]=info_cumulative_model_Calculus_MCJK(P_YgivenS_local, P_YgivenS_JK_local, 'MCParameter',10^6);
    % Monte Carlo estimation with full memory
    [Icum_EstMonteCarlo4(tt),~, ~, Icum_EstMonteCarlo4_JK(:,tt),~,~]=info_cumulative_model_Calculus_MCJK(P_YgivenS_local,P_YgivenS_JK_local, 'MCParameter',10^4);
    % Monte Carlo estimation with full memory
    [Icum_EstMonteCarlo3(tt),~, ~, Icum_EstMonteCarlo3_JK(:,tt),~,~]=info_cumulative_model_Calculus_MCJK(P_YgivenS_local,P_YgivenS_JK_local, 'MCParameter',10^3);
    % Monte Carlo estimation with full memory
    [Icum_EstMonteCarlo2(tt),~, ~, Icum_EstMonteCarlo2_JK(:,tt),~,~]=info_cumulative_model_Calculus_MCJK(P_YgivenS_local,P_YgivenS_JK_local, 'MCParameter',10^2);
    telapsed = toc(tstart);
    fprintf('Elapsed time: %d s\n', telapsed)
end
end