function [Icum_EstMonteCarlo6,Icum_EstMonteCarlo4,Icum_EstMonteCarlo3,Icum_EstMonteCarlo2,Icum_ExactMem0_5,Icum_EstMarkov2,Icum_EstMarkov3,Icum_EstMarkov4,Icum_EstMarkov5] = cumInfo_4AN(P_YgivenS)

Nb_Win = length(P_YgivenS);
Icum_EstMonteCarlo6_local = nan(2,Nb_Win);
Icum_EstMonteCarlo4_local = Icum_EstMonteCarlo6_local;
Icum_EstMonteCarlo3_local = Icum_EstMonteCarlo6_local;
Icum_EstMonteCarlo2_local = Icum_EstMonteCarlo6_local;
Icum_ExactMem0_5 = nan(1,Nb_Win);
Icum_EstMarkov2 = Icum_ExactMem0_5;
Icum_EstMarkov3 = Icum_ExactMem0_5;
Icum_EstMarkov4 = Icum_ExactMem0_5;
Icum_EstMarkov5 = Icum_ExactMem0_5;
HY_Markov2 = nan(1,Nb_Win);
HY_Markov3 = nan(1,Nb_Win);
HY_Markov4 = nan(1,Nb_Win);
HY_Markov5 = nan(1,Nb_Win);


for tt=2:Nb_Win
    tstart = tic;
    fprintf('Time point %d/%d\n', tt, Nb_Win);
    P_YgivenS_local = P_YgivenS(1:tt);
    % Monte Carlo estimation with full memory
    [Icum_EstMonteCarlo6_local(:,tt),~]=info_cumulative_model_Calculus(P_YgivenS_local,'CalMode','MonteCarlo', 'MCParameter',10^6);
    % Monte Carlo estimation with full memory
    [Icum_EstMonteCarlo4_local(:,tt),~]=info_cumulative_model_Calculus(P_YgivenS_local,'CalMode','MonteCarlo', 'MCParameter',10^4);
    % Monte Carlo estimation with full memory
    [Icum_EstMonteCarlo3_local(:,tt),~]=info_cumulative_model_Calculus(P_YgivenS_local,'CalMode','MonteCarlo', 'MCParameter',10^3);
    % Monte Carlo estimation with full memory
    [Icum_EstMonteCarlo2_local(:,tt),~]=info_cumulative_model_Calculus(P_YgivenS_local,'CalMode','MonteCarlo', 'MCParameter',10^2);
    % Exact calculation with 50 ms memory
    [Icum_ExactMem0_5(tt),~]=info_cumulative_model_Calculus(P_YgivenS_local,'CalMode','Exact_Mem', 'Exact_history',5);
    if tt==2
        % Markov chain estimation 20 ms
        [Icum_EstMarkov2(tt),HY_Markov2(tt),~]=info_cumulative_model_Calculus(P_YgivenS_local,'CalMode','MarkovChain', 'MarkovParameters',[2,1]);
        % Markov chain estimation 30 ms
        [Icum_EstMarkov3(tt),HY_Markov3(tt),~]=info_cumulative_model_Calculus(P_YgivenS_local,'CalMode','MarkovChain', 'MarkovParameters',[3,1]);
        % Markov chain estimation 40 ms
        [Icum_EstMarkov4(tt),HY_Markov4(tt),~]=info_cumulative_model_Calculus(P_YgivenS_local,'CalMode','MarkovChain', 'MarkovParameters',[4,1]);
        % Markov chain estimation 50 ms
        [Icum_EstMarkov5(tt),HY_Markov5(tt),~]=info_cumulative_model_Calculus(P_YgivenS_local,'CalMode','MarkovChain', 'MarkovParameters',[5,1]);
    else
        % Markov chain estimation 20 ms
        [Icum_EstMarkov2(tt),HY_Markov2(tt),~]=info_cumulative_model_Calculus(P_YgivenS_local,'CalMode','MarkovChain', 'MarkovParameters',[2,1], 'HY_old', HY_Markov2(tt-1));
        % Markov chain estimation 30 ms
        [Icum_EstMarkov3(tt),HY_Markov3(tt),~]=info_cumulative_model_Calculus(P_YgivenS_local,'CalMode','MarkovChain', 'MarkovParameters',[3,1], 'HY_old', HY_Markov3(tt-1));
        % Markov chain estimation 40 ms
        [Icum_EstMarkov4(tt),HY_Markov4(tt),~]=info_cumulative_model_Calculus(P_YgivenS_local,'CalMode','MarkovChain', 'MarkovParameters',[4,1], 'HY_old', HY_Markov4(tt-1));
        % Markov chain estimation 50 ms
        [Icum_EstMarkov5(tt),HY_Markov5(tt),~]=info_cumulative_model_Calculus(P_YgivenS_local,'CalMode','MarkovChain', 'MarkovParameters',[5,1], 'HY_old', HY_Markov5(tt-1));
    end
    telapsed = toc(tstart);
    fprintf('Elapsed time: %d s\n', telapsed)
end
Icum_EstMonteCarlo6 = Icum_EstMonteCarlo6_local(1,:);
Icum_EstMonteCarlo4 = Icum_EstMonteCarlo4_local(1,:);
Icum_EstMonteCarlo3 = Icum_EstMonteCarlo3_local(1,:);
Icum_EstMonteCarlo2 = Icum_EstMonteCarlo2_local(1,:);
end