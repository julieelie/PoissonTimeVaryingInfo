function [Icum, HY, HYgivenS, Icum_corr_mean, Icum_corr_std, N_MC_tot]=cumulative_info_poisson_model_calculus_MCJK(P_YgivenS_Full,P_YgivenS_JK, NTrials, varargin)
%% Calculates the cumulative information in a neural response given conditional probability distribution of responses given stims
% at each time point.  Although it is assumed that this stimulus conditional probability is time
% independent (i.e know that stim identity and the rate at a previous time does not tell you 
% about the rate at the current time - as in a poisson distribution), the joint probability
% distribution is not making the cumulative information calculation
% difficult.
    % To calculate the cumulative information at a given time bin T, one
    % needs to estimate the distribution of joint probabilities of reponse
    % sequences. This distribution can be very very large as the number of
    % time bins in the neural response sequence increases. This algorithm
    % uses a Monte Carlo approximation to estimate
    % this joint distribution and calculate the entropy of the response and the
    % conditional entropy (entropy of the response given the event
    % (stimulus)), and finally the cumulative mutual information. The
    % number of samples used in the Monte Carlo is upper bounded by
    % MaxMCParameter.
    % The algorithm expects as input two cell arrays (1 x T) of T matrices
    % (c x N) containing the conditional probability of observing a reponse
    % c given the stimulus N at time t.
    % The matrices from the first cell array correspond to the conditional
    % probabilities obtained using a spike rate estimated with all the
    % trials (NTrials stimulus presentations). The matrices from the second cell
    % array correspond to the conditional probabilities obtained using a
    % spike rate estimated through a Jack-knife procedure (e.i. with all-1
    % the trials or stimulus presentations). This second set of matrices is
    % used to estimate the error on the calculation of cumulative
    % information by a Jack-knife procedure.
    % To see how to obtain matrices of conditional probabilities at each
    % time bin, please refer to info_model_Calculus.m
    
    
% Requires:
% P_YgivenS_Full a cell array of length equals to the number of time points
%               T in the neural responses. Each cell contains a c*N matrix
%               of probabilities,the conditional probability of a given
%               neural spike count c at the particular time point t for the
%               stimulus N. These probabilities are obtained by running
%               info_model_Calculus.m on the spike rate estimated using all
%               stimulus presentations NTrials.
%
% P_YgivenS_JK  a cell array with # of rows equals to the number of Jack-knife
%               estimations of the conditional probability distributions
%               and the # of columns equals the number of time points T.
%               Each cell of P_YgivenS_JK contains a c*N matrix of
%               probabilities,the conditional probability of a given
%               neural spike count c at the particular time point t for the
%               stimulus N. These probabilities are obtained by running
%               info_model_Calculus.m on the spike rate estimated using all
%               but one stimulus presentations NTrials, following a
%               Jack-knife procedure.
%
% NTrials       Number of stimulus presentations (trials) used to estimate
%               the conditional response probabilities reported in
%               P_YgivenS_Full. NTrials can be the number of rows
%               in P_YgivenS_JK depending on your setting.



% Optional parameters
%   [...] =
%   cumulative_info_poisson_model_calculus_MCJK(...,'PARAM1',VAL1,'PARAM2',VAL2,...)
%   specifies various parameters of the calculation.
% Valid parameters are the following:
%  Parameter      Value
%         
%   'MaxMCParameter' is the maximum number of samples in the Monte Carlo estimation
%                    of the entropy of the response defaulted to 5*10^6
%   'IncrMCParameter' is the increment by which the number of samples is
%                   increased in the MonteCarlo estimation of the entropies
%   'ConvThresh'  is the targeted maximum error in bits on the calculation of
%                   cumulative information at each bin
%   'DebugFig'      set to 1 to see some debugging figures
%   'ScaleY'        Set to 1 to scale the distribution of conditional
%                   probabilities given the stimulus the same way for both
%                   entropy calculations
%   'Verbose'       set to 1 to see a lot of output, 0 for none default is
%                   none



% Returns
% Icum          The cumulative mutual information in bits between N events 
%               (stimuli, categories of stimuli, motor output...) and an
%               inhomogenous poisson response of length T bins defined by
%               the probability distributions of neural response sequences
%               P_YgivenS_Full. Icum_Full is a 1x1 scalar.
%
% H_y           a scalar that specifies the entropy of the neural response 
%
% H_y_mu        a scalar that specifies the entropy of the neural response
%               given the stimulus.
%
% Icum_corr_mean  The Jack-knife biais corrected cumulative mutual
%               information in bits between N events (stimuli, categories
%               of stimuli, motor output...) and an inhomogenous poisson
%               response of length T bins defined by the probability
%               distributions of neural response sequences P_YgivenS_Full
%               and P_YgivenS_JK. Icum_corr_mean is a 1x1 scalar.
%
% Icum_corr_std  The error on Icum_cor_mean as estimated by the Jack-knife
%               procedure. Icum_corr_sd is a 1x1 scalar
%
% N_MC_tot     The number of samples used in the Monte Carlo.

tic
%% Sorting input arguments
pnames = {'MaxMCParameter','IncrMCParameter','ConvThresh','ScaleY','DebugFig', 'Verbose'};

% Get input arguments
dflts  = {5*10^6 10^5 0.2 1 0 0};
[N_MC_max, N_MC, ConvThresh, ScaleY, DebugFig, Verbose] = internal.stats.parseArgs(pnames,dflts,varargin{:});

%% Set some parameters
win = length(P_YgivenS_Full);
Nb_Boot = length(P_YgivenS_JK);
Nb_JKSet = size(P_YgivenS_JK{1},1);

%% Construct the data set of individual p for all JK Bootstraps and also for the full dataset
if ScaleY  % Scale probability distributions so that they sum to 1 for each stimulus
    for ww=1:(win)
        P_YgivenS_Full{ww}=P_YgivenS_Full{ww} ./ repmat(sum(P_YgivenS_Full{ww},1),size(P_YgivenS_Full{ww},1),1);
    end
    
    parfor bb=1:(Nb_Boot)
        for jk=1:(Nb_JKSet)
            for ww=1:(win)
                P_YgivenS_JK{bb}{jk,ww} = P_YgivenS_JK{bb}{jk,ww} ./ repmat(sum(P_YgivenS_JK{bb}{jk,ww},1),size(P_YgivenS_JK{bb}{jk,ww},1),1);
            end
        end
    end
end


%% Calculating Cumulative Mutual Information
% The startegy used here is to progressively increase the number fo samples
% used in the Monte Carlo to estimate the information until the error on
% the calculation drops below a reasonnable threshold set by ConvThresh or
% the maximum number of samples set by MaxMCParameter is reached.

% Initialize output variables
N_MC_tot = 0;
Icum_corr_std = ConvThresh+1; % set the error on information to an arbitrary value above the threshold

% cdf of the marginal probabilities (probabilities of
% responses at each time point) using the probability distributions
% estimated with all the stimulus presentations/trials (full dataset).
P_Yt=cell(win,1);
for ww=1:win
    P_Yt{ww} = cumsum(mean(P_YgivenS_Full{ww},2)./sum(mean(P_YgivenS_Full{ww},2)));
end

% Begin the calculations
while Icum_corr_std>ConvThresh && N_MC_tot<N_MC_max
    % Calculating the entropies of the response following a Monte Carlo
    % estimation, using the number of samples requested N_MC. The function
    % is defined below.
    [HY_Full_local,HYgivenS_Full_local,HY_JK_local,HYgivenS_JK_local]= cuminfo_MC(P_YgivenS_Full,P_YgivenS_JK, P_Yt, N_MC);
    
    % update the new total number of samples used to estimate entropies
    N_MC_tot_new = N_MC + N_MC_tot;
    
    % Update values of entropies with the values obtained for the last
    % batch of samples
    if N_MC_tot == 0 % initialization of entropy values
        HY = HY_Full_local;
        HYgivenS = HYgivenS_Full_local;
        HYgivenS_JK = HYgivenS_JK_local;
        HY_JK = HY_JK_local;
    else
        HY = (HY * N_MC_tot + HY_Full_local * N_MC) / N_MC_tot_new;
        HYgivenS = (HYgivenS * N_MC_tot + HYgivenS_Full_local * N_MC) / N_MC_tot_new;
        HYgivenS_JK = (HYgivenS_JK .* N_MC_tot + HYgivenS_JK_local .* N_MC) ./ N_MC_tot_new;
        HY_JK = (HY_JK .* N_MC_tot + HY_JK_local.*N_MC) ./ N_MC_tot_new;
    end
    
    % update the past total number of samples used in the Monte Carlo
    % approximation
    N_MC_tot = N_MC_tot_new;
    
    % Calculate the information given the current total number of samples
    Icum = HY-HYgivenS;
    Icum_JK = HY_JK-HYgivenS_JK;
    
    % Calculate the Jack-knife biais corrected value of information and its error 
    Icum_JKcorrected = NTrials * repmat(Icum,Nb_Boot,Nb_JKSet) - (NTrials-1) * Icum_JK;
    Icum_corr_std = (mean(var(Icum_JKcorrected,0,2)))^0.5;
    Icum_corr_mean = mean(mean(Icum_JKcorrected));
end


if Verbose
    fprintf(1,'Cumulative information = %f bits\n',Icum);
    ElapsedTime = toc;
    fprintf(1,'cumulative_info_poisson_model_calculus_MCJK run for %d seconds with %d samples\n',ElapsedTime, N_MC_tot);
end




%% Strategy: Monte Carlo Estimate of the entropy of the response
function[HY_Full,HYgivenS_Full,HY_JK,HYgivenS_JK]= cuminfo_MC(P_YgivenS_Full,P_YgivenS_JK, P_Yt, N_MC)
% Calculates the entropies of neural response sequences
% modeled as inhomogenous Poisson functions as defined by their
% probability distributions P_YgivenS_all, using a Monte Carlo
% approximation with N_MC samples.

if Verbose
    fprintf(1, 'Calculation of Exact entropy (entropy of the neural response)\nusing the Monte Carlo estimation with %d samples\n',N_MC);
end

%% Set the values of MonteCarlo samples used for all datasets
MC_Randu = rand(N_MC, win);
Resp_MC = nan(N_MC, win);


%% First calculate entropies for the full dataset as we need QY_MC of that
% dataset to calculate the entropies of JK datasets

% Find the number of stimuli and initialize output variables
NStim_local = size(P_YgivenS_Full{end},2);
PYgivenS_MC = nan(N_MC,NStim_local); %exact joint probability of MC sequences of responses for each stimulus
PY_MC = nan(N_MC,1); %exact joint probability of MC sequences of responses averaged over the stimuli
QY_MC = nan(N_MC,1); %estimation of the joint probability of MC sequences of responses using the marginals

% Loop through the number of samples
parfor ss=1:N_MC
    if sum(ss==(1:10)*N_MC/10) && Verbose
        fprintf('%d/%d MC samples Full dataset\n',ss,N_MC);
    end
    
    % Probabilistically determine a sequence of responses using the marginal
    % probabilities and extract the conditional probabilities corresponding
    % to that response sequence for all the stimuli:
    % For each time point, take a random number from the uniform
    % distribution 0-1 and find the neural response that corresponds to
    % that probability in the full dataset. Then retrieve the conditional
    % probability of that neural response given each stimuli.
    P_YgivenS_local_Resp = nan(win,NStim_local);
    for www=1:win
        u=MC_Randu(ss,www);
        Resp_MC(ss, www) = sum(P_Yt{www}<u)+1;
        P_YgivenS_local_Resp(www,:) = P_YgivenS_Full{www}(Resp_MC(ss,www),:);
    end

    % Calculate the exact joint probability of that sequence of responses
    % and store it.
    % PYgivenS_MC is the conditional joint probabilities of the neural
    % response sequence defined by the Monte Carlo sample for all
    % stimuli.
    %
    % PY is the actual joint probability (the unconditional probability)
    % calculated by averaging the conditional joint probabilities over stimuli.
    % 
    % QY is the joint probability of the neural response sequence in the 
    % proposal distribution used to find Monte Carlo estimates - it assumes
    % that the joint probability disribution is independent across time and
    % can be obtained from the product at each time bin.
    PYgivenS_MC(ss,:) = prod(P_YgivenS_local_Resp,1);
    PY_MC(ss) = mean(PYgivenS_MC(ss,:));           
    QY_MC(ss) = prod(mean(P_YgivenS_local_Resp,2));

end

% Calculate the MC estimate of the conditional response entropy to the
% stimulus
HYgivenS_Full = - sum(sum(PYgivenS_MC./repmat(QY_MC,1,NStim_local).*log2(PYgivenS_MC),1))/(N_MC*NStim_local);
if Verbose
    fprintf(1, 'Conditional entropy (entropy of the neural response given the stimulus): %f\n', HYgivenS_Full);
end

% Calculate the MC estimate of the response entropy weight corrected
HY_Full = sum(-PY_MC./QY_MC.*log2(PY_MC))/N_MC;

if DebugFig
    figure()
    plot(1:100, QY_MC(1:100), 1:100,PY_MC(1:100))
    legend('Product of probas','exact joint probas')
    xlabel('Monte Carlo samples')
    ylabel('probability')
end

%% Now calculate for JackKnife sets using the distribution of probabilities calculated for the full dataset to choose the MC samples

% Initialize output variables
HYgivenS_JK = nan(Nb_Boot, Nb_JKSet);
HY_JK = nan(Nb_Boot,Nb_JKSet);

% Loop through bootstrap versions of Jack-Knife estimation of distributions
% of conditional response probablities 
parfor bout=1:(Nb_Boot)
    for jkk=1:Nb_JKSet
        % initialize variables
        PYgivenS_MC_bb = nan(N_MC,NStim_local); %exact joint probability of MC sequences of responses for each stimulus
        PY_MC_bb = nan(N_MC,1); %exact joint probability of MC sequences of responses averaged over the stimuli
    
        % Loop through the number of samples
        for ss=1:N_MC
            if sum(ss==(1:10)*N_MC/10) && Verbose
                fprintf('%d/%d MC samples JK bootstrap %d/%d\n',ss,N_MC,bout,Nb_Boot);
            end
        
            % The sequence of responses previously determine for that MC sample
            % ss is Resp_MC(ss, :). Here, extract, for all the stimuli,the 
            % JK conditional probabilities corresponding to that response sequence:
            P_YgivenS_local_Resp = nan(win,NStim_local);
            for www=1:win
                P_YgivenS_local_Resp(www,:) = P_YgivenS_JK{bout}{jkk,www}(Resp_MC(ss, www),:);
            end
        
        
            % Calculate the exact joint probability of that sequence of responses
            % and store it
            % Conditional joint probability of that neural sequence for each
            % stimulus
            PYgivenS_MC_bb(ss,:) = prod(P_YgivenS_local_Resp,1);
        
            % actual joint probability (the unconditional probability)
            PY_MC_bb(ss) = mean(PYgivenS_MC_bb(ss,:));
        
        end


        % Calculate the MC estimate of the JK conditional response entropy to the
        % stimulus
        HYgivenS_JK(bout,jkk) = - sum(sum(PYgivenS_MC_bb./repmat(QY_MC,1,NStim_local).*log2(PYgivenS_MC_bb),1))/(N_MC*NStim_local);
        if Verbose
            fprintf(1, 'JKSet %d Bootstrap %d: Conditional entropy (entropy of the neural response given the stimulus): %f\n',jkk,bout, HYgivenS_JK(bout,jkk));
        end

        % Calculate the JK MC estimate of the response entropy weight corrected
        HY_JK(bout,jkk) = sum(-PY_MC_bb./QY_MC.*log2(PY_MC_bb))/N_MC;

        if Verbose
            fprintf(1, 'JKSet %d Bootstrap %d: Exact entropy (entropy of the neural response)\nusing the Monte Carlo estimation with %d samples: %f\n',jkk,bout,N_MC, HY_JK(bout,jkk));
        end
    end
end

end

end
