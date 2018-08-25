function [I,P_YgivenS,H_Y, H_YgivenS]=info_poisson_model_calculus(Mu_in, varargin)
%% Calculates the information obtain from N Poisson with different means.

% This functions is used to calculate the mutual information between n
% events (stimuli, categories of stimuli, motor output...) and
% a poisson response defined by its n means. The maximum value is log2(n).
% The minum value is zero and will be obtained when the n means are identical.
% This function assumes a Poisson noise for the neural response with a mean
% between 0 and 50.
% Write your own routine for values higher than 50 using a the Normal approximation
% for the Poisson distribution.


% Requires:
% Mu_in         vector 1xn that specifies the mean (in count/bin i.e. no
%               units here) of the n Poisson distributions.

% Optional parameters
%   [...] =
%   info_model_Calculus(...,'PARAM1',VAL1,'PARAM2',VAL2,...)
%   specifies various parameters of the calculation.
% Valid parameters are the following:
%  Parameter      Value
%         
%   'MuLims'    This parameter attends the case of mu=0, which is higly
%               unlikely for neurons and lead to abberant calculations
%               because log(0) = -Inf. A better estimation of the lowest
%               spike rate should be 1/2 * 1/NTrials. By default the value
%               is set to 1/20, so works well for 10 trials datasets.

%    'Ymax'     Y_range defines the range of values that are used to
%               calculate the entropies. Ymax corresponds to the maximum
%               number of spikes you can observed per bin given the input
%               mean. We can indeed reasonably limit the calculation
%               of entropy from 0 up to the maximum number of spikes
%               expected instead of summing to + infinity.
%               Ymax is three standard deviations above the maximum
%               variance (note that the mean equals the variance for
%               Poisson distributions).


% Returns: 
% I             the mutual information in bits between n events (stimuli,
%               categories of stimuli, motor output...) and a poisson
%               response defined by its means mu_in.

% P_YgivenS     a matrix ym x n where each column gives the conditional
%               probability distribution of the response given the stimulus
%               i.e. the input mean mu_in(n). ym corresponds to the number 
%               of values of y (0:Ymax) that are used to estimate entropies.

% H_Y           a scalar that specifies the entropy of the response 

% H_YgivenS        a scalar that specifies the entropy of the response given
%               the stimulus i.e., the input means Mu_in.

%% Sanitary checks and parameters setting
% Check for max value of Mu_in to be below upper bound for a Poisson
% distribution assumption.
MAXMEANVALUE = 50;
if ~sum(Mu_in <= MAXMEANVALUE)
    ERRMSG= 'This function assumes Poisson distribution means inferior or equal to %d spikes per bin.\nWrite your own routine for values higher than 50 using a the Normal approximation for the Poisson distribution.\n';
    error(ERRMSG, MAXMEANVALUE)
end


% Sorting input arguments
Pnames = {'MuLims','Ymax'};

% Calculating default values of input arguments
MuLims=1/20;
Ymax_def = fix(3.0*sqrt(max(Mu_in)) + max(Mu_in));

% Get input arguments
Dflts  = {MuLims Ymax_def};
[MuLims, Ymax] = internal.stats.parseArgs(Pnames,Dflts,varargin{:});
if fix(Ymax)<Ymax_def
    ERRMSG= 'The range of neural response values that are used to calculate mutual information is too small.\nThe maximum value investigated is %d when it should be at least %d given the maximum mean entered as input to the function.\n';
    error(ERRMSG, Ymax, Ymax_def)
end
    


% Set the range of values that are used to calculate the entropies.
Y_range = 0:Ymax;

% To save time, calculations are only done for unique values of Mu_in    
[Mu_unique, ~, Ind_mu_unique] = unique(Mu_in);


%% Calculate the conditional entropy H(y/mu) assuming a Poisson distribution
% Initialize output variables
N_mu = length(Mu_in);
P_YgivenS = nan(length(Y_range),N_mu);
P_YgivenS_rescaled = nan(length(Y_range),N_mu);

% To save time calculation of conditional entropy is only done for unique values of mu.
H_YgivenS_unique = nan(size(Mu_unique));

% Y factorial
Fac_y = factorial(Y_range);   % Saved to save computational time.

% Now loop through unique mu and calculate the conditional entropy for each
% of them
for mm=1:length(Mu_unique)
    
    % Local value of mean
    Mu = Mu_unique(mm);
    
    % Calculate the log of the probability for numerical reasons.
    % Deal with cases where Mu is below MuLims
    if Mu < MuLims
        Mu = MuLims;
    end
    Log_P = Y_range.*log2(Mu) - log2(Fac_y) - Mu/log(2);
    
    % From logP to P
    P = 2.^(Log_P);
    % Fill in P_YgivenS for every column that shares the same mean (Mu_in)
    P_YgivenS(:,(Ind_mu_unique==mm))=repmat(P',1,sum(Ind_mu_unique==mm));
   
    % Check that sum(P) sum to 1 and rescale if not to have a more
    % realistic distribution.
    P = P ./sum(P);
    % Fill in P_YgivenS_rescaled for every column that shares the same mean (Mu_in)
    P_YgivenS_rescaled(:,(Ind_mu_unique==mm))=repmat(P',1,sum(Ind_mu_unique==mm));
    % conditional entropy value for each unique mu (each stim) across all y of the world
    H_YgivenS_unique(mm) = sum(-P.*log2(P + (P==0))); % Note that adding P==0 in the logarithmic expression ensures that P*log2(P) = 0 when P=0.
end

% Given that all stimuli have equal probability: H(y|s) = expectation(H(y|si)) = sum(p(si)*H(y|si)) = sum(H(y|si))/(number of stims)
H_YgivenS = sum(H_YgivenS_unique(Ind_mu_unique))/N_mu;


%% Calculate the entropy of the response y based on Poisson distributions
% Probability of responses (marginal)
P_Y = sum(P_YgivenS_rescaled,2)./N_mu;

% Entropy of the response
H_Y = sum(-P_Y.*log2(P_Y + (P_Y==0)));% entropy of the response = entropy of the marginals = entropy of the average response probabilities over all Mu_in

%% Calculate mutual information
I = (H_Y - H_YgivenS);

end