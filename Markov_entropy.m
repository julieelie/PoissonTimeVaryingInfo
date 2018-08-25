function [entropy,Mat_Num_stim] = Markov_entropy(P_YgivenS_AllWin, MarkovHistory,MinProbThresh, MinProb,entropy_old, Verbose)

if nargin<2
    MarkovHistory=[1 1];% first element is the number of past events used, the second is the sampling, second has to be smaller than the first(e.g. 6 and 2 would be looking at t-2 t-4 and t-6)
end
if nargin<3
    MinProbThresh=1;
end
if nargin<4
    MinProbThresh=1;
    MinProb = 1/(8*365*24*60*60*5); %1/number of 20ms bins in a life of a zebra finch in the lab
end
if nargin < 5
    entropy_old = [];
end
if nargin < 6
    Verbose=0; % set to 1 to see more outputs
end

if size(P_YgivenS_AllWin,1)==1
    %% initilize entropy
    P_YgivenS_local = P_YgivenS_AllWin{1};
    P_y_i = sum(P_YgivenS_local,2)./size(P_YgivenS_local,2);
    % Check that sum(P) sum to 1 and rescale if not to have a more
    % realistic distribution.
    if sum(P_y_i)~=1
        P_y_i = P_y_i ./sum(P_y_i);
        if Verbose
            fprintf('rescaled the distribution of y for Hy0 calculation in Markov_entropy\n');
        end
    end
    
    % Calculate log probability
    entropy = sum(-P_y_i.*log2(P_y_i + (P_y_i==0)));% entropy of the response for each unique y of the world across all mu of the model (across all stimuli)
    Mat_Num_stim={};
else
    %% get the sum of the entropies conditional to previous time points...
    ...(e.g. for MarkovHistory=1 H(t-1)/H(t-2), H(t-2)/H(t-3)) of the...
        ...previous time points
    if isempty(entropy_old)||isnan(entropy_old)
        [entropy_old,Mat_Num_stim_old] = Markov_entropy(P_YgivenS_AllWin(1:end-1),MarkovHistory,MinProbThresh, MinProb,[], Verbose);
    end
    
    if MarkovHistory(1)>0
        %% get the entropy conditional to previous time points for that last time point
        % Construct the matrix of conditional probabilities of the response
        % across stimuli. This matrix has MarkovHistory+1 dimensions and
        % represents p(y(t)/y(t-1),y(t-2)...y(t-MarkovHistory))
        pastpointmax = min(MarkovHistory(1), length(P_YgivenS_AllWin)-1);%Oldest time point that we can investigate
        %stepmax = min(MarkovHistory(2), length(P_YgivenS_AllWin)-1);% we want to investigate at least the immediate last time point if no other point further back in time is available
        pastpoints = flip([0 1:MarkovHistory(2):pastpointmax]);%0 is for pyt
        
        if MarkovHistory(1)>=length(P_YgivenS_AllWin)-1 && MarkovHistory(2)==1 && exist('Mat_Num_stim_old', 'var') % Use the previous calculations of the product of probabilities Mat_Num_stim_old since the oldest time point is still the same as the previous entropy calculation
            [Mat_Num_stim,Mat_Den_stim] = insideMultiMat_Dim(P_YgivenS_AllWin(end-pastpoints),Mat_Num_stim_old,MinProbThresh, MinProb);
        else
            [Mat_Num_stim,Mat_Den_stim] = insideMultiMat_Dim(P_YgivenS_AllWin(end-pastpoints),MinProbThresh, MinProb);
        end
        
        % Average probabilities over stimuli
        Mat_PresentPastPoints = squeeze(sum(Mat_Num_stim{1},1))./size(Mat_Num_stim{1},1);
        Mat_PastPoints = squeeze(sum(Mat_Den_stim{1},1))./size(Mat_Num_stim{1},1);
        if size(Mat_PastPoints,1)==1 && ismatrix(Mat_PastPoints)
            Mat_PastPoints=Mat_PastPoints';
        end
        clear Mat_Den_stim
        %clear Mat_Num_stim
        
        % Conditional probability based on MarkovHistory previous measurements
        PY_Cond = bsxfun(@rdivide, Mat_PresentPastPoints, Mat_PastPoints);
        
        % Take the p*log2(p) of the conditional proba
        PY_Cond = PY_Cond.*log2(PY_Cond + (PY_Cond==0));
        PY_Cond = sum(PY_Cond,ndims(PY_Cond)); % here we sum over probabilities over yt which are along the last dimension
        if sum(size(PY_Cond)>1) ~= (length(pastpoints)-1)
            fprintf('Something wrong in Markov_entropy: the dimension of PY_cond is %d when it should be %d the same as or smaller than the markov history of %d with a step of%d\n', ndims(PY_Cond),length(pastpoints), MarkovHistory(1),MarkovHistory(2))
        end
        
        if sum(size(PY_Cond) ~= size(Mat_PastPoints))
            fprintf('in Markov_entropy PY_Cond and Mat_PastPoints should have the same sizes and dimensions\nTheir dimensions %d (PY_Cond) and %d (Mat_PastPoints) should be that of Markov history %d\nTheir sizes should be that of yt-MarkovHistory... yt-1\n', ndims(PY_Cond), ndims(Mat_PastPoints), MarkovHistory);
        end
        
        % Multiply the entropy of the conditional probability by the joint
        % probability of past events and sum up all values to get the entropy
        entropy_local= - sum(reshape(PY_Cond.*Mat_PastPoints, numel(PY_Cond),1));
        
        % Sum up local entropy with previous steps entropy
        entropy = entropy_local + entropy_old;
    
    elseif MarkovHistory(1)==0
        % No history for the markov chain, the entropy is just the sum of the
        % entropy of the response at each time point, as if there was no
        % dependencies between time points
        %% initilize entropy
        P_YgivenS_local = P_YgivenS_AllWin{end};
        P_y_i = sum(P_YgivenS_local,2)./size(P_YgivenS_local,2);
        % Check that sum(P) sum to 1 and rescale if not to have a more
        % realistic distribution.
        if sum(P_y_i)~=1
            P_y_i = P_y_i ./sum(P_y_i);
            if Verbose
                fprintf('rescaled the distribution of y for Hy calculation in Markov_entropy\n');
            end
        end
        
        % Calculate log probability
        entropy_local = sum(-P_y_i.*log2(P_y_i + (P_y_i==0)));% entropy of the response for each unique y of the world across all mu of the model (across all stimuli)
        Mat_Num_stim={};
        % Sum up local entropy with previous steps entropy
        entropy = entropy_local + entropy_old;

    else
        fprintf(1,'in Markov_entropy: The markov history value is not valid it needs to be a positive integer');
        return 
    end
end