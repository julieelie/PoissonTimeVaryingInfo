function [HalfTrain1, HalfTrain2, Min_Duration]=organiz_SpikeArray4coherence(Spike_array,ParamModel)
%% define new dataset depending on the size of the Max window
% loop through the stims and only keep the Win first ms of the response when
% they are longer than MaxWin ms


NbStim=length(Spike_array);
Duration = nan(length(Spike_array),1);

MaxWin_local = ParamModel.MaxWin* ParamModel.Response_samprate/1000;
HalfTrain1 = [];
HalfTrain2 = [];


% Check if all responses have the same duration and take the smallest one
% if not
for dd = 1:NbStim
    Duration(dd) = size(Spike_array{dd},2);
end
Min_Duration = min(Duration);
if MaxWin_local>Min_Duration
    fprintf('Coherence will be calculated on responses of %dms to respect the duration of the shortest reponse\n', Min_Duration*1000/ParamModel.Response_samprate);
    MaxWin_local = Min_Duration;
end

for dd = 1:NbStim
    NTrials = size(Spike_array{dd},1);
    
    if mod(NTrials,2)==0
        Trialset1=1:2:NTrials;
        Trialset2=2:2:NTrials;
    else
        Trialset1=1:2:(NTrials-1);
        Trialset2=2:2:(NTrials-1);
    end
%    fprintf('Stim %d\n',dd)
    HalfTrain1 = [HalfTrain1 mean(Spike_array{dd}(Trialset1,1:MaxWin_local),1)];
    HalfTrain2 = [HalfTrain2 mean(Spike_array{dd}(Trialset2,1:MaxWin_local),1)];
    
end
end
