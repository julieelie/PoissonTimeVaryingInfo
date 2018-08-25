function [HalfTrain1, HalfTrain2, NumTrials,SpikeTrains]=organiz_data4coherence(Trials,PSTH,ParamModel)
%% define new dataset depending on the size of the window of the model
% loop through the stims and only keep the Win first ms of them when
% they are longer than Win ms
NbStim=length(Trials);
HalfTrain1=[];
HalfTrain2=[];
NumTrials=nan(NbStim,1);
SpikeTrains = cell(1,NbStim);
for dd = 1:NbStim
    if strcmp(ParamModel.NeuroRes, 'count')
        NTrials=length(Trials{dd});
        duration=size(PSTH{dd},2);% This is the duration of the neural response in ms
        if duration<=(ParamModel.MaxWin + ParamModel.ResDelay)
            PSTH_local=zeros(NTrials,duration);
        else
            PSTH_local=zeros(NTrials,ParamModel.MaxWin + ParamModel.ResDelay);
        end

        for tt=1:length(Trials{dd})
            Spike_times=Trials{dd}{tt}(find(Trials{dd}{tt}<= (ParamModel.MaxWin + ParamModel.ResDelay)));
            ns = length(Spike_times);
            for nn=1:ns
                SpikeIdx=ceil(Spike_times(nn));
                if (SpikeIdx < 1 || SpikeIdx > duration)
                    fprintf(1, 'Warning time index out of bounds for stim# %s trial %d: time_ind = %d ntimebins = %d\n', dd, tt, SpikeIdx, duration);
                    continue;
                end
                PSTH_local(tt,SpikeIdx)=1+PSTH_local(tt,SpikeIdx);
            end
        end
        SpikeTrains{dd} = PSTH_local;
    elseif strcmp(ParamModel.NeuroRes, 'count_gaussfiltered')
        NTrials=size(Trials{dd},1);
        duration=size(PSTH{dd},2);% This is the duration of the neural response in number of points
        if duration<=(ParamModel.MaxWin + ParamModel.ResDelay)*ParamModel.Response_samprate/1000
            PSTH_local=zeros(NTrials,duration);
        else
            PSTH_local=zeros(NTrials,(ParamModel.MaxWin + ParamModel.ResDelay)*ParamModel.Response_samprate/1000);
        end


        for tt=1:size(Trials{dd},1)
            PSTH_local(tt,:) = Trials{dd}(tt,1:size(PSTH_local,2));
        end
    end
    
    if mod(NTrials,2)==0
        Trialset1=1:2:NTrials;
        Trialset2=2:2:NTrials;
        NumTrials(dd)=NTrials;
    else
        Trialset1=1:2:(NTrials-1);
        Trialset2=2:2:(NTrials-1);
        NumTrials(dd)=NTrials-1;
    end
    HalfTrain1 = [HalfTrain1 mean(PSTH_local(Trialset1,:),1)];
    HalfTrain2 = [HalfTrain2 mean(PSTH_local(Trialset2,:),1)];
    
end

end