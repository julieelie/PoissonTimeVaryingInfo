function [stimCC] = stimulusCorrelations(Res,ParamModel)


% Arguments that are passed.
if nargin<2
    ParamModel = struct();
end
if  ~isfield(ParamModel,'MinTime') || isempty(ParamModel.MinTime)
    ParamModel.MinTime = 10; % beginning point of analysis
end
if ~isfield(ParamModel,'MaxTime') || isempty(ParamModel.MaxTime)
    ParamModel.MaxTime = 600; %end point of analysis
end
if ~isfield(ParamModel,'Increment') || isempty(ParamModel.Increment)
    ParamModel.Increment = 10; %increase the size of the spectro window with a Xms pace
end

%% Setting up working variables
% Select subset of stimuli - we are only looking at vocalizations
DataSel = selVoc(Res);

% Analysis windows.
Wins = ParamModel.MinTime:ParamModel.Increment:ParamModel.MaxTime;

% Number of time windows analysed
WinNum = length(Wins);

% Number of stims in the data set
NbStim = length(DataSel);

% Cumulative variables
rateStimAvg = zeros(1, WinNum);
stimCorr = 0.0;
ystimVar = 0.0;

nCorr = 0;

%% Loop through stims and windows

% First find average rate
for is = 1:NbStim
    rateJN = Res.JackKnife_KDE_Filtered{DataSel(is)};   % JN rate data
    indLast = min([length(Res.PSTH{DataSel(is)}), length(Res.PSTH_KDE_Filtered{DataSel(is)}), length(rateJN)]);
    
    for ww = 1:WinNum          
       tStart = Wins(ww);
       tEnd = tStart + ParamModel.Increment;
        
       if (tEnd > indLast) 
           break;
       end
       rateStimAvg(ww) = rateStimAvg(ww) + sum(Res.PSTH{DataSel(is)}(tStart:tEnd));
    end
end

rateStimAvg = rateStimAvg./NbStim;

% Calculate stimulus correlations
for is = 1:NbStim
    
    rateJN = Res.JackKnife_KDE_Filtered{DataSel(is)};   % JN rate data   
    indLast = min([length(Res.PSTH{DataSel(is)}), length(Res.PSTH_KDE_Filtered{DataSel(is)}), length(rateJN)]);
    
    for ww = 1:WinNum          
        tStart = Wins(ww);
        tEnd = tStart + ParamModel.Increment;
        
        if (tEnd > indLast) 
            break;
        end
        
        rateStim = sum(Res.PSTH{DataSel(is)}(tStart:tEnd));  
               
        if ww > 1
            tStartLast = Wins(ww-1);
            tEndLast = tStartLast + ParamModel.Increment;
            rateStimLast = sum(Res.PSTH{DataSel(is)}(tStartLast:tEndLast));
            stimCorr = stimCorr + (rateStim - rateStimAvg(ww))*(rateStimLast - rateStimAvg(ww-1));
            ystimVar = ystimVar + (rateStim - rateStimAvg(ww))^2;
            
            nCorr = nCorr + 1;
        end
    end

    fprintf(1, '\t StimCC = %.3f\n', stimCorr/ystimVar);
end

stimCC = stimCorr/ystimVar;                % Noise Correlation with JN mean rate (excluding the trial)

    
end

