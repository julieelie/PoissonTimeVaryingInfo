function [Poisson_PL] = PoissonTestNeuralResponses(Res,ParamModel)

% Global argurments
figFlg=1;     % Flag for plotting
nBoot = 1000; % Number of boot points

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
Wins = ParamModel.MinWin:ParamModel.Increment:ParamModel.MaxWin;

% Number of time windows analysed
WinNum = length(Wins);

% Number of stims in the data set
NbStim = length(DataSel);

% Initialize output variables
Poisson_PL = nan(NbStim, WinNum);


%% Loop through stims and windows

for is = 1:NbStim
    trial = Res.Trials{DataSel(is)};                  % Raw trial data
    rateJN = Res.JackKnife_KDE_Filtered{DataSel(is)};   % JN rate data
    
    for ww = 1:WinNum  
        Win = Wins(ww); 
        tStart = Win;
        tEnd = Win + ParamModel.Increment;
        
        LL = 0;   % Log likelihood of data
        LLBoot = zeros(1,nBoot); % Log Likelihood bootstrapped
        rateBoot = sum(Res.PSTH_KDE_Filtered{DataSel(is)}(tStart:tEnd));
        
        for it = 1:length(trial)
            
            % Number of spikes observed
            y=sum((trial{it}>=tStart).*(trial{it}<tEnd));
            
            % rate from JN
            rateExpected = sum(rateJN(it, tStart:tEnd));
            
            % LL
            LL = LL + log(poisspdf(y, rateExpected));
            
            % Repeat for bootrap data
            for ib = 1:nBoot
                yboot = poissrnd(rateBoot);
                LLBoot(ib) = LLBoot(ib) +  log(poisspdf(yboot, rateBoot));
            end
        end
        nDiscovery = sum(LLBoot > LL);
        Poisson_PL(is, ww) = 1.0 - nDiscovery/nBoot;
        if figFlg > 1
            figure(1);
            histogram(LLBoot);
            l = axis;
            hold on;
            plot([LL LL], [l(3) l(4)], 'k--');
            title(sprintf('Spike rate %f', rateBoot*1000.0/ParamModel.Increment));
            pause();
            hold off;
        end
            
    end
    fprintf(1, 'Done with Stim %d/%d - %d/%d Windows Reject Poisson at 5%%\n', is, NbStim, sum(Poisson_PL(is,:)<= 0.05), WinNum);
end
%% Plot if asked 
if figFlg
    figure();
    
    rateUsed = nan(NbStim, WinNum);
    for is = 1:NbStim
        for ww = 1:WinNum  
           Win = Wins(ww); 
           tStart = Win;
           tEnd = Win + ParamModel.Increment;
       
           rateUsed(is, ww) = sum(Res.PSTH_KDE_Filtered{DataSel(is)}(tStart:tEnd));
        end
    end
       

    plot(reshape(Poisson_PL, 1,  NbStim*WinNum), reshape(rateUsed, 1, NbStim*WinNum), '+');
    xlabel('Prob is Poisson');
    ylabel('Rate in 10 ms');
    
    figure();
    histogram(reshape(Poisson_PL, 1,  NbStim*WinNum));
    
end
    

end

