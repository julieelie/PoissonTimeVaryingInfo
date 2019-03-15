function [Poisson_P, noiseCC, noiseCC_P] = PoissonTestNeuralResponses(Res,ParamModel)

% Global argurments
figFlg=0;     % Flag for plotting - 1 summary plots, 2 debugging plots + summary plots
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
Wins = ParamModel.MinTime:ParamModel.Increment:ParamModel.MaxTime;

% Number of time windows analysed
WinNum = length(Wins);

% Number of stims in the data set
NbStim = length(DataSel);

% Initialize output variables
Poisson_P = nan(NbStim, WinNum);
noiseCC = nan;
noiseCC_P = nan;

% Cumulative variables
noiseCorr = 0.0;
yVar = 0.0;
noiseCorrBS = zeros(1, nBoot);
yVarBS = zeros(1, nBoot);
nCorr = 0;


%% Loop through stims and windows

for is = 1:NbStim
    trial = Res.Trials{DataSel(is)};                  % Raw trial data
    rateJN = Res.JackKnife_KDE_Filtered{DataSel(is)};   % JN rate data
    yBootLast = zeros(length(trial), nBoot);            % to store the values of the last random data
    
    for ww = 1:WinNum          
        tStart = Wins(ww);
        tEnd = tStart + ParamModel.Increment;

        LL = 0;   % Log likelihood of data
        LLBoot = zeros(1,nBoot); % Log Likelihood bootstrapped
        rateBoot = sum(Res.PSTH_KDE_Filtered{DataSel(is)}(tStart:tEnd));
        
        if ww > 1
            tStartLast = Wins(ww-1);
            tEndLast = tStartLast + ParamModel.Increment;
            rateBootLast = sum(Res.PSTH_KDE_Filtered{DataSel(is)}(tStartLast:tEndLast));
        end
         
        
        for it = 1:length(trial)
            
            % Number of spikes observed
            y=sum((trial{it}>=tStart).*(trial{it}<tEnd));
            
            % rate from JN
            rateExpected = sum(rateJN(it, tStart:tEnd));
            
            % LL
            LL = LL + log(poisspdf(y, rateExpected));
            
            if ww > 1
                yLast = sum((trial{it}>=tStartLast).*(trial{it}<tEndLast));
                rateExpectedLast = sum(rateJN(it, tStartLast:tEndLast));
                noiseCorr = noiseCorr + (yLast-rateExpectedLast)*(y-rateExpected);
                yVar = yVar + (y-rateExpected)^2;
                nCorr = nCorr + 1;
            end
                
            
            % Repeat for boostrap data
            for ib = 1:nBoot
                yBoot = poissrnd(rateBoot);
                LLBoot(ib) = LLBoot(ib) +  log(poisspdf(yBoot, rateBoot));
                if ww > 1
                    noiseCorrBS(ib) = noiseCorrBS(ib) + (yBoot-rateBoot)*(yBootLast(it,ib)-rateBootLast);
                    yVarBS(ib) = yVarBS(ib) + (yBoot-rateBoot)^2;
                    yBootLast(it, ib) = yBoot;
                end
            end
        end
        nDiscovery = sum(LLBoot > LL);
        Poisson_P(is, ww) = 1.0 - nDiscovery/nBoot;
        if figFlg > 1
            figure(1);
            subplot(2,1,1);
            histogram(LLBoot);
            l = axis;
            hold on;
            plot([LL LL], [l(3) l(4)], 'k--');
            title(sprintf('Spike rate %f', rateBoot*1000.0/ParamModel.Increment));
            pause();
            hold off;
            subplot(2,1,2);
            histogram(noiseCorrBS./yVarBS);
            l = axis;
            hold on;
            plot([noiseCorr/yVar noiseCorr/yVar], [l(3) l(4)], 'k--');
            title('Noise Correlations');
            hold off;
        end
            
    end
    fprintf(1, 'Done with Stim %d/%d - %d/%d Windows Reject Poisson at 5%%\n', is, NbStim, sum(Poisson_P(is,:)<= 0.05), WinNum);
    fprintf(1, '\t Running NoiseCC = %.3f\n', noiseCorr/yVar);
end

noiseCC = noiseCorr/yVar;
noiseCCBS = noiseCorrBS./yVarBS;

% Two tail test
nDiscovery = sum(abs(noiseCC) < noiseCCBS) + sum(-abs(noiseCC) > noiseCCBS);  % number of times BS value is more extreme than actual value
noiseCC_P = nDiscovery/nBoot;

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
       

    plot(reshape(Poisson_P, 1,  NbStim*WinNum), reshape(rateUsed, 1, NbStim*WinNum), '+');
    xlabel('Prob is Poisson');
    ylabel('Rate in 10 ms');
    
    figure();
    histogram(reshape(Poisson_P, 1,  NbStim*WinNum));
    
end
    

end
