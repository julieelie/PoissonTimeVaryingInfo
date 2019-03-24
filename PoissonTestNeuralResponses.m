function [Poisson_P, noiseCC, noiseCC_P, noiseCCA, noiseCCA_P, noiseCCC, noiseCCC_P, noiseCCCA, noiseCCCA_P] = PoissonTestNeuralResponses(Res,ParamModel)

% Global argurments
figFlg=1;     % Flag for plotting - 1 summary plots, 2 debugging plots + summary plots
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

% Cumulative variables
noiseCorr = 0.0;
yVar = 0.0;
noiseCorrA = 0.0;
yVarA = 0.0;

% Same for corrected for trial variations
noiseCorrC = 0.0;
yVarC = 0.0;
noiseCorrAC = 0.0;
yVarAC = 0.0;

noiseCorrBS = zeros(1, nBoot);
yVarBS = zeros(1, nBoot);
nCorr = 0;


%% Loop through stims and windows

for is = 1:NbStim
    trial = Res.Trials{DataSel(is)};                  % Raw trial data
    rateJN = Res.JackKnife_KDE_Filtered{DataSel(is)};   % JN rate data
    yBootLast = zeros(length(trial), nBoot);            % to store the values of the last random data
    
    indLast = min([length(Res.PSTH{DataSel(is)}), length(Res.PSTH_KDE_Filtered{DataSel(is)}), length(rateJN)]);
    
    % Rate for entire period
    tBeg = Wins(1);
    tLast = Wins(WinNum) + ParamModel.Increment;
    if (tLast > indLast) 
        tLast = indLast;
    end
    
    rateActualTrial = sum(Res.PSTH{DataSel(is)}(tBeg:tLast));
    yTrial = zeros(1, length(trial));
    rateJNTrial = zeros(1, length(trial));
    
    for it = 1:length(trial)
        yTrial(it) = sum((trial{it}>=tBeg).*(trial{it}<tLast));
        rateJNTrial(it) = sum(rateJN(it, tBeg:tLast));
    end
    yCorrTrial = (yTrial-rateJNTrial)./rateJNTrial;            % Correction factor for trial effects in % using the JN rate
    yCorrATrial = (yTrial-rateActualTrial)/rateActualTrial;    % using actual rate
    
    
    for ww = 1:WinNum          
        tStart = Wins(ww);
        tEnd = tStart + ParamModel.Increment;
        
        if (tEnd > indLast) 
            break;
        end
        
        LL = 0;   % Log likelihood of data
        LLBoot = zeros(1,nBoot); % Log Likelihood bootstrapped
        rateBoot = sum(Res.PSTH_KDE_Filtered{DataSel(is)}(tStart:tEnd));
        rateActual = sum(Res.PSTH{DataSel(is)}(tStart:tEnd));
        
        if ww > 1
            tStartLast = Wins(ww-1);
            tEndLast = tStartLast + ParamModel.Increment;
            rateBootLast = sum(Res.PSTH_KDE_Filtered{DataSel(is)}(tStartLast:tEndLast));
            rateActualLast = sum(Res.PSTH{DataSel(is)}(tStartLast:tEndLast));
        end
         
        
        for it = 1:length(trial)
            
            % Number of spikes observed
            y=sum((trial{it}>=tStart).*(trial{it}<tEnd));
                        
            % rate from JN
            rateExpected = sum(rateJN(it, tStart:tEnd));
            
            % Corrected values of spike count
            yC = y - yCorrTrial(it)*rateExpected;
            yCA = y - yCorrATrial(it)*rateActual;
            
            % LL
            LL = LL + log(poisspdf(y, rateExpected));
            
            if ww > 1                
                yLast = sum((trial{it}>=tStartLast).*(trial{it}<tEndLast));                                              
                
                % Noise correlation using JN KDE               
                rateExpectedLast = sum(rateJN(it, tStartLast:tEndLast)); 
                yCLast = yLast - yCorrTrial(it)*rateExpectedLast;                
                noiseCorr = noiseCorr + (yLast-rateExpectedLast)*(y-rateExpected);
                yVar = yVar + (y-rateExpected)^2;
                noiseCorrC = noiseCorrC + (yCLast-rateExpectedLast)*(yC-rateExpected);
                yVarC = yVarC + (yC-rateExpected)^2;
                
                % Noise correlation using actual rate
                yCALast = yLast - yCorrATrial(it)*rateActualLast;
                noiseCorrA = noiseCorrA + (yLast-rateActualLast)*(y-rateActual);
                yVarA = yVarA + (y-rateActual)^2;
                noiseCorrAC = noiseCorrAC + (yCALast-rateActualLast)*(yCA-rateActual);
                yVarAC = yVarAC + (yC-rateActual)^2;
                
                nCorr = nCorr + 1;
            end
                
            
            % Repeat for boostrap data
            yBoot = poissrnd(rateBoot, 1, nBoot);
            LLBoot = LLBoot + log(poisspdf(yBoot, rateBoot));
            if ww > 1
                noiseCorrBS = noiseCorrBS + (yBoot-rateBoot).*(yBootLast(it,:)-rateBootLast);
                yVarBS = yVarBS + (yBoot-rateBoot).^2;
            end
            yBootLast(it, :) = yBoot;
            
%             for ib = 1:nBoot               
%                 LLBoot(ib) = LLBoot(ib) +  log(poisspdf(yBoot(ib), rateBoot));
%                 if ww > 1
%                     noiseCorrBS(ib) = noiseCorrBS(ib) + (yBoot(ib)-rateBoot)*(yBootLast(it,ib)-rateBootLast);
%                     yVarBS(ib) = yVarBS(ib) + (yBoot(ib)-rateBoot)^2;                   
%                 end
%                 yBootLast(it, ib) = yBoot(ib);
%             end
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
            hold off;
            
            subplot(2,1,2);
            histogram(noiseCorrBS./yVarBS);
            l = axis;
            hold on;
            plot([noiseCorr/yVar noiseCorr/yVar], [l(3) l(4)], 'k--');
            plot([noiseCorrA/yVarA noiseCorrA/yVarA], [l(3) l(4)], 'r--');
            plot([noiseCorrC/yVarC noiseCorrC/yVarC], [l(3) l(4)], 'g--');
            plot([noiseCorrAC/yVarAC noiseCorrAC/yVarAC], [l(3) l(4)], 'b--');
            title('Noise Correlations');
            hold off;
            
            pause();
        end
            
    end
    fprintf(1, 'Stim %d/%d - %d/%d Reject Poisson\n', is, NbStim, sum(Poisson_P(is,:)<= 0.05), WinNum);
    fprintf(1, '\t NoiseCC = %.3f NoiseCCC = %.3f\n', noiseCorr/yVar, noiseCorrC/yVarC);
end

noiseCC = noiseCorr/yVar;                % Noise Correlation with JN mean rate (excluding the trial)
noiseCCA = noiseCorrA/yVarA;             % Noise Correlation with Actual rate (using all trial)
noiseCCC = noiseCorrC/yVarC;             % Noise Correlation with JN rate + correction for trial variability
noiseCCCA = noiseCorrAC/yVarAC;          % Noise Correlation with Actual rate  + correction for trial variability

noiseCCBS = noiseCorrBS./yVarBS;

% Two tail tests
nDiscovery = sum(abs(noiseCC) < noiseCCBS) + sum(-abs(noiseCC) > noiseCCBS);  % number of times BS value is more extreme than actual value
noiseCC_P = nDiscovery/nBoot;

nDiscovery = sum(abs(noiseCCA) < noiseCCBS) + sum(-abs(noiseCCA) > noiseCCBS); 
noiseCCA_P = nDiscovery/nBoot;

nDiscovery = sum(abs(noiseCCC) < noiseCCBS) + sum(-abs(noiseCCC) > noiseCCBS);  % number of times BS value is more extreme than actual value
noiseCCC_P = nDiscovery/nBoot;

nDiscovery = sum(abs(noiseCCCA) < noiseCCBS) + sum(-abs(noiseCCCA) > noiseCCBS); 
noiseCCCA_P = nDiscovery/nBoot;
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
    
    figure();
    histogram(noiseCCBS);
    l = axis;
    hold on;
    plot([noiseCC noiseCC], [l(3) l(4)], 'k--');
    plot([noiseCCA noiseCCA], [l(3) l(4)], 'r--');
    plot([noiseCCC noiseCCC], [l(3) l(4)], 'g--');
    plot([noiseCCCA noiseCCCA], [l(3) l(4)], 'b--');
    title(sprintf('Noise Correlations P=%.4f P=%.4f', noiseCC_P, noiseCCA_P));
    hold off;
    
end
    

end

