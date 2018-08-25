function [PG_Index, FanoFactor_Index, Wins, Ymean, Yvar] = PoissonGaussianNeuralResponses(Trials,ParamModel,SWITCH, Cellname, Spectro)
FIG=0;
if nargin<2
    ParamModel = struct();
end
if  ~isfield(ParamModel,'MinWin') || isempty(ParamModel.MinWin)
    ParamModel.MinWin = 10; % end point of the first analysis window (spectrogram and neural response)
end
if ~isfield(ParamModel,'MaxWin') || isempty(ParamModel.MaxWin)
    ParamModel.MaxWin = 600; %end point of the last anaysis window for...
    ... neural response and end point of the largest analysis window for...
        ... spectrogram
end
if ~isfield(ParamModel,'Increment') || isempty(ParamModel.Increment)
    ParamModel.Increment = 10; %increase the size of the spectro window with a Xms pace
end
if ~isfield(ParamModel,'NeuroBin') || isempty(ParamModel.NeuroBin)
    ParamModel.NeuroBin = 10; % size of the window (ms) within which the neural response is analyzed
                               % The end of the window of analysis is
                               % determined by the Increment and ResDelay (see below).
end
if ~isfield(ParamModel,'ResDelay') || isempty(ParamModel.ResDelay)
    ParamModel.ResDelay = 0; % Delay in ms between the end of the...
    ... spectrogram window and the end of the neural response window
end

if nargin<3
    SWITCH = struct();
end
if ~isfield(SWITCH,'Models') || isempty(SWITCH.Models)
    SWITCH.Models=0;
end

if nargin<4
    Cellname = 'Your Cell';
end

if nargin<5
    if SWITCH.Models
        fprintf(1,'Please specify the spectrograms of the sound stimuli\n');
        return
    else
        Spectro = [];
    end
end

% define the list of end points of spectrogram and neural responses windows
Wins = ParamModel.MinWin:ParamModel.Increment:ParamModel.MaxWin;

% # of models to run on the data
WinNum = length(Wins);

% Number of stims in the data set
NbStim = length(Trials);

%% Initialize output variables
PG_Index = nan(WinNum,1);
FanoFactor_Index = nan(WinNum,1);
Ymean = nan(WinNum,NbStim);
Yvar = nan(WinNum,NbStim);

%% Now loop through window sizes and look at spike rate distributions
for ww = 1:WinNum
    %fprintf(1,'%d/%d models\n', mm, modNum);
    Win = Wins(ww);
    
    %% define new dataset depending on the size of the window of the model
    % loop through the stims and only keep the Win first ms of them when
    % they are longer than Win ms or disgard
    if SWITCH.Models
        duration = nan(NbStim,1);
        for ss = 1:NbStim
            duration(ss)=Spectro.to{ss}(end)*1000; %converting s in ms here
        end
        Stim_local = find(duration >= (Win+ParamModel.ResDelay));% here we add ResDelay because we need to get sounds with corresponding psth that go ResDelay beyond the spectrogram of size Win
        NbStim_local = length(Stim_local);
        if NbStim_local<20
            sprintf('Only %d stims long enough to run the model: no model is run with window size %dms\n', NbStim_local, Win);
            break
        end
    else
        Stim_local = 1:NbStim;
        NbStim_local = NbStim;
    end
    
    % Initializing outputs for spike rate/count
    y=cell(NbStim_local,1);
    ymean=nan(NbStim_local,1);
    yvar=ymean;
        
    
    for ss = 1:NbStim_local
        dd=Stim_local(ss);
        
        % Values of max spike rate(y), mean spike rate (y) and exact number of spike per trial (y) within the window
        FirstTimePoint = Win - ParamModel.NeuroBin+ ParamModel.ResDelay;
        LastTimePoint = Win + ParamModel.ResDelay;
        y{ss}=nan(length(Trials{dd}),1);
        for tt=1:length(Trials{dd})
            y{ss}(tt)=sum((Trials{dd}{tt}>=FirstTimePoint).*(Trials{dd}{tt}<LastTimePoint));
        end
        ymean(ss)=mean(y{ss});
        yvar(ss)=var(y{ss});
    end
    Ymean(ww,:) = ymean;
    Yvar(ww,:) = yvar;
    
    % Investigate how poisson or gaussian neural responses are for this
    % neuron
    MAX=max(max(yvar),max(ymean));
    PG_Index(ww) = sum(power(yvar-repmat(mean(yvar),length(yvar),1),2))/sum(power(yvar-ymean,2));
    yvar_local = yvar;
    ymean_local = ymean;
    yvar_local(yvar==ymean)=1;
    ymean_local(yvar==ymean)=1;
    FanoFactor_Index(ww) = mean(yvar_local./ymean_local);
    if FIG>0
        figure(1)
        plot(ymean,yvar,'r.', 'MarkerSize',40)
        ylabel('Variance spike counts per stim')
        xlabel('mean spike count per stim')
        title(sprintf('Win=%d PoissonGaussian Index=%f\n FanoFactor=%f\n',Win,PG_Index(ww),FanoFactor_Index(ww)))
        hold on
        line([0 MAX], [0 MAX]);
        hold off
        pause()
    end    
         
end
if FIG>0
    figure()
    subplot(1,2,1)
    plot(1:WinNum,log2(PG_Index))
    set(gca,'XTick',1:WinNum)
    set(gca,'XTickLabel',Wins)
    ylabel('log2(PG_Index) >0 Poisson <0 Non-Poisson')
    xlabel('windows in ms')
    title(sprintf('SS Errors ratio Non-Poisson/Poisson\n%s',Cellname));
    line([0 WinNum], [0 0])
    subplot(1,2,2)
    plot(1:WinNum,FanoFactor_Index)
    set(gca,'XTick',1:WinNum)
    set(gca,'XTickLabel',Wins)
    ylabel('FanoFactor')
    xlabel('windows in ms')
    title(sprintf('Fano Factor var/mean\n%s',Cellname));
    pause(1)
end
end

