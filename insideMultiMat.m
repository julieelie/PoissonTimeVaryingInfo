function [MMD_local] = insideMultiMat(ThreeDDat,MinProbThresh, MinProb,Verbose)

if nargin<2
    MinProbThresh=1;
end
if nargin<3
    % Define a minimum probability under which the proba should be considered 0
    MinProb = 1/(8*365*24*60*60*5); %1/number of 20ms bins in a life of a zebra finch in the lab
end

if nargin<4
    Verbose=0; %Set to 1 to see more outputs
end

if size(ThreeDDat,1)==1
    MMD_local = ThreeDDat{1};
else
    MultiMultDat_old = insideMultiMat(ThreeDDat(1:end-1), MinProbThresh, MinProb);
    
    % Identify rows of ThreeDDat{end} that have all values above MinProb
    % calculus for the other rows is unnecessary
    NbCol = size(ThreeDDat{end},2);
    if MinProbThresh
        Rows2keep = find(sum((ThreeDDat{end}> MinProb), 2)>=NbCol);
        MultiMultDat_new = ThreeDDat{end}(Rows2keep,:);
        NbRow_end = length(Rows2keep);
    else
        MultiMultDat_new = ThreeDDat{end};
        NbRow_end = size(ThreeDDat{end},1);
    end
    NbRow_old = size(MultiMultDat_old,1);
    
    if Verbose
        fprintf('Keeping %d/%d paths from the new matrix', NbRow_end, size(ThreeDDat{end},1));
        [~, memAvail]=system('top -l 1 | head -n 10 | grep PhysMem');
        memAvail = strsplit(memAvail,' ');
        memAvail = str2double(memAvail{6}(1:end-1))/1024;
        fprintf(1, '\n\nAvailable Memory %f\n', memAvail);
        fprintf(1, 'Number of old Paths %d Asking for %d combinations for a total of %d paths\n',  NbRow_old, NbRow_end, NbRow_end * NbRow_old);
        fprintf(1, 'Asking for %d MBytes\n', fix(NbRow_end*NbRow_old*NbCol*8/10^6));
    end
    
    MMD_local = nan(NbRow_end * NbRow_old ,NbCol);

    Dimstep =  NbRow_old;
    for yy=1:NbRow_end
        MMD_local((yy-1)*Dimstep+1 : yy*Dimstep , :) = MultiMultDat_old .* repmat(MultiMultDat_new(yy,:),Dimstep,1);
    end
end

%% Only keep paths (rows) that give mean row>0 (average product of probabilities over stims)
if MinProbThresh
    Paths2keep = find(mean(MMD_local,2) > MinProb);
    MMD_local = MMD_local(Paths2keep,:);
end

end