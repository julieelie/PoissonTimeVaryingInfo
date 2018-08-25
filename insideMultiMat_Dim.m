function [MultiMultDat,MultiMultDatprev] = insideMultiMat_Dim(ThreeDDat,MultiMultDatprev,MinProbThresh, MinProb)
% The output is organize such that first dimension corresponds to stimuli,
% second dimension to p(y(t-size(ThreeDDat,1))), third dimension to
% p(y(t-size(ThreeDDat,1)+1))...


if nargin<2
    MultiMultDatprev={};
end
if nargin<3
    % Define a minimum probability under which the proba should be considered 0
    MinProb = 1/(8*365*24*60*60*5); %1/number of 20ms bins in a life of a zebra finch in the lab
    if iscell(MultiMultDatprev)
        MinProbThresh=1;
    else
        MinProbThresh=MultiMultDatprev;
        MultiMultDatprev={};
    end
end
    
if nargin<4
    if iscell(MultiMultDatprev)
        MinProb = MinProbThresh;
        MinProbThresh=1;
    else
        MinProb = MinProbThresh;
        MinProbThresh=MultiMultDatprev;
        MultiMultDatprev={};
    end
end

if size(ThreeDDat,1)==1
    MultiMultDat{1} = ThreeDDat{1}';%here the column and rows of the first matrix are swap because rows should be the stimuli
    MultiMultDatprev = {};
else
    if isempty(MultiMultDatprev)
        [MultiMultDatprev,~] = insideMultiMat_Dim(ThreeDDat(1:end-1),MinProbThresh, MinProb);
    end
    
    NbStim = size(ThreeDDat{end},2);
    if MinProbThresh
        % Identify rows of ThreeDDat{end} that have all values above MinProb
        % calculus for the other rows is unnecessary
        Rows2getrid = find(sum((ThreeDDat{end}< MinProb), 2)>=NbStim);
        Rows2keep = setdiff(1:size(ThreeDDat{end},1),Rows2getrid);
        MultiMultDat_new = ThreeDDat{end}(Rows2keep,:)';% swaping rows and column to get stims as rows
        SizeDim_end = length(Rows2keep);%Now this is the number of column of MultiMultDat_new
    else
        MultiMultDat_new = ThreeDDat{end}';% swaping rows and column to get stims as rows
        SizeDim_end = size(ThreeDDat{end},1);%Now this is the number of column of MultiMultDat_new
    end
    SizeDim_old = size(MultiMultDatprev{1});
    if SizeDim_old(1) ~= NbStim
        fprintf('problem here the nb of stim (rows) in MultiMultDat_old is %d and %d in MultiMultDat_new', SizeDim_old(1), NbStim);
    end
    
    %     fprintf('Keeping %d/%d paths from the new matrix', SizeDim_end, size(ThreeDDat{end},1));
    %     [~, memAvail]=system('top -l 1 | head -n 10 | grep PhysMem');
    %     memAvail = strsplit(memAvail,' ');
    %     memAvail = str2double(memAvail{6}(1:end-1))/1024;
    %     fprintf(1, '\n\nAvailable Memory %f\n', memAvail);
    %     fprintf(1, 'Number of old Paths %d Asking for %d combinations for a total of %d paths\n',  SizeDim_old, SizeDim_end, SizeDim_end * SizeDim_old);
    %     fprintf(1, 'Asking for %d MBytes\n', fix(SizeDim_end*SizeDim_old*NbStim*8/10^6));
    %
    % Now calculate the new matrix
    MultiMultDat = {nan([SizeDim_old SizeDim_end])};
    StepSize = prod(SizeDim_old);
    Indices_yy = flip(1:SizeDim_end);
    for yy=1:SizeDim_end
        yy_local = Indices_yy(yy);
        MultiMultDat{1}((end-yy_local*StepSize+1) :(end - (yy_local-1)*StepSize)) = MultiMultDatprev{1} .* repmat(MultiMultDat_new(:,yy),[1,SizeDim_old(2:end)]);
    end
    
end

% %% Only keep paths (columns) that give mean row>0 (average product of probabilities over stims)
% Paths2keep=cell(ndims(MMD_local)-1,1);
% [Paths2keep{:}] = ind2sub([SizeDim_old(2:end) SizeDim_end],find(mean(MMD_local,1) > MinProb));
% MultiMultDat = MMD_local(:,Paths2keep{:});
% clear MMD_local

end