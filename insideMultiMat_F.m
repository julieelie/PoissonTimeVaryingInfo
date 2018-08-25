function [MatDim,File_new] = insideMultiMat_F(ThreeDDat,Firstwin,StimIndices_AllWin, Stim_local,modnb,InfMat_Storage_location,MinProbThresh, MinProb)
% Work with all the files in a temp folder on the computer machine so that there is no traffic Jam!!!
if nargin <6
    InfMat_Storage_location = '/tmp/LocalTempStorageInfo';
end
if nargin<7
    MinProbThresh=1;
end

if nargin<8
    % Define a minimum probability under which the proba should be considered 0
    MinProb = 1/(8*365*24*60*60*5); %1/number of 20ms bins in a life of a zebra finch in the lab
end

    
%% First step
if size(ThreeDDat,1)==1
    MultiMultDat = ThreeDDat{1};
    File_new.name = sprintf('Temp_CumInf%d_%d_1',modnb, Firstwin);
    File_new.path = InfMat_Storage_location;
    fid=fopen(fullfile(InfMat_Storage_location, File_new.name),'w');
    fwrite(fid, MultiMultDat','double')
    MatDim=size(MultiMultDat');
     % Write down the size of the matrix as the two last digits
    fseek(fid, 0, 'eof');
    fwrite(fid, MatDim,'double');
    fclose(fid);
    
else
    %% We are in the recurrent loop
    % Check if the file we need is already created
    Files_old = dir([InfMat_Storage_location sprintf('/Temp_CumInf%d_%d*',modnb, Firstwin)]);
    Nminus1File=0;
    if ~isempty(Files_old)
        for ff=1:length(Files_old)
            if strcmp(Files_old(ff).name, sprintf('Temp_CumInf%d_%d_%d', modnb, Firstwin, size(ThreeDDat,1)-1))
                Nminus1File=1;
                break
            end
        end
    end
    if Nminus1File
        File_old = Files_old(ff);
        % Open the file from the previous matrix
        fid_old = fopen(fullfile(InfMat_Storage_location,File_old.name));
        fseek(fid_old, -2*8, 'eof');
        MatDim_old = fread(fid_old, 2, 'double');
        % Identify the Stimuli for that window in the old file if it
        % contains more stims than the new
        if MatDim_old(1)>length(Stim_local)
            ID_ind_old = nan(length(Stim_local),1);
            for ss=1:length(Stim_local)
                st = Stim_local(ss);
                ID_ind_old(ss) = find(StimIndices_AllWin{size(ThreeDDat,1)-1}==st);
            end
        else
            ID_ind_old = 1:MatDim_old(1);
        end
    else
    
        [MatDim_old, File_old] = insideMultiMat_F(ThreeDDat(1:end-1),Firstwin,StimIndices_AllWin(1:end-1), Stim_local,modnb,InfMat_Storage_location,MinProbThresh, MinProb);
        % Open the file from the previous matrix
        fid_old = fopen(fullfile(File_old.path,File_old.name));
        ID_ind_old = 1:MatDim_old(1);
    end
    % Open the new file for the new matrix
    File_new.name = sprintf('Temp_CumInf%d_%d_%d', modnb,Firstwin,length(ThreeDDat));
    File_new.path = InfMat_Storage_location;
    % rm the file if it already exists at that point
    system(sprintf('rm %s', fullfile(File_new.path, File_new.name)));
    fid_new = fopen(fullfile(InfMat_Storage_location,File_new.name), 'w');
    
    % Change the orientation of the matrix to have the stimuli as rows and
    % the probability distribution of neural responses as columns
    MultiMultDat_new = ThreeDDat{end}';
    
    % Chunk the old matrix in set of columns of 10^6 or smaller if smaller
    N = floor(MatDim_old(2)/10^6);
    NbCol_chunks = [repmat(10^6, 1,N) MatDim_old(2) - N*10^6];
    NbStim = MatDim_old(1);
    Offsets = [0 cumsum(NbCol_chunks(1:end-1))].*NbStim.*8; % each double number is coded by 8 bytes in the file
            
    
    % Loop through line of new matrix and chuncks of old matrix to form the
    % new matrix
    MatDim1 = length(Stim_local);
    NbCol_new = size(MultiMultDat_new,2);
    MatDim2 = nan(NbCol_new*length(NbCol_chunks),1);
    jj=0;
    for yy=1:NbCol_new
        for cc=1:length(NbCol_chunks)
            % correctly position into the old file and read out chunk of
            % matrix
            fseek(fid_old,Offsets(cc),'bof');
            MultiMultDat_old_local = fread(fid_old,[NbStim NbCol_chunks(cc)], 'double');
            if cc==length(NbCol_chunks)
                % check that we read all the old file
                Current_p = ftell(fid_old);
                fseek(fid_old,-2*8,'eof');%Here minus 2 because we use two doubles to code the size of the matrix
                Endoffile_p = ftell(fid_old); 
                if Current_p ~=Endoffile_p
                    fprintf('Issue here!! in insideMuliMat_F\nthe columns if the matrix %s where not all used for the calculation\n',File_old.name);
                end
            end
            % multiply each column of the chunk of the old matrix with a column
            % from the new matrix
            MMD_local = MultiMultDat_old_local(ID_ind_old,:) .* repmat(MultiMultDat_new(:,yy),1,NbCol_chunks(cc));
            
            % Discard paths (columns) where all stims have a very low proba
            if MinProbThresh
                BadPaths = find(sum((MMD_local < MinProb),1)==NbStim);
                Paths2keep=setdiff(1:size(MMD_local,2), BadPaths);
                MultiMultDat = MMD_local(:,Paths2keep);
            else
                MultiMultDat = MMD_local;
            end
            
            % save the output by appending to the new output file
            fseek(fid_new, 0, 'eof');
            NBE=fwrite(fid_new, MultiMultDat,'double');
            if NBE~=size(MultiMultDat,1)*size(MultiMultDat,2)
                fprintf('Error in write here! in insideMuliMat_F\nNb of elements written=%d when %d should be written %d\n',NBE,size(MultiMultDat,1)*size(MultiMultDat,2));
                fprintf('The last Error message of the file is %s\n', ferror(fid_new));
            end
            % Keep count of the size of the growing matrix in the new
            % output file
            jj=jj+1;
            MatDim2(jj) = size(MultiMultDat,2);
%             if NbStim~=size(MultiMultDat,1)
%                 fprintf('Problem here, the number of rows should always be the same between all matrices!!!\n');
%             end
            %fprintf('In insideMultiMat_F computing line %d/%d of new matrix\nwith chunck %d/%d of old one\n',yy,NbCol_new,cc,length(NbCol_chunks)); 
        end
    end
%     % check that the matrix saved in the new file is of the expected size
%     if NbCol_new*MatDim_old(2)~=MatDim(2)
%         fprintf('Problem of size here, did we forget some combination on the run?!')
%     end
    MatDim=[MatDim1 sum(MatDim2)];
    % Write down the size of the matrix as the two last digits
    fseek(fid_new, 0, 'eof');
    fwrite(fid_new, MatDim,'double');

    % Close old and new files and if wanted, delete old file
    fclose(fid_new);
    fclose(fid_old);
    %system(sprintf('rm %s',fullfile(InfMat_Storage_location,File_old.name)), '-echo')
    
    
end



end