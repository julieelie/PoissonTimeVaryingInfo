function slurm_savio_RunCumInfo_SemanticCells(CIType)
%% Set up the paths, retrieve list of cells already done, list of all cells
addpath('/auto/k1/queued/')
addpath(genpath('/auto/fhome/julie/Code/SingleUnitModels'));
addpath(genpath('/auto/fhome/julie/Code/GeneralCode'));
addpath('/auto/fhome/julie/Code/tlab/src/slurmbot/matlab')
if nargin==0
    CIType = 'CIS'; % CIS=Cumulative information for stimuli, CIC = Cumulative Information for Semantic Categories, CISR= Cumulative information for stimuli with fixed rates
% CICR= Cumulative information for semantic categories with fixed rates,
% CICRand = Cumulative information for random categories % Info =
% Instantaneous Information
end

cd /auto/tdrive/julie/k6/julie/matfile/ModMatInfo/
%DoneFile=dir('InfoPoissonGF_*');

%system('ssh TheunissenLab')
% cd /auto/tdrive/julie/k6/julie/matfile/ModMat
% matlab
% DoneFile=dir('Models_GLMPoisson*');
% save('/auto/tdrive/julie/k6/julie/matfile/DoneFile.mat','DoneFile')
% quit
% logout
% system('scp TheunissenLab:/auto/tdrive/julie/k6/julie/matfile/DoneFile.mat /global/home/users/jelie/MatFiles/DoneFile.mat')
% load('/global/home/users/jelie/MatFiles/DoneFile.mat');

if strcmp(CIType, 'Info') || strcmp(CIType, 'CIS')
    load('/auto/tdrive/julie/NeuralData/SemanticGLMModel/FanoFactor_CoherenceOptPSTHBin_SemCell.mat','List_SemanticCellspath');
    Local_list = List_SemanticCellspath;
elseif strcmp(CIType, 'CIC') || strcmp(CIType, 'CICR') || strcmp(CIType, 'CISR') || strcmp(CIType, 'CICRand') 
    load('/auto/tdrive/julie/k6/julie/matfile/ModMatInfo/600msInfoCumInfoMCJKRRand/ListFilesCumInfoSignif.mat','ValidStimCumInfoFilenames')
    Local_list = ValidStimCumInfoFilenames;
end
% FID = fopen('/auto/tdrive/julie/k6/julie/matfile/ModMatInfo/JobToDoSavio/DoneCells07032017.txt');
% DoneFiles = textscan(FID, '%s');
% DoneFiles = DoneFiles{1};

%% Set up the variables for slurm
JobParams = struct;
JobParams.Name = CIType;
JobParams.Partition = 'savio';
JobParams.Account = 'ac_birdpow';
JobParams.Qos = 'savio_normal';
JobParams.NTasks = 1;
JobParams.CPU = 20;
JobParams.Type = CIType;
if strcmp(CIType, 'Info')
    jobParams.TimeLimit = '00:05:00';
else
    jobParams.TimeLimit = '72:00:00';
end
SlurmParams.cmd = 'Semantic_NeuroInfo_Poisson_savio(''%s'', ''%s'');';
SlurmParams.resultsDirectory='/global/scratch/jelie/MatFiles/ModMatInfo';

%% Set up variables to identify cells to run and order
MatfileToDo = cell(length(Local_list),1);
MatNameToDo = cell(length(Local_list),1);

%% Create jobs' files
cd /auto/tdrive/julie/k6/julie/matfile/ModMatInfo/JobToDoSavio
if strcmp(CIType, 'CIS')
    for ff=1:length(Local_list)
        fprintf(1,'checking file %d/%d\n',ff,length(Local_list));
        [P,TheFile,ext]=fileparts(Local_list{ff});
        if ~isempty(strfind(TheFile, 'InfoPoissonKDEF'))
            MatfileToDo{ff}= fullfile('/auto/tdrive/julie/k6/julie/matfile/FirstVoc1sMat',['FirstVoc1s' TheFile(16:end) ext]);
            MatNameToDo{ff}=['FirstVoc1s' TheFile(16:end) ext];
        else
            BegPath=strfind(P,'k6');
            MatfileToDo{ff}= fullfile('/auto/tdrive/julie/k6/julie/matfile/FirstVoc1sMat',['FirstVoc1s' TheFile(8:end) ext]);
            MatNameToDo{ff}=['FirstVoc1s' TheFile(8:end) ext];
        end
    %     DC = 0;
    %     for dd=1:length(DoneFiles)
    %         Var_local = strfind(DoneFiles(dd), MatNameToDo{ff});
    %         if ~isempty(Var_local{1})
    %             DC=1;
    %         end
    %     end
    %     if ~DC
            JobParams.Name = MatNameToDo{ff};
            JobParams.out = fullfile(SlurmParams.resultsDirectory,sprintf('slurm_out_%s_%s_%%j.txt', JobParams.Name, JobParams.Type));
            JobParams.err = JobParams.out;
            icmd = sprintf(SlurmParams.cmd, MatfileToDo{ff},  JobParams.Type);
            fprintf(1,'creating file slurm_sbatch with command %s\n',icmd);
            slurm_sbatch_savio(icmd,JobParams);
    %     else
    %         fprintf(1, '!!!! File already done!!!!');
    %     end
    end
elseif strcmp(CIType, 'CIC') || strcmp(CIType, 'CISR') || strcmp(CIType, 'CICR') || strcmp(CIType, 'CICRand')
     for ff=1:length(Local_list)
        fprintf(1,'checking file %d/%d\n',ff,length(Local_list));
        TheFile=Local_list(ff).name;
        Ind = strfind(TheFile, '_CIS');
        MatfileToDo{ff}= fullfile('/auto/tdrive/julie/k6/julie/matfile/FirstVoc1sMat',['FirstVoc1s' TheFile(16:(Ind-1)) '.mat']);
        MatNameToDo{ff}=['FirstVoc1s' TheFile(16:(Ind-1)) '.mat'];
        
        JobParams.Name = MatNameToDo{ff};
        JobParams.out = fullfile(SlurmParams.resultsDirectory,sprintf('slurm_out_%s_%s_%%j.txt', JobParams.Name, JobParams.Type));
        JobParams.err = JobParams.out;
        icmd = sprintf(SlurmParams.cmd, MatfileToDo{ff},  JobParams.Type);
        fprintf(1,'creating file slurm_sbatch with command %s\n',icmd);
        slurm_sbatch_savio(icmd,JobParams);
     end
end
    
fprintf(1,'DONE Creating all Jobs'' files!!!\n');

% %% Find the jobs of the interesting cells
% IndicesHipCells = nan(4,1);
% IndicesHipCells(1) = IndOrd(1); %Highest spike rate cell Site2_L1100R1450_e14_s0_ss1 %Ac Cell
% IndicesHipCells(2) = IndOrd(end); % lowest spike rate cell
% for ff=1:length(MatNameToDo)
% %     if strcmp(MatNameToDo{ff}(10:end-4), 'Site2_L1000R900_e13_s0_ss1') %SYN Cell cannot be found in the dataset of single semantic units...
% %         IndicesHipCells(5)=ff;
%     if strcmp(MatNameToDo{ff}(10:end-4), 'Site3_L2500R2300_e22_s1_ss1')%Ag127
%         IndicesHipCells(3)=ff;
%     elseif strcmp(MatNameToDo{ff}(10:end-4), 'Site4_L1500R1900_e23_s0_ss2') %DC128
%         IndicesHipCells(4)=ff;
%     end
% end
% 
% 
% 
% for jj=1:length(IndicesHipCells)
%     IndexFile=IndicesHipCells(jj)
%     JobParams.out = fullfile(SlurmParams.resultsDirectory,sprintf('slurm_out_%s_%%j.txt', MatNameToDo{IndexFile}));
%     JobParams.err = JobParams.out;
%     icmd = sprintf(SlurmParams.cmd, MatfileToDo{IndexFile});
%     fprintf(1,'creating file slurm_sbatch with command %s\n',icmd);
%     slurm_sbatch_savio(icmd,JobParams);
% end
% fprintf(1,'DONE Creating all Jobs'' files!!!\n');

   