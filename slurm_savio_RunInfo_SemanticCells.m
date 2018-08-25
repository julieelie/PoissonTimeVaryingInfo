%% Set up the paths, retrieve list of cells already done, list of all cells
addpath('/auto/k1/queued/')
addpath(genpath('/auto/fhome/julie/Code/SingleUnitModels'));
addpath(genpath('/auto/fhome/julie/Code/GeneralCode'));
addpath('/auto/fhome/julie/Code/tlab/src/slurmbot/matlab')

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

load('/auto/tdrive/julie/NeuralData/SemanticGLMModel/FanoFactor_CoherenceOptPSTHBin_SemCell.mat','List_SemanticCellspath');

%% Set up the variables for slurm
JobParams = struct;
JobParams.Name = 'Info';
JobParams.Partition = 'savio2';
JobParams.Account = 'fc_birdpow';
JobParams.Qos = 'savio_normal';
JobParams.NTasks = 1;
JobParams.CPU = 20;
SlurmParams.cmd = 'Semantic_NeuroInfo_Poisson_savio(''%s'');';
SlurmParams.resultsDirectory='/global/scratch/jelie/MatFiles/ModMatInfo';

%% Set up variables to identify cells to run and order
MatfileToDo = cell(length(List_SemanticCellspath),1);
MatNameToDo = cell(length(List_SemanticCellspath),1);

%% Create jobs' files
cd /auto/tdrive/julie/k6/julie/matfile/ModMatInfo/JobToDoSavio
for ff=1:length(List_SemanticCellspath)
    fprintf(1,'checking file %d/%d\n',ff,length(List_SemanticCellspath));
    [P,TheFile,ext]=fileparts(List_SemanticCellspath{ff});
    BegPath=strfind(P,'k6');
    MatfileToDo{ff}= fullfile('/auto/tdrive/julie/k6/julie/matfile/FirstVoc1sMat',['FirstVoc1s' TheFile(8:end) ext]);
    MatNameToDo{ff}=['FirstVoc1s' TheFile(8:end) ext];
    JobParams.out = fullfile(SlurmParams.resultsDirectory,sprintf('slurm_out_%s_%%j.txt', MatNameToDo{ff}));
    JobParams.err = JobParams.out;
    icmd = sprintf(SlurmParams.cmd, MatfileToDo{ff});
    fprintf(1,'creating file slurm_sbatch with command %s\n',icmd);
    slurm_sbatch_savio(icmd,JobParams);
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

   