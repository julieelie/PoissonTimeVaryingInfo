%load('/Users/elie/Documents/CODE/data/matfile/ModMatInfo/Models_InfoPoisson_Site4_L1500R1900_e23_s0_ss2_win16_exactCumInfo.mat')
DataMkv1=load('/Users/elie/Documents/CODE/data/matfile/ModMatInfo/Models_InfoPoisson_Site4_L1500R1900_e23_s0_ss2_MarkovHistory1_1.mat','Data');
%R=load('Models_GLMPoisson_Site4_L1500R1900_e23_s0_ss2.mat')
load('/Users/elie/Documents/CODE/data/matfile/ModMatInfo/Models_GLMPoisson_Site4_L1500R1900_e23_s0_ss2_InfoOnly_ExactCumInfo.mat')
% Find the number of windows that run
WinMax=sum(~isnan(Model.Ceiling.info));
%WinFirsts=[1 6 11];
% Calculate the mean spike rate
Mean_SR = nan(19,1);
for ww=1:19
    Mean_SR(ww)=mean(DataMkv1.Data.y_wholeset{ww}*50);%data are the number of spike in 20ms windows so *50 to obtain nb spikes per s
end

figure(1)
% Information
subplot(1,3,1)
Info_MatrixPlot_y_left = [Model.Ceiling.info(1:WinMax), Model.Semantic.info(1:WinMax), Model.Floor.info(1:WinMax)];
Info_MatrixPlot_x_left = repmat((1:WinMax)',1,3);
SPA_Plot_y_right = Mean_SR(1:WinMax);
SPA_Plot_x_right = 1:WinMax;
LegendInfo.CellType = 'One Cell';
LegendInfo.YleftAxis = 'Information (bits)';
myplotyyInfoSPAsemOnly(Info_MatrixPlot_x_left,Info_MatrixPlot_y_left,SPA_Plot_x_right,SPA_Plot_y_right,LegendInfo,Wins,[-0.1 0.7])

% Cumulative Information
subplot(1,3,2)
Info_MatrixPlot_y_left = [Model.Ceiling.cum_info(1:WinMax), Model.Semantic.cum_info(1:WinMax), Model.Floor.cum_info(1:WinMax)];
Info_MatrixPlot_x_left = repmat((1:WinMax)',1,3);
SPA_Plot_y_right = Mean_SR(1:WinMax);
SPA_Plot_x_right = 1:WinMax;
LegendInfo.CellType = 'One Cell';
LegendInfo.YleftAxis = 'Cumulative Information (bits)';
myplotyyInfoSPAsemOnly(Info_MatrixPlot_x_left,Info_MatrixPlot_y_left,SPA_Plot_x_right,SPA_Plot_y_right,LegendInfo,Wins,[-0.1 2])

% Cumulative sum of information
subplot(1,3,3)
Info_MatrixPlot_y_left = [cumsum(Model.Ceiling.info(1:WinMax)), cumsum(Model.Semantic.info(1:WinMax)), cumsum(Model.Floor.info(1:WinMax))];
Info_MatrixPlot_x_left = repmat((1:WinMax)',1,3);
SPA_Plot_y_right = Mean_SR(1:WinMax);
SPA_Plot_x_right = 1:WinMax;
LegendInfo.CellType = 'One Cell';
LegendInfo.YleftAxis = 'Cumulative Sum of Information (bits)';
myplotyyInfoSPAsemOnly(Info_MatrixPlot_x_left,Info_MatrixPlot_y_left,SPA_Plot_x_right,SPA_Plot_y_right,LegendInfo,Wins,[-0.1 2])


figure(3)
plot(1:WinMax, Model.Ceiling.info(1:WinMax), 'g', 1:WinMax, Model.Floor.info(1:WinMax), 'k', 1:WinMax, Model.Semantic.info(1:WinMax), 'r')
legend('Ceiling', 'Floor', 'Semantic', 'Location','NorthWest')
xlabel('Time (ms)')
ylabel('Information (bits)')
IndLab=get(gca,'XTickLabel');
XTickLab = nan(length(IndLab),1);
for ii=1:length(XTickLab)
    XTickLab(ii) = (str2num(IndLab{ii})+1)*10;
end
set(gca,'XTickLabel', XTickLab)

figure(4)
plot(1:WinMax, Model.Ceiling.cum_info(1:WinMax), 'g', 1:WinMax, Model.Floor.cum_info(1:WinMax), 'k', 1:WinMax, Model.Semantic.cum_info(1:WinMax), 'r')
hold on
plot(1:WinMax, cumsum(Model.Ceiling.info(1:WinMax)), 'g--', 1:WinMax, cumsum(Model.Floor.info(1:WinMax)), 'k--', 1:WinMax, cumsum(Model.Semantic.info(1:WinMax)), 'r--')
legend('CI Ceiling', 'CI Floor', 'CI Semantic', 'CSI Ceiling', 'CSI Floor', 'CSI Semantic', 'Location','NorthWest')
xlabel('Time (ms)')
ylabel('Cumulative Information (CI, bits) and Cumulative Sum of Information (CSI, bits)')
IndLab=get(gca,'XTickLabel');
XTickLab = nan(length(IndLab),1);
for ii=1:length(XTickLab)
    XTickLab(ii) = (str2num(IndLab{ii})+1)*10;
end
set(gca,'XTickLabel', XTickLab)

Mkv1=load('/Users/elie/Documents/CODE/data/matfile/ModMatInfo/Models_InfoPoisson_Site4_L1500R1900_e23_s0_ss2_MarkovHistory1_1.mat','Model');
Mkv2=load('/Users/elie/Documents/CODE/data/matfile/ModMatInfo/Models_InfoPoisson_Site4_L1500R1900_e23_s0_ss2_MarkovHistory2_1.mat','Model');
Mkv3=load('/Users/elie/Documents/CODE/data/matfile/ModMatInfo/Models_InfoPoisson_Site4_L1500R1900_e23_s0_ss2_MarkovHistory3_1.mat','Model');
Mkv4=load('/Users/elie/Documents/CODE/data/matfile/ModMatInfo/Models_InfoPoisson_Site4_L1500R1900_e23_s0_ss2_MarkovHistory4_1.mat','Model');
Mkv4_2=load('/Users/elie/Documents/CODE/data/matfile/ModMatInfo/Models_InfoPoisson_Site4_L1500R1900_e23_s0_ss2_MarkovHistory4_2.mat','Model');
Mkv5_2=load('/Users/elie/Documents/CODE/data/matfile/ModMatInfo/Models_InfoPoisson_Site4_L1500R1900_e23_s0_ss2_MarkovHistory5_2.mat','Model');
Mkv6_2=load('/Users/elie/Documents/CODE/data/matfile/ModMatInfo/Models_InfoPoisson_Site4_L1500R1900_e23_s0_ss2_MarkovHistory6_2.mat','Model');
Mkv7_2=load('/Users/elie/Documents/CODE/data/matfile/ModMatInfo/Models_InfoPoisson_Site4_L1500R1900_e23_s0_ss2_MarkovHistory7_2.mat','Model');
Mkv8_2=load('/Users/elie/Documents/CODE/data/matfile/ModMatInfo/Models_InfoPoisson_Site4_L1500R1900_e23_s0_ss2_MarkovHistory8_2.mat','Model');
Mkv9_2=load('/Users/elie/Documents/CODE/data/matfile/ModMatInfo/Models_InfoPoisson_Site4_L1500R1900_e23_s0_ss2_MarkovHistory9_2.mat','Model');

figure(5)
plot(1:WinMax, Model.Ceiling.cum_info(1:WinMax), 'g',1:WinMax, cumsum(Model.Ceiling.info(1:WinMax)), 'g--')
hold on
plot(1:length(Mkv1.Model.Ceiling.cum_info1), [Mkv1.Model.Ceiling.cum_info1 Mkv2.Model.Ceiling.cum_info1  Mkv3.Model.Ceiling.cum_info1  Mkv4.Model.Ceiling.cum_info1])
legend('CI Ceiling', 'CSI Ceiling', 'CI Ceiling M1', 'CI Ceiling M2', 'CI Ceiling M3', 'CI Ceiling M4', 'Location','NorthWest')
xlabel('Time (ms)')
ylabel('Cumulative Information (CI, bits) and Cumulative Sum of Information (CSI, bits)')
IndLab=get(gca,'XTickLabel');
XTickLab = nan(length(IndLab),1);
for ii=1:length(XTickLab)
    XTickLab(ii) = (str2num(IndLab{ii})+1)*10;
end
set(gca,'XTickLabel', XTickLab)

figure(6)
plot(1:WinMax, Model.Ceiling.cum_info(1:WinMax), 'g',1:WinMax, cumsum(Model.Ceiling.info(1:WinMax)), 'g--')
hold on
plot(1:length(Mkv4_2.Model.Ceiling.cum_info1), [Mkv4_2.Model.Ceiling.cum_info1 Mkv5_2.Model.Ceiling.cum_info1  Mkv6_2.Model.Ceiling.cum_info1  Mkv7_2.Model.Ceiling.cum_info1 Mkv8_2.Model.Ceiling.cum_info1])
legend('CI Ceiling', 'CSI Ceiling', 'CI Ceiling M4_2', 'CI Ceiling M5_2', 'CI Ceiling M6_2', 'CI Ceiling M7_2', 'CI Ceiling M8_2', 'Location','NorthWest')
xlabel('Time (ms)')
ylabel('Cumulative Information (CI, bits) and Cumulative Sum of Information (CSI, bits)')
IndLab=get(gca,'XTickLabel');
XTickLab = nan(length(IndLab),1);
for ii=1:length(XTickLab)
    XTickLab(ii) = (str2num(IndLab{ii})+1)*10;
end
set(gca,'XTickLabel', XTickLab)

% figure(4)
% plot(1:WinMax, Model.Ceiling.cum_info1(1:WinMax), 'g', 1:WinMax, Model.Floor.cum_info1(1:WinMax), 'k', 1:WinMax, Model.Semantic.cum_info1(1:WinMax), 'r')
% hold on
% plot(1:WinMax, cumsum(Model.Ceiling.info(1:WinMax)), 'g--', 1:WinMax, cumsum(Model.Floor.info(1:WinMax)), 'k--', 1:WinMax, cumsum(Model.Semantic.info(1:WinMax)), 'r--')
% legend('CI Ceiling', 'CI Floor', 'CI Semantic', 'CSI Ceiling', 'CSI Floor', 'CSI Semantic', 'Location','NorthWest')
% hold on
% for ww=2:length(WinFirsts)
%     WinFirst=(WinFirsts(ww));
%     plot(WinFirst:WinMax, Model.Ceiling.(sprintf('cum_info%d',WinFirst))(1:(WinMax-WinFirst+1)), 'g', WinFirst:WinMax, Model.Floor.(sprintf('cum_info%d',WinFirst))(1:(WinMax-WinFirst+1)), 'k', WinFirst:WinMax, Model.Semantic.(sprintf('cum_info%d',WinFirst))(1:(WinMax-WinFirst+1)), 'r')
%     hold on
% end
% xlabel('Time (ms)')
% ylabel('Cumulative Information (CI, bits) and Cumulative Sum of Information (CSI, bits)')
% IndLab=get(gca,'XTickLabel');
% XTickLab = nan(length(IndLab),1);
% for ii=1:length(XTickLab)
%     XTickLab(ii) = (str2num(IndLab{ii})+1)*10;
% end
% set(gca,'XTickLabel', XTickLab)

SpeedInfo.Semantic=[0.0604 0.0969 0.3815 2.195 11.74 49.56 181.21 552.04 4181.1 10927.3 1520.2 3193.6 18046]
SpeedInfo.Floor=[0.0616 0.0915 0.223 0.974 3.855 16.2 48.32 2500.2 1453.1 9849.6 5628.9 9744.1 6261.2]
SpeedInfo.Ceiling=[0.06466 0.1237 0.736 6.323 35.59 220.61 828.74 1938 5108.5 734.9 20707.5 34805.8 67347.3]
SpeedInfo.AR=[0.06634 0.1632 1.269 5.981 27.93 228.68 783.04 144.09  318.5 2734.4 17373.5 41155.3 0]
WinMax_local = length(SpeedInfo.Semantic);
figure()
plot(1:WinMax_local,SpeedInfo.Semantic, 1:WinMax_local,SpeedInfo.Floor,1:WinMax_local,SpeedInfo.Ceiling,1:WinMax_local,SpeedInfo.AR)
hold on
plot(1:WinMax_local, max([SpeedInfo.Semantic; SpeedInfo.Floor; SpeedInfo.Ceiling;SpeedInfo.AR],[],1), '--r')
legend('Semantic','Floor','Ceiling','AR','Max','Location','NorthWest')