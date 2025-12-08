
clear all;

spikeFile = './data_ephy/exampleSpikes_Fig6ac.xlsx';
[n,t,r] = xlsread(spikeFile,'Fig.4a');

h=figure('units','centimeters','papersize',[17,12],'Position',[25,10,15,20]);

t_VEN = n(:,1);
venFirstAP_data = n(:,2);
t_PC = n(:,5);
pyrFirstAP_data = n(:,6);
subplot(3,1,1);
plot(t_VEN,venFirstAP_data,'r');
hold on; plot(t_PC,pyrFirstAP_data,'b');
set(gca,'xlim',[0,15]);% example data sampleFreq = 50000, 1000/50000 = 0.02sec, same as heatmap creatted by 'plotAPwaveHeatmap.m'

% plot the heatmat figures.
AP_DataMat1_pc = n(1:239,11:end);
AP_DataMat1_ven = n(245:511,11:end);
subplot(3,1,2);
heatmap(AP_DataMat1_pc,'GridVisible','off','ColorLimits',[0,120],'XLabel','','YLabel','','Colormap', summer);
title('Spikes heatmaps - PCs');
subplot(3,1,3);
heatmap(AP_DataMat1_ven,'GridVisible','off','ColorLimits',[0,120],'XLabel','','YLabel','','Colormap', summer);%parula);
title('Spikes heatmaps - VENs');





% ---------------plot the normalized first action potentials and compare between VEN and PC---------
[n,t,r] = xlsread(spikeFile,'Fig.4b');

h=figure('units','centimeters','papersize',[17,12],'Position',[25,10,15,20]);

t_VEN = n(:,1);
venFirstAP_data = n(:,2);
t_PC = n(:,5);
pyrFirstAP_data = n(:,6);
subplot(3,1,1);
plot(t_PC,pyrFirstAP_data,'b');
hold on;
plot(t_VEN,venFirstAP_data,'r');
set(gca,'xlim',[0,11],'ylim',[0,1]);

% plot the heatmat figures.
AP_rise_pc = n(1:237,11:261);
AP_decay_pc = n(242:478,11:511);
AP_rise_ven = n(484:749,11:261);
AP_decay_ven = n(756:1021,11:511);

subplot('Position',[0.1,0.4,0.2,0.25])
hh = heatmap(AP_rise_pc,'GridVisible','off','ColorLimits',[0,1],'XLabel','','YLabel','','Colormap', summer);
subplot('Position',[0.45,0.4,0.5,0.25])
heatmap(AP_decay_pc,'GridVisible','off','ColorLimits',[0,1],'XLabel','','YLabel','','Colormap', summer);

subplot('Position',[0.1,0.05,0.2,0.25])
hh = heatmap(AP_rise_ven,'GridVisible','off','ColorLimits',[0,1],'XLabel','','YLabel','','Colormap', summer);
subplot('Position',[0.45,0.05,0.5,0.25])
heatmap(AP_decay_ven,'GridVisible','off','ColorLimits',[0,1],'XLabel','','YLabel','','Colormap', summer);




% -------------plot the suprathreshold stim evoked spikes comparasion-------- 
[n,t,r] = xlsread(spikeFile,'Fig.4c');

h=figure('units','centimeters','papersize',[17,12],'Position',[25,10,15,20]);

t_VEN = n(:,1);
venFirstAP_data = n(:,2);
venLastAP_data = n(:,3);
t_PC = t_VEN;
pyrFirstAP_data = n(:,4);
pyrLastAP_data = n(:,5);

subplot(5,1,1);
plot(venFirstAP_data,'r');
hold on; 
plot(pyrFirstAP_data,'b');
plot(venLastAP_data,'r--');
plot(pyrLastAP_data,'b--');
set(gca,'xlim',[0,700]);  % example data sampleFreq = 50000, 1300/50000 = 0.026sec, same as heatmap creatted by 'plotAPwaveHeatmap.m'


AP_DataMat2_pc = n(2:214,9:end);
AP_DataMat2_ven = n(221:455,9:end);

subplot(5,1,2);
hh = heatmap(AP_DataMat1_pc,'GridVisible','off','ColorLimits',[0,120],'XLabel','','YLabel','','Colormap', summer);
subplot(5,1,3);
hh = heatmap(AP_DataMat2_pc,'GridVisible','off','ColorLimits',[0,120],'XLabel','','YLabel','','Colormap', summer);
subplot(5,1,4);
hh = heatmap(AP_DataMat1_ven,'GridVisible','off','ColorLimits',[0,120],'XLabel','','YLabel','','Colormap', summer);
subplot(5,1,5);
hh = heatmap(AP_DataMat2_ven,'GridVisible','off','ColorLimits',[0,120],'XLabel','','YLabel','','Colormap', summer);





