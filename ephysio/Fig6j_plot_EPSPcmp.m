clear all;
close all;


%% =============plot example eEPSPs comparison ========================
load('data_ephy\eEPSPs_examples.mat');
freqSample = PC.stimEPSP_WM(1).freqSample;
data_c1_stimInWM_Min = PC.stimEPSP_WM(1).data;
data_c1_stimInWM_2Min = PC.stimEPSP_WM(2).data;
data_c1_stimInWM_4Min = PC.stimEPSP_WM(3).data;
data_c2_stimInWM_Min = VEN.stimEPSP_WM(1).data;
data_c2_stimInWM_2Min = VEN.stimEPSP_WM(2).data;
data_c2_stimInWM_4Min = VEN.stimEPSP_WM(3).data;

figure('position',[700,600,300,500]); 

% -----------------PC VS.VEN (stim in WM, low stim)----------------------------
data = mean(data_c1_stimInWM_Min,2);
numTrial = size(data,2);
refData = data(1:4000,:);
refY = mean( mean(refData) );
meanY = mean(refData,1);
for i = 1:length(meanY)
    dif = meanY(i)-refY;
    data(:,i) = data(:,i)-dif;
    data(:,i) = smooth(data(:,i),5);
end
aveData_c1_Min = mean(data,2);

data = mean(data_c2_stimInWM_Min,2);
numTrial = size(data,2);
refData = data(1:4000,:);
refY = mean( mean(refData) );
meanY = mean(refData,1);
for i = 1:length(meanY)
    dif = meanY(i)-refY;
    data(:,i) = data(:,i)-dif;
    data(:,i) = smooth(data(:,i),5);
end
aveData_c2_Min = mean(data,2);

% normalize traces aligning the rest potentials between c2 and c1.
sampleData(:,1) = aveData_c1_Min(1:4000);
sampleData(:,2) = aveData_c2_Min(1:4000);
refY = mean( mean(sampleData) );
meanSample_c1 = mean(sampleData(:,1),1);
aveData_c1_Min = aveData_c1_Min -(meanSample_c1-refY);
meanSample_c2 = mean(sampleData(:,2),1);
aveData_c2_Min = aveData_c2_Min -(meanSample_c2-refY);
% cut the first EPSPs.
d_c1 = aveData_c1_Min(3000:10000); % cut the first evoked EPSP. Sample frequency:20000 points/s
d_c2 = aveData_c2_Min(3000:10000);
% resample data
numPoint = length(d_c1);
idx = 1:20:numPoint;
d_c1 = d_c1(idx);
d_c2 = d_c2(idx);
t = ( 1:length(d_c1) ).*20/freqSample;
%plot the traces
subplot(3,1,1);
plot(t,d_c1,'b-');
hold on;plot(t,d_c2,'r-');
set(gca,'ylim',[-0.095,-0.045]);
title('VEN Vs PC (0.05 mA)');

% -----------------PC VS.VEN (stim in WM, middle stim)----------------------------
data = mean(data_c1_stimInWM_2Min,2);
numTrial = size(data,2);
refData = data(1:4000,:);
refY = mean( mean(refData) );
meanY = mean(refData,1);
for i = 1:length(meanY)
    dif = meanY(i)-refY;
    data(:,i) = data(:,i)-dif;
    data(:,i) = smooth(data(:,i),5);
end
aveData_c1_mid = mean(data,2);

data = mean(data_c2_stimInWM_2Min,2);
numTrial = size(data,2);
refData = data(1:4000,:);
refY = mean( mean(refData) );
meanY = mean(refData,1);
for i = 1:length(meanY)
    dif = meanY(i)-refY;
    data(:,i) = data(:,i)-dif;
    data(:,i) = smooth(data(:,i),5);
end
aveData_c2_mid = mean(data,2);

% normalize traces aligning the rest potentials between c2 and c3.
sampleData(:,1) = aveData_c1_mid(1:4000);
sampleData(:,2) = aveData_c2_mid(1:4000);
refY = mean( mean(sampleData) );
meanSample_c2 = mean(sampleData(:,1),1);
aveData_c1_mid = aveData_c1_mid -(meanSample_c2-refY);
meanSample_c3 = mean(sampleData(:,2),1);
aveData_c2_mid = aveData_c2_mid -(meanSample_c3-refY);
% cut the first EPSPs.
d_c1 = aveData_c1_mid(3000:10000); % cut the first evoked EPSP. Sample frequency:20000 points/s
d_c2 = aveData_c2_mid(3000:10000);
% resample data
numPoint = length(d_c1);
idx = 1:20:numPoint;
d_c1 = d_c1(idx);
d_c2 = d_c2(idx);
t = ( 1:length(d_c1) ).*20/freqSample;
%plot the traces
subplot(3,1,2);
plot(t,d_c1,'b-');
hold on;plot(t,d_c2,'r-');
set(gca,'ylim',[-0.095,-0.045]);
title('VEN Vs PC (0.1 mA)');
ylabel("Membrane potential (Voltage)");

% -----------------PC VS.VEN (stim in WM, high stim)----------------------------
data = mean(data_c1_stimInWM_4Min,2);
numTrial = size(data,2);
refData = data(1:4000,:);
refY = mean( mean(refData) );
meanY = mean(refData,1);
for i = 1:length(meanY)
    dif = meanY(i)-refY;
    data(:,i) = data(:,i)-dif;
    data(:,i) = smooth(data(:,i),5);
end
aveData_c1_max = mean(data,2);

data = mean(data_c2_stimInWM_4Min,2);
numTrial = size(data,2);
refData = data(1:4000,:);
refY = mean( mean(refData) );
meanY = mean(refData,1);
for i = 1:length(meanY)
    dif = meanY(i)-refY;
    data(:,i) = data(:,i)-dif;
    data(:,i) = smooth(data(:,i),5);
end
aveData_c2_max = mean(data,2);

% normalize traces aligning the rest potentials between c2 and c3.
sampleData(:,1) = aveData_c1_max(1:4000);
sampleData(:,2) = aveData_c2_max(1:4000);
refY = mean( mean(sampleData) );
meanSample_c2 = mean(sampleData(:,1),1);
aveData_c1_max = aveData_c1_max -(meanSample_c2-refY);
meanSample_c3 = mean(sampleData(:,2),1);
aveData_c2_max = aveData_c2_max -(meanSample_c3-refY);
% cut the first EPSPs.
d_c1 = aveData_c1_max(3000:10000); % cut the first evoked EPSP. Sample frequency:20000 points/s
d_c2 = aveData_c2_max(3000:10000);
% resample data
numPoint = length(d_c1);
idx = 1:20:numPoint;
d_c1 = d_c1(idx);
d_c2 = d_c2(idx);
t = ( 1:length(d_c1) ).*20/freqSample;
%plot the traces
subplot(3,1,3);
plot(t,d_c1,'b-');
hold on;plot(t,d_c2,'r-');
h_ax = gca;
h_ax.XLabel.String = 'Time (s)';
set(h_ax,'ylim',[-0.095,-0.045]);
title('VEN Vs PC (0.2 mA)');











%% ================== plot population comparison ====================
dataFile = './data_ephy/sumData_StimWM.xlsx';
[num1,txt,raw] = xlsread(dataFile,'re-orgnized_data');


% -----------sort data and plot figure (All VENs vs. PCs)-------------
aveAmp_PC = num1(81:83,7).*1000;
seAmp_PC = num1(85:87,7).*1000;
aveAmp_VEN = num1(81:83,8).*1000;
seAmp_VEN = num1(85:87,8).*1000;

I_in = [0.5,1,2];
figure('position',[700,600,750,250]); 
subplot(1,3,1); hold on;
errorbar(I_in,aveAmp_PC,seAmp_PC,'blue','marker','o','LineWidth',1);
errorbar(I_in,aveAmp_VEN,seAmp_VEN,'red','marker','o','LineWidth',1);
set(gca,'xlim',[0,2.5],'ylim',[0,20],'xTick',[0,0.5,1,2],'yTick',[0:10:40],'box','off');
xlabel('Input current (mA)');
ylabel('EPSP amp (mV)');
legend('PC','VEN');
title("PC vs VEN");


% -------------sort data and plot figure (VEN-L vs. PCs)---------------
aveAmp_PC = num1(39:41,13).*1000;
seAmp_PC = num1(43:45,13).*1000;
aveAmp_VENL = num1(39:41,14).*1000;
seAmp_VENL = num1(43:45,14).*1000;

subplot(1,3,2); hold on;
errorbar(I_in,aveAmp_PC,seAmp_PC,'blue','marker','o','LineWidth',1);
errorbar(I_in,aveAmp_VENL,seAmp_VENL,'red','marker','o','LineWidth',1);
set(gca,'xlim',[0,2.5],'ylim',[0,22],'xTick',[0,0.5,1,2],'yTick',[0:10:40],'box','off');
xlabel('Input current (mA)');
ylabel('EPSP amp (mV)');
legend('PC','VEN-L');
title("PC vs VEN-L");


% ---------------sort data and plot figure (VEN-S vs. PCs)---------------
aveAmp_PC = num1(87:89,13).*1000;
seAmp_PC = num1(92:94,13).*1000;
aveAmp_VENS = num1(87:89,14).*1000;
seAmp_VENS = num1(92:94,14).*1000;

subplot(1,3,3); hold on
errorbar(I_in,aveAmp_PC,seAmp_PC,'blue','marker','o','LineWidth',1);
errorbar(I_in,aveAmp_VENS,seAmp_VENS,'color','#FFA500','marker','o','LineWidth',1);
set(gca,'xlim',[0,2.5],'ylim',[0,20],'xTick',[0,0.5,1,2],'yTick',[0:10:40],'box','off');
xlabel('Input current (mA)');
ylabel('EPSP amp (mV)');
legend('PC','VEN-S');
title("PC vs VEN-S");




