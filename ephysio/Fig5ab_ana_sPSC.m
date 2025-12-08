clear

% ====================================New analysis================================
xlsxFile = './data_ephy/sPSCs_Fig7ab.xlsx';
nums =xlsread(xlsxFile,'summary'); 

EPSC_PC = nums(1:11,1:5);
EPSC_VEN = nums(1:8,10:11);
IPSC_PC = nums(19:26,1:5);
IPSC_VEN = nums(19:25,10:11);

f1 = figure();
% ========================Frequencies comparision======================
% EPSC Frequencies comparision.
EPSCs = nan(20, 2);
EPSCs(1:length(EPSC_PC),1) = EPSC_PC(:,1);
EPSCs(1:length(EPSC_VEN),2) = EPSC_VEN(:,1);
subplot(1,2,1);%figure; 
hold on; 
h = boxplot(EPSCs,'Labels',{'PCs','VENs'},'symbol','');
scatter(ones(size(EPSC_PC,1),1),EPSC_PC(:,1),50,'jitter','on','markerEdgeColor','k','markerFaceColor',[21,21,211]/255);
scatter(2*ones(size(EPSC_VEN,1),1),EPSC_VEN(:,1),50,'jitter','on','markerEdgeColor','k','markerFaceColor',[255,56,56]/255);
% plot(1,imean(EPSCs(:,1)),'bo');
% plot(2,imean(EPSCs(:,2)),'ro');
box off;
set(gca,'xlim',[0.5,2.5],'ylim',[0,10],'YTick',[0,2,4,6,8,10]);
title('EPSC Freq comprision');
[p1,h,stat1] = ranksum(EPSC_PC(:,1),EPSC_VEN(:,1));
fprintf('Comparision of EPSC frequencies bewteen PCs and VENs, and p = %6.4f.\n',p1);

% IPSC Frequencies comparision.
IPSCs = nan(20, 2);
IPSCs(1:length(IPSC_PC),1) = IPSC_PC(:,1);
IPSCs(1:length(IPSC_VEN),2) = IPSC_VEN(:,1);
subplot(1,2,2);%figure; 
hold on; 
h = boxplot(IPSCs,'Labels',{'PCs','VENs'},'symbol','');
scatter(ones(size(IPSC_PC,1),1),IPSC_PC(:,1),50,'jitter','on','markerEdgeColor','k','markerFaceColor',[21,21,211]/255);
scatter(2*ones(size(IPSC_VEN,1),1),IPSC_VEN(:,1),50,'jitter','on','markerEdgeColor','k','markerFaceColor',[255,56,56]/255);
% plot(1,imean(IPSCs(:,1)),'bo');
% plot(2,imean(IPSCs(:,2)),'ro');
box off;
set(gca,'xlim',[0.5,2.5],'ylim',[0,10],'YTick',[0,2,4,6,8,10]);
title('IPSC Freq comprision');
[p2,h,stat2] = ranksum(IPSC_PC(:,1),IPSC_VEN(:,1));
fprintf('Comparision of IPSC frequencies bewteen PCs and VENs, and p = %6.4f.\n',p2);





% ========================Amplitude comparision======================
% EPSC amplitude comparision between PCs and VENs.
EPSCs = nan(20, 2);
EPSCs(1:length(EPSC_PC),1) = EPSC_PC(:,2);
EPSCs(1:length(EPSC_VEN),2) = EPSC_VEN(:,2);
figure;
subplot(1,2,1);%
h = boxplot(EPSCs,'Labels',{'PCs','VENs'},'symbol','');
hold on; 
scatter(ones(size(EPSCs,1),1),EPSCs(:,1),50,'jitter','on','markerEdgeColor','k','markerFaceColor',[21,21,211]/255);
scatter(2*ones(size(EPSCs,1),1),EPSCs(:,2),50,'jitter','on','markerEdgeColor','k','markerFaceColor',[255,56,56]/255);
% plot(1,imean(EPSCs(:,1)),'bo');
% plot(2,imean(EPSCs(:,2)),'ro');
box off;
set(gca,'xlim',[0.5,2.5],'ylim',[0,0.8],'YTick',[0,0.2,0.4,0.6,0.8],'YTickLabel',{'0','20','40','60','80'});  % Here, the unit of vertical axis is pA, the values were multiplicated 1E+10 when transfer from raw data to txt format, Please see 'exportWave2txt.m'.
title('EPSC Amp comprision');
[p,h] = ranksum(EPSC_PC(:,2),EPSC_VEN(:,2));
fprintf('Comparision of EPSC amplitudes bewteen PCs and VENs, and p = %6.4f.\n',p);

% IPSC amplitude comparision between PCs and VENs.
IPSCs = nan(20, 2);
IPSCs(1:length(IPSC_PC),1) = IPSC_PC(:,2);
IPSCs(1:length(IPSC_VEN),2) = IPSC_VEN(:,2);
subplot(1,2,2);%figure;
h = boxplot(IPSCs,'Labels',{'PCs','VENs'},'symbol','');
hold on; 
scatter(ones(size(IPSCs,1),1),IPSCs(:,1),50,'jitter','on','markerEdgeColor','k','markerFaceColor',[21,21,211]/255);
scatter(2*ones(size(IPSCs,1),1),IPSCs(:,2),50,'jitter','on','markerEdgeColor','k','markerFaceColor',[255,56,56]/255);
% plot(1,imean(IPSCs(:,1)),'bo');
% plot(2,imean(IPSCs(:,2)),'ro');
box off;
set(gca,'xlim',[0.5,2.5],'ylim',[0,2],'YTick',[0,0.5,1.0,1.5,2],'YTickLabel',{'0','50','100','150','200'});  % Here, the unit of vertical axis is pA, the values were multiplicated 1E+10 when transfer from raw data to txt format, Please see 'exportWave2txt.m'.
title('IPSC Amp comprision');
[p,h] = ranksum(IPSC_PC(:,2),IPSC_VEN(:,2));
fprintf('Comparision of IPSC amplitudes bewteen PCs and VENs, and p = %6.4f.\n',p);


%% plot example cell data
% % -----------------Exmaple sEPSC traces------------------------------
nums =xlsread(xlsxFile,'exmaples_sEPSP');
t_PC = nums(:,1)/1000;  % Unit: s
d_PC = nums(:,2);
t_VEN = nums(:,5)/1000;  % Unit: s
d_VEN = nums(:,6);

figure;subplot(4,1,1); %PC
plot(t_PC,d_PC,'b');
set(gca,'ylim',[-20,20])
title('EPSC PC');
 
subplot(4,1,2);% VEN
plot(t_VEN,d_VEN,'r');
set(gca,'ylim',[-120,-80]);
hold on;plot([106,107],[-100,-100],'k','LineWidth',2); % plot Scale Bar.
hold on;plot([107,107],[-120,-100],'k','LineWidth',2); % plot Scale Bar.
hold off;
title('EPSC VEN');



% % -----------------Exmaple sIPSC traces------------------------------
% ---Example IPSC of a PC----
nums =xlsread(xlsxFile,'exmaples_sIPSP');
t_PC = nums(:,1);  % Unit: s
d_PC = nums(:,2);
t_VEN = nums(:,9);  % Unit: s
d_VEN = nums(:,10);

% plot data.
subplot(4,1,3);%figure;
plot( t_PC,d_PC ,'b');
set(gca,'ylim',[-150,-50])
title('IPSC PC');
hold on;plot([81,82],[-100,-100],'k','LineWidth',2); % plot Scale Bar.
hold on;plot([82,82],[-150,-100],'k','LineWidth',2); % plot Scale Bar

% ----Example IPSC of a VEN----
subplot(4,1,4);%figure;
plot(t_VEN,d_VEN,'r');
set(gca,'ylim',[-260,-160])
title('IPSC VEN');
hold off;