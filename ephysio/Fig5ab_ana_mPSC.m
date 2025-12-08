clear all;
close all;

% import data from the new analysis result file. (X. Tie did)
xlsxFile = './data_ephy/miniEPSCs_Fig7a-b.xlsx';
nums =xlsread(xlsxFile,'mEPSC_sum'); 
mEPSC_VEN_freq = nums(1:8,1);
mEPSC_PC_freq = nums(1:8,5);
mEPSC_VEN_amp = nums(1:8,2);
mEPSC_PC_amp = nums(1:8,6);
mEPSC_freq(1:8,1) = mEPSC_PC_freq;
mEPSC_freq(1:8,2) = mEPSC_VEN_freq;
mEPSC_amp(1:8,1) = mEPSC_PC_amp;  
mEPSC_amp(1:8,2) = mEPSC_VEN_amp;
mEPSC = [mEPSC_freq,mEPSC_amp];


%=================================================================
% plot the figure in box style
figure('Position',[600,500,500,300]);
subplot(1,2,1)
h1 = boxplot(mEPSC_freq,'Labels',{'PC_freq','VEN_freq'},'symbol','','Notch','off','Whisker',2);
set(gca,'xTickLabel',{'Freq','Amp'},'ylim',[0,3])
hold on; 
scatter(ones(numel(mEPSC(:,1)),1),mEPSC(:,1),50,'jitter','on','jitterAmount',0.2 ,'markerFaceColor',[21,21,211]/255 ,'markerEdgeColor','k' );
scatter(2*ones(numel(mEPSC(:,2)),1),mEPSC(:,2),50,'jitter','on','jitterAmount',0.2 ,'markerFaceColor',[255,56,56]/255 ,'markerEdgeColor','k' );
% plot(1,imean(mEPSC_freq(:,1)),'bo');
% plot(2,imean(mEPSC_freq(:,2)),'ro');
box off;
set(gca,'xlim',[0.5,2.5],'ylim',[0,4],'YTick',[0,1,2,3,4,5]);
title('mEPSC Freq comprision');
ylabel("Frequency (events/s)");
[p,h] = ranksum(mEPSC_PC_freq(:,1),mEPSC_VEN_freq(:,1));
fprintf('Comparision of mEPSC frequencies bewteen PCs and VENs, and p = %6.4f.\n',p);

subplot(1,2,2)
h2 = boxplot(mEPSC_amp,'Labels',{'PC_amp','VEN_amp'},'symbol','','Notch','off','Whisker',2);
set(gca,'xTickLabel',{'Freq','Amp'})
hold on; 
scatter(ones(numel(mEPSC(:,3)),1),mEPSC(:,3),50,'jitter','on','jitterAmount',0.2 ,'markerFaceColor',[21,21,211]/255 ,'markerEdgeColor','k' );
scatter(2*ones(numel(mEPSC(:,4)),1),mEPSC(:,4),50,'jitter','on','jitterAmount',0.2 ,'markerFaceColor',[255,56,56]/255 ,'markerEdgeColor','k' );
% plot(1,imean(mEPSC_amp(:,1)),'bo');
% plot(2,imean(mEPSC_amp(:,2)),'ro');
box off;
set(gca,'xlim',[0.5,2.5],'ylim',[0,0.15],'YTick',[0,0.05,0.10,0.15,0.20],'YTickLabel',{'0','5','10','15','20'}); % Here, the unit of vertical axis is pA, the values were multiplicated 1E+10 when transfer from raw data to txt format, Please see 'exportWave2txt.m'.
title('mEPSC Amp comprision');
ylabel("Amplitude (pA)");
[p,h] = ranksum(mEPSC_PC_amp(:,1),mEPSC_VEN_amp(:,1));
fprintf('Comparision of mEPSC amplitude bewteen PCs and VENs, and p = %6.4f.\n',p);


%% plot the miniEPSCs of example cells
nums =xlsread(xlsxFile,'mIPSC_exmples'); 
exmPC = nums(1:5000,2);
exmVEN = nums(1:5000,6);

figure;
subplot(2,1,1);
plot(exmPC,'color',[0.0824,0.0824,0.8275]);
hold on; plot([3500,4500],[-0.4,-0.4],'k','LineWidth',2); % plot scale bar.
set(gca,'ylim',[-1.2,-0.2]);
text(500,-0.4,'miniEPSC of Example PC');
subplot(2,1,2);
plot(exmVEN,'color',[1,0.2196,0.2196]);hold on; plot([4500,4500],[-1,-0.8],'k','LineWidth',2)
set(gca,'ylim',[-1.7,-0.7]);
text(500,-0.9,'miniEPSC of Example VEN');
xlabel('Time (ms)');
ylabel('Voltage (mV)');



%%  =========================Analysis of mIPSPs =========================

% ------------- plot the miniIPSCs of example cells----------------
nums =xlsread(xlsxFile,'mIPSC_exmples'); 
t_PC = nums(1:5000,1)/1000;
exmPC = nums(1:5000,2);
t_VEN = nums(1:5000,4)/1000;
exmVEN = nums(1:5000,5);

figure;
subplot(2,1,1);hold on; 
plot(t_PC,exmPC,'color',[0.0824,0.0824,0.8275]);
plot([21.5,22.5],[-50,-50],'k','LineWidth',2); % plot scale bar.
plot([22.5,22.5],[-50,-20],'k','LineWidth',2);
set(gca,'ylim',[-60,0]);
xlabel("Time (s)");
ylabel("mIPSC (pA)");
title('mIPSC of example PC');
subplot(2,1,2);hold on;
plot(t_VEN,exmVEN,'color',[1,0.2196,0.2196]);
xlabel("Time (s)");
ylabel("mIPSC (pA)");
set(gca,'ylim',[-100,-40]);
title('mIPSC of Example VEN');



% ---------------- plot the boxplot for Frequency---------------------
[nums,txt,~] =xlsread(xlsxFile,'mIPSC_sum'); 

cellType = txt(2:end,2);
freq_mIPSC = nums(:,3);
amp_mIPSC = nums(:,4);
freq_PC = freq_mIPSC( strcmp(cellType,'PC') );
freq_VEN = freq_mIPSC( strcmp(cellType,'VEN') );
amp_PC = amp_mIPSC( strcmp(cellType,'PC') );
amp_VEN = amp_mIPSC( strcmp(cellType,'VEN') );

nPC = length(freq_PC);
nVEN = length(freq_VEN);

N = max(nPC,nVEN);
mIPSC_freq = nan(N,2); 
mIPSC_amp = mIPSC_freq;

mIPSC_freq(1:nPC,1) = freq_PC; 
mIPSC_freq(1:nVEN,2) = freq_VEN;
mIPSC_amp(1:nPC,1) = amp_PC; 
mIPSC_amp(1:nVEN,2) = amp_VEN;






%=================================================================
% plot the figure in box style
markersize = 50;
figure('Position',[600,500,500,300]);

subplot(1,2,1);
h1 = boxplot(mIPSC_freq,'Labels',{'PC_freq','VEN_freq'},'symbol','','Notch','off','Whisker',2);
set(gca,'xTickLabel',{'PC','VEN'},'ylim',[0,10])
hold on; 
scatter(ones(numel(mIPSC_freq(:,1)),1),mIPSC_freq(:,1),markersize,'jitter','on','jitterAmount',0.2 ,'markerEdgeColor',[21,21,211]/255 ,'markerFaceColor','none' );
scatter(2*ones(numel(mIPSC_freq(:,2)),1),mIPSC_freq(:,2),markersize,'jitter','on','jitterAmount',0.2 ,'markerEdgeColor',[255,56,56]/255 ,'markerFaceColor','none' );
box off;
set(gca,'xlim',[0.5,2.5],'ylim',[0,10],'YTick',0:2.5:10);
title('mIPSC Freq');
ylabel("Freauency (events/s)");
[p,h] = ranksum(mIPSC_freq(:,1),mIPSC_freq(:,2));
fprintf('Comparision of mEPSC frequencies bewteen PCs and VENs, and p = %6.4f.\n',p);
hold off;

subplot(1,2,2);
h2 = boxplot(mIPSC_amp,'Labels',{'PC_freq','VEN_freq'},'symbol','','Notch','off','Whisker',2);
set(gca,'xTickLabel',{'PC','VEN'},'ylim',[0,10])
hold on; 
scatter(ones(numel(mIPSC_amp(:,1)),1),mIPSC_amp(:,1),markersize,'jitter','on','jitterAmount',0.2 ,'markerEdgeColor',[21,21,211]/255 ,'markerFaceColor','none' );
scatter(2*ones(numel(mIPSC_amp(:,2)),1),mIPSC_amp(:,2),markersize,'jitter','on','jitterAmount',0.2 ,'markerEdgeColor',[255,56,56]/255 ,'markerFaceColor','none' );
box off;
set(gca,'xlim',[0.5,2.5],'ylim',[0,50],'YTick',0:12.5:50);
title('mIPSC Amp');
ylabel("Amplitude (pA)");
[p,h] = ranksum(mIPSC_amp(:,1),mIPSC_amp(:,2));
fprintf('Comparision of mEPSC frequencies bewteen PCs and VENs, and p = %6.4f.\n',p);
hold off;















