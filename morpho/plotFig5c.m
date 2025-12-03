clear all; 
close all;

%
[n,t,r] = xlsread('./dataMorph/AIScmp_PCvsVEN_Fig5c.xlsx','Sheet1');

AIS_L(:,1) = n(:,2);   % VEN
AIS_L(:,2) = n(:,1);   % PC
AIS_D(:,1) = n(:,5);   % VEN
AIS_D(:,2) = n(:,4);   % PC

figure('Position',[700,500,300,600]);
markerSize = 50; colorTable = {[21,21,211]/255; [255,56,56]/ 255};
subplot(2,1,1); hold on;
scatter(ones(size(AIS_L,1),1),AIS_L(:,1),markerSize,'marker','.','jitter','on','jitterAmount',0.1,'markerEdgeColor',colorTable{1});
scatter(2*ones(size(AIS_L,1),1),AIS_L(:,2),markerSize,'marker','.','jitter','on','jitterAmount',0.1,'markerEdgeColor',colorTable{2});
h=boxplot(AIS_L,'Notch','off','Labels',{'PC','VEN'},'whisker',1,'symbol','','colors','k');
set(h,{'linew'},{1});
set(gca,'xlim',[0,2.9],'ylim',[-0,80],'xTick', 1:10, 'yTick', 0:20:80, 'box','off');
ylabel("AIS L (um)")

subplot(2,1,2); hold on;
scatter(ones(size(AIS_D,1),1),AIS_D(:,1),markerSize,'marker','.','jitter','on','jitterAmount',0.1,'markerEdgeColor',colorTable{1});
scatter(2*ones(size(AIS_D,1),1),AIS_D(:,2),markerSize,'marker','.','jitter','on','jitterAmount',0.1,'markerEdgeColor',colorTable{2});
h=boxplot(AIS_D,'Notch','off','Labels',{'PC','VEN'},'whisker',1,'symbol','','colors','k');
set(h,{'linew'},{1});
set(gca,'xlim',[0,2.9],'ylim',[-0,100],'xTick', 1:10, 'yTick', 0:25:80, 'box','off');
ylabel("AIS D (um)")









