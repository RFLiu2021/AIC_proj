clear all;
close all;

% load peak indices data from an excel file
[num,txt,raw] = xlsread('./data_ephy/PeakIdx_Fig5e.xlsx','Sheet1');

markerSize = 150;
colorTable = {'#008000','#98A14E','#4666A6','#E50012','#00CED1','#6FBC1E'};

% sorting data according to cell types.
L23IT_peakIdx = num(1:30,1);
ITrorb_peakIdx = num(1:22,4);
L5IT_peakIdx = num(1:122,7);
VEN_peakIdx = num(1:29,10);
L5ET_peakIdx = num(1:17,13);
L6PC_peakIdx = num(1:13,16);

% build the matrix for plotting.
mat_peakIdx = nan(122,6);
mat_peakIdx(1:30,1) = L23IT_peakIdx;
mat_peakIdx(1:22,2) = ITrorb_peakIdx;
mat_peakIdx(1:122,3) = L5IT_peakIdx;
mat_peakIdx(1:29,4) = VEN_peakIdx;
mat_peakIdx(1:17,5) = L5ET_peakIdx;
mat_peakIdx(1:13,6) = L6PC_peakIdx;

% plotting...
figure('Position',[800,600,600,300]);
hold on;
scatter(ones(size(L23IT_peakIdx,1),1),L23IT_peakIdx,markerSize,'marker','.','jitter','on','jitterAmount',0.1,'markerEdgeColor',colorTable{1});
scatter(2*ones(size(ITrorb_peakIdx,1),1),ITrorb_peakIdx,markerSize,'marker','.','jitter','on','jitterAmount',0.1,'markerEdgeColor',colorTable{2});
scatter(3*ones(size(L5IT_peakIdx,1),1),L5IT_peakIdx,markerSize,'marker','.','jitter','on','jitterAmount',0.1,'markerEdgeColor',colorTable{3});
scatter(4*ones(size(VEN_peakIdx,1),1),VEN_peakIdx,markerSize,'marker','.','jitter','on','jitterAmount',0.1,'markerEdgeColor',colorTable{4});
scatter(5*ones(size(L5ET_peakIdx,1),1),L5ET_peakIdx,markerSize,'marker','.','jitter','on','jitterAmount',0.1,'markerEdgeColor',colorTable{5});
scatter(6*ones(size(L6PC_peakIdx,1),1),L6PC_peakIdx,markerSize,'marker','.','jitter','on','jitterAmount',0.1,'markerEdgeColor',colorTable{6});
h=boxplot(mat_peakIdx,'Notch','off','Labels',{'L23IT','IT_RORB','L5IT','VEN','L5ET','L6PC'},'whisker',1,'symbol','','colors','k');
set(h,{'linew'},{0.5});
set(gca,'xlim',[0,6.9],'ylim',[-0.5,1.0],'xTick', 1:10, 'yTick', -0.5:0.5:1.0, 'box','off');%,'yscale','log');
ylabel("Peak index");
title("Peak index comparison");
hold off;

