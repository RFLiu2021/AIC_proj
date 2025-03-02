clear all;

cellCountsFile = './dataMorph/cellCounts_AIC_vs_V1.xlsx';
[num,txt,~] = xlsread( cellCountsFile, 'Sheet1');

% =====calculate the cell density====
% the numbers in th excel file are number of cells in each ROI region (100um X 100um) of 
% 30-um-thick slices. Therefore, the average numbers of each collumn are
% the estimated cell numbers (N) in a 100um X 100um X 30um space. The number of
% cells in a 1 cubic mm should be N*(100/30)*1000.
nCell_AIC_L23 = num(:,1)*3.3*1000;
nCell_AIC_L5 = num(:,2)*3.3*1000;
nCell_V1_L23 = num(:,5)*3.3*1000;
nCell_V1_L5 = num(:,6)*3.3*1000;

Y = num(:,[1,5,2,6]).*3.3.*1000;
colorLib = {'r','b'};
markerSize = 30;

figure('position',[700,500,500,250]);
subplot(1,2,1);hold on;
scatter(ones(size(nCell_AIC_L23)),nCell_AIC_L23,markerSize,'marker','o','jitter','on','jitterAmount',0.1,'markerEdgeColor',colorLib{1});
scatter(2*ones(size(nCell_V1_L23)),nCell_V1_L23,markerSize,'marker','o','jitter','on','jitterAmount',0.1,'markerEdgeColor',colorLib{2});
scatter(3*ones(size(nCell_AIC_L5)),nCell_AIC_L5,markerSize,'marker','^','jitter','on','jitterAmount',0.1,'markerEdgeColor',colorLib{1});
scatter(4*ones(size(nCell_V1_L5)),nCell_V1_L5,markerSize,'marker','^','jitter','on','jitterAmount',0.1,'markerEdgeColor',colorLib{2});
boxplot(Y,'Notch','off','Labels',{'AIC_L2/3','V1_L2/3','AIC_L5','V1_L5'},'whisker',1,'symbol','','colors','k');
set(gca,'xlim',[0.5,4.5],'ylim',[0,500000])
ylabel('# of cells per mm3');

% output the results
[aveN_AIC_L23,sdN_AIC_L23] = imean(nCell_AIC_L23);
[aveN_V1_L23,sdN_V1_L23] = imean(nCell_V1_L23);
[aveN_AIC_L5,sdN_AIC_L5] = imean(nCell_AIC_L5);
[aveN_V1_L5,sdN_V1_L5] = imean(nCell_V1_L5);
fprintf('---------------------Neuron density----------------------\n');
fprintf('Average cell density in AIC Layer 2/3: %d±% per mm^3.\n',aveN_AIC_L23, sdN_AIC_L23);
fprintf('Average cell density in V1 Layer 2/3: %d±%d mm^3.\n',aveN_V1_L23, sdN_V1_L23);
fprintf('Average cell density in AIC Layer 5: %d±%d mm^3.\n',aveN_AIC_L5, sdN_AIC_L5);
fprintf('Average cell density in V1 Layer 5: %d±%d mm^3.\n',aveN_V1_L5, sdN_V1_L5);

% ======calculate the cell diameter ============
cell_AIC_L23 = [num(4:38,19);num(4:38,20);num(4:38,21)];
cell_AIC_L5 = [num(4:38,23);num(4:38,24);num(4:38,25)];
cell_V1_L23 = [num(4:38,27);num(4:38,28);num(4:38,29)];
cell_V1_L5 = [num(4:38,31);num(4:38,32);num(4:38,33)];

Y = nan(35,4);
Y(1:length(cell_AIC_L23),1) = cell_AIC_L23;
Y(1:length(cell_V1_L23),2) = cell_V1_L23;
Y(1:length(cell_AIC_L5),3) = cell_AIC_L5;
Y(1:length(cell_V1_L5),4) = cell_V1_L5;

colorLib = {'b','r'};
markerSize = 30;

subplot(1,2,2);hold on;
scatter(ones(size(cell_AIC_L23)),cell_AIC_L23,markerSize,'marker','o','jitter','on','jitterAmount',0.1,'markerEdgeColor',colorLib{1});
scatter(2*ones(size(cell_V1_L23)),cell_V1_L23,markerSize,'marker','o','jitter','on','jitterAmount',0.1,'markerEdgeColor',colorLib{2});
scatter(3*ones(size(cell_AIC_L5)),cell_AIC_L5,markerSize,'marker','^','jitter','on','jitterAmount',0.1,'markerEdgeColor',colorLib{1});
scatter(4*ones(size(cell_V1_L5)),cell_V1_L5,markerSize,'marker','^','jitter','on','jitterAmount',0.1,'markerEdgeColor',colorLib{2});
boxplot(Y,'Notch','off','Labels',{'AIC_L2/3','V1_L2/3','AIC_L5','V1_L5'},'whisker',1,'symbol','','colors','k');
set(gca,'xlim',[0.5,4.5],'ylim',[0,30])
ylabel('Diameter of cell body (um)');

% output the results
[aveD_AIC_L23,sdD_AIC_L23] = imean(cell_AIC_L23);
[aveD_V1_L23,sdD_V1_L23] = imean(cell_AIC_L5);
[aveD_AIC_L5,sdD_AIC_L5] = imean(cell_V1_L23);
[aveD_V1_L5,sdD_V1_L5] = imean(cell_V1_L5);
fprintf('--------------------Soma diameter--------------------------\n');
fprintf('Average oma diameter in AIC Layer 2/3: %6.2f±%6.2f um.\n',aveD_AIC_L23,sdD_AIC_L23);
fprintf('Average oma diameter in V1 Layer 2/3: %6.2f±%6.2f um.\n',aveD_V1_L23,sdD_V1_L23);
fprintf('Average oma diameter in AIC Layer 5: %6.2f±%6.2f um.\n',aveD_AIC_L5,sdD_AIC_L5);
fprintf('Average soma diameter in V1 Layer 5: %6.2f±%6.2f um.\n',aveD_V1_L5,sdD_V1_L5);
fprintf('--------------------End---------------------------\n');


