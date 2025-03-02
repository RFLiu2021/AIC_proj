clear all;
close all;

% [numInfo1, txtInfo1, table_all1] = xlsread('../patchSeqMappingTo10X\data\mappingRes/results_Ext.csv');
[numInfo1, txtInfo1, table_all1] = xlsread('../ephysio/data_ephy/Ins_patchseq_ephs.xlsx','Ins_patchseq_ephs');
% [numInfo2, txtInfo2, table_all2] = xlsread('./dataMorph/morphFeatMatrix20240414.csv');
[numInfo2, txtInfo2, table_all2] = xlsread('./dataMorph/morphFeatMatrix.csv');

cellNames1 = txtInfo1(2:end,2);
predictClass = numInfo1(:,1);
cellNames2 = txtInfo2(2:end,1);
morphData = numInfo2(:,2:end);
morphParaNames = txtInfo2(1,3:end);expertCellTypeID = numInfo2(:,1);

% remove the cells without seuencing.
cellNames1 = cellNames1(~isnan(predictClass));
predictClass = predictClass(~isnan(predictClass));
% expertCellTypeID = expertCellTypeID(~isnan(predictClass));

% combine the predict cluster id with morph data.
[cellNames, idx1, idx2] = intersect(cellNames1, cellNames2);
predictClass = predictClass(idx1);
morphData = morphData(idx2,:);
expertCellTypeID = expertCellTypeID(idx2);
cluster_VEN_L5ET =  {16};
cluster_VEN_L56CT = {8};

% assign cells to corresponding clusters.
morphData_VENL= [];  VENL_names={}; tem_cellType  ={};
morphData_VENS = [];  VENS_names={};
morphData_L5ET = []; L5ET_names = {};
morphData_L56CT = []; L56CT_names = {}; count_16 = 0;count_8=0;
for i = 1:length(predictClass)
    switch predictClass(i)
        case cluster_VEN_L5ET
            if expertCellTypeID(i)==4   % cell type is VENL
                morphData_VENL(end+1,:) = morphData(i,:);
                VENL_names{end+1} = cellNames{i};
            else
                morphData_L5ET(end+1,:) = morphData(i,:);
                L5ET_names{end+1} = cellNames{i};
            end
            count_16 = count_16+1;
        case cluster_VEN_L56CT
            if expertCellTypeID(i)==5    % cell type is VENS
                morphData_VENS(end+1,:) = morphData(i,:);
                VENS_names{end+1} = cellNames{i};
            else
                morphData_L56CT(end+1,:) = morphData(i,:);
                L56CT_names{end+1} = cellNames{i};
                
%                 tem_cellType{end+1} = expertCellTypeID{i};
            end
            count_8 = count_8+1;
    end
end 

colorDict = {[255,56,56]/255,[255,41,251]/255,[21,21,211]/255};

examParas = {'axonDist','allBasalDenLen','soma_radius','thickRatio','maxOrder_basal','maxApicalY',};%'apicalLen'};
% plot the high variable parameters (example parameters)..
figure('position',[600,500,600,350]);
markerSize = 50;
s(1)=subplot(2,3,1);
hold on;
paraIdx = find(strcmp(morphParaNames, examParas{1}));
data_tem1 = morphData_VENL(:,paraIdx);%
data_tem2 = morphData_VENS(:,paraIdx);
data_tem3 = morphData_L56CT(:,paraIdx);
nVENL = size(data_tem1,1);
nVENS = size(data_tem2,1);
nL56CT = size(data_tem3,1);
nMax = max([nVENL,nVENS,nL56CT]);
Y = nan(nMax,3);
Y(1:nVENL,1) = reshape(data_tem1,[nVENL,1]);
Y(1:nVENS,2) = reshape(data_tem2,[nVENS,1]);
Y(1:nL56CT,3) = reshape(data_tem3,[nL56CT,1]);
scatter(ones(length(data_tem1),1),data_tem1,markerSize,'marker','.','jitter','on','jitterAmount',0.1,'markerEdgeColor',colorDict{1});
scatter(2*ones(length(data_tem2),1),data_tem2,markerSize,'marker','.','jitter','on','jitterAmount',0.1,'markerEdgeColor',colorDict{2});
scatter(3*ones(length(data_tem3),1),data_tem3,markerSize,'marker','.','jitter','on','jitterAmount',0.1,'markerEdgeColor',colorDict{3});
h=boxplot(Y,'Notch','off','Labels',{'VENL','VENS','L56CT'},'whisker',1,'symbol','','color','k');
set(h,{'linew'},{0.5});
set(gca,'xlim',[0,3.9],'ylim',[-10,150],'xTick', [1,2,3], 'yTick',0:50:150,'box','off');
ylabel(examParas{1});
title(s(1),examParas{1});
hold off;

s(2)=subplot(2,3,2);
hold on;
paraIdx = find(strcmp(morphParaNames, examParas{2}));
data_tem1 = morphData_VENL(:,paraIdx);%
data_tem2 = morphData_VENS(:,paraIdx);%
data_tem3 = morphData_L56CT(:,paraIdx);%
nVENL = size(data_tem1,1);
nVENS = size(data_tem2,1);
nL56CT = size(data_tem3,1);
nMax = max([nVENL,nVENS,nL56CT]);
Y = nan(nMax,3);
Y(1:nVENL,1) = reshape(data_tem1,[nVENL,1]);
Y(1:nVENS,2) = reshape(data_tem2,[nVENS,1]);
Y(1:nL56CT,3) = reshape(data_tem3,[nL56CT,1]);
scatter(ones(length(data_tem1),1),data_tem1,markerSize,'marker','.','jitter','on','jitterAmount',0.1,'markerEdgeColor',colorDict{1});
scatter(2*ones(length(data_tem2),1),data_tem2,markerSize,'marker','.','jitter','on','jitterAmount',0.1,'markerEdgeColor',colorDict{2});
scatter(3*ones(length(data_tem3),1),data_tem3,markerSize,'marker','.','jitter','on','jitterAmount',0.1,'markerEdgeColor',colorDict{3});
h=boxplot(Y,'Notch','off','Labels',{'VENL','VENS','L56CT'},'whisker',1,'symbol','','color','k');
set(h,{'linew'},{0.5});
set(gca,'xlim',[0,4],'ylim',[0,5000],'xTick', [1,2,3], 'yTick',0:1000:5000,'box','off');
ylabel(examParas{2});
title(s(2),examParas{2});
hold off;

s(3)=subplot(2,3,3);
hold on;
paraIdx = find(strcmp(morphParaNames, examParas{3}));
data_tem1 = morphData_VENL(:,paraIdx);%
data_tem2 = morphData_VENS(:,paraIdx);%
data_tem3 = morphData_L56CT(:,paraIdx);%
nVENL = size(data_tem1,1);
nVENS = size(data_tem2,1);
nL56CT = size(data_tem3,1);
nMax = max([nVENL,nVENS,nL56CT]);
Y = nan(nMax,3);
Y(1:nVENL,1) = reshape(data_tem1,[nVENL,1]);
Y(1:nVENS,2) = reshape(data_tem2,[nVENS,1]);
Y(1:nL56CT,3) = reshape(data_tem3,[nL56CT,1]);
scatter(ones(length(data_tem1),1),data_tem1,markerSize,'marker','.','jitter','on','jitterAmount',0.1,'markerEdgeColor',colorDict{1});
scatter(2*ones(length(data_tem2),1),data_tem2,markerSize,'marker','.','jitter','on','jitterAmount',0.1,'markerEdgeColor',colorDict{2});
scatter(3*ones(length(data_tem3),1),data_tem3,markerSize,'marker','.','jitter','on','jitterAmount',0.1,'markerEdgeColor',colorDict{3});
h=boxplot(Y,'Notch','off','Labels',{'VENL','VENS','L56CT'},'whisker',1,'symbol','','color','k');
set(h,{'linew'},{0.5});
set(gca,'xlim',[0,4],'ylim',[5,25],'xTick', [1,2,3], 'yTick',5:5:25,'box','off');
ylabel(examParas{3});
title(s(3),examParas{3});
hold off;

s(4)=subplot(2,3,4);
hold on;
paraIdx = find(strcmp(morphParaNames, examParas{4}));
data_tem1 = morphData_VENL(:,paraIdx);%
data_tem2 = morphData_VENS(:,paraIdx);%
data_tem3 = morphData_L56CT(:,paraIdx);%
nVENL = size(data_tem1,1);
nVENS = size(data_tem2,1);
nL56CT = size(data_tem3,1);
nMax = max([nVENL,nVENS,nL56CT]);
Y = nan(nMax,3);
Y(1:nVENL,1) = reshape(data_tem1,[nVENL,1]);
Y(1:nVENS,2) = reshape(data_tem2,[nVENS,1]);
Y(1:nL56CT,3) = reshape(data_tem3,[nL56CT,1]);
scatter(ones(length(data_tem1),1),data_tem1,markerSize,'marker','.','jitter','on','jitterAmount',0.1,'markerEdgeColor',colorDict{1});
scatter(2*ones(length(data_tem2),1),data_tem2,markerSize,'marker','.','jitter','on','jitterAmount',0.1,'markerEdgeColor',colorDict{2});
scatter(3*ones(length(data_tem3),1),data_tem3,markerSize,'marker','.','jitter','on','jitterAmount',0.1,'markerEdgeColor',colorDict{3});
h=boxplot(Y,'Notch','off','Labels',{'VENL','VENS','L56CT'},'whisker',1,'symbol','','color','k');
set(h,{'linew'},{0.5});
set(gca,'xlim',[0,4],'ylim',[0,4],'xTick', [1,2,3], 'yTick',0:1:4,'box','off');
ylabel(examParas{4});
title(s(4),examParas{4});
hold off;

s(5)=subplot(2,3,5);
hold on;
paraIdx = find(strcmp(morphParaNames, examParas{5}));
data_tem1 = morphData_VENL(:,paraIdx);%
data_tem2 = morphData_VENS(:,paraIdx);%
data_tem3 = morphData_L56CT(:,paraIdx);%
nVENL = size(data_tem1,1);
nVENS = size(data_tem2,1);
nL56CT = size(data_tem3,1);
nMax = max([nVENL,nVENS,nL56CT]);
Y = nan(nMax,3);
Y(1:nVENL,1) = reshape(data_tem1,[nVENL,1]);
Y(1:nVENS,2) = reshape(data_tem2,[nVENS,1]);
Y(1:nL56CT,3) = reshape(data_tem3,[nL56CT,1]);
scatter(ones(length(data_tem1),1),data_tem1,markerSize,'marker','.','jitter','on','jitterAmount',0.1,'markerEdgeColor',colorDict{1});
scatter(2*ones(length(data_tem2),1),data_tem2,markerSize,'marker','.','jitter','on','jitterAmount',0.1,'markerEdgeColor',colorDict{2});
scatter(3*ones(length(data_tem3),1),data_tem3,markerSize,'marker','.','jitter','on','jitterAmount',0.1,'markerEdgeColor',colorDict{3});
h=boxplot(Y,'Notch','off','Labels',{'VENL','VENS','L56CT'},'whisker',1,'symbol','','color','k');
set(h,{'linew'},{0.5});
set(gca,'xlim',[0,4],'ylim',[0,15],'xTick', [1,2,3], 'yTick',0:5:15,'box','off');
ylabel(examParas{5});
title(s(5),examParas{5});
hold off;

s(6)=subplot(2,3,6);
hold on;
paraIdx = find(strcmp(morphParaNames, examParas{6}));
data_tem1 = morphData_VENL(:,paraIdx);%d
data_tem2 = morphData_VENS(:,paraIdx);%
data_tem3 = morphData_L56CT(:,paraIdx);%
nVENL = size(data_tem1,1);
nVENS = size(data_tem2,1);
nL56CT = size(data_tem3,1);
nMax = max([nVENL,nVENS,nL56CT]);
Y = nan(nMax,3);
Y(1:nVENL,1) = reshape(data_tem1,[nVENL,1]);
Y(1:nVENS,2) = reshape(data_tem2,[nVENS,1]);
Y(1:nL56CT,3) = reshape(data_tem3,[nL56CT,1]);
scatter(ones(length(data_tem1),1),data_tem1,markerSize,'marker','.','jitter','on','jitterAmount',0.1,'markerEdgeColor',colorDict{1});
scatter(2*ones(length(data_tem2),1),data_tem2,markerSize,'marker','.','jitter','on','jitterAmount',0.1,'markerEdgeColor',colorDict{2});
scatter(3*ones(length(data_tem3),1),data_tem3,markerSize,'marker','.','jitter','on','jitterAmount',0.1,'markerEdgeColor',colorDict{3});
h=boxplot(Y,'Notch','off','Labels',{'VENL','VENS','L56CT'},'whisker',0.5,'symbol','','color','k');
set(h,{'linew'},{0.5});
set(gca,'xlim',[0,4],'ylim',[0,2000],'xTick', [1,2,3], 'yTick',0:1000:4000,'box','off');
ylabel(examParas{6});
title(s(6),examParas{6});
hold off;

[p,tbl,stats] = kruskalwallis(Y);
c=multcompare(stats)

% figure;
% paraIdx = find(strcmp(morphParaNames, 'maxApicalY'));
% data_tem1 = morphData_VENL(:,paraIdx);%d
% data_tem2 = morphData_VENS(:,paraIdx);%
% data_tem3 = morphData_L56CT(:,paraIdx);%
% nVENL = size(data_tem1,1);
% nVENS = size(data_tem2,1);
% nL56CT = size(data_tem3,1);
% nMax = max([nVENL,nVENS,nL56CT]);
% Y = nan(nMax,3);
% Y(1:nVENL,1) = reshape(data_tem1,[nVENL,1]);
% Y(1:nVENS,2) = reshape(data_tem2,[nVENS,1]);
% Y(1:nL56CT,3) = reshape(data_tem3,[nL56CT,1]);
% scatter(ones(length(data_tem1),1),data_tem1,markerSize,'marker','.','jitter','on','jitterAmount',0.1,'markerEdgeColor',colorDict{1});
% scatter(2*ones(length(data_tem2),1),data_tem2,markerSize,'marker','.','jitter','on','jitterAmount',0.1,'markerEdgeColor',colorDict{2});
% scatter(3*ones(length(data_tem3),1),data_tem3,markerSize,'marker','.','jitter','on','jitterAmount',0.1,'markerEdgeColor',colorDict{3});
% h=boxplot(Y,'Notch','off','Labels',{'VENL','VENS','L56CT'},'whisker',0.5,'symbol','','color','k');
% set(h,{'linew'},{0.5});
% hold off;






%% plot all morphological features comparison.  For testing...
figure;
nFeat = length(morphParaNames);
for i = 1:nFeat-3
    subplot(7,8,i); hold on;
    data_tem1 = morphData_VENL(:,i+2);%d
    data_tem2 = morphData_VENS(:,i+2);%
    data_tem3 = morphData_L56CT(:,i+2);%
    nVENL = size(data_tem1,1);
    nVENS = size(data_tem2,1);
    nL56CT = size(data_tem3,1);
    nMax = max([nVENL,nVENS,nL56CT]);
    Y = nan(nMax,3);
    Y(1:nVENL,1) = reshape(data_tem1,[nVENL,1]);
    Y(1:nVENS,2) = reshape(data_tem2,[nVENS,1]);
    Y(1:nL56CT,3) = reshape(data_tem3,[nL56CT,1]);
    scatter(ones(length(data_tem1),1),data_tem1,markerSize,'marker','.','jitter','on','jitterAmount',0.1,'markerEdgeColor',colorDict{1});
    scatter(2*ones(length(data_tem2),1),data_tem2,markerSize,'marker','.','jitter','on','jitterAmount',0.1,'markerEdgeColor',colorDict{2});
    scatter(3*ones(length(data_tem3),1),data_tem3,markerSize,'marker','.','jitter','on','jitterAmount',0.1,'markerEdgeColor',colorDict{3});
    h=boxplot(Y,'Notch','off','Labels',{'VENL','VENS','L56CT'},'whisker',2,'symbol','','color','k');
    set(h,{'linew'},{0.5});
    % set(gca,'xlim',[0,4],'ylim',[0,6],'xTick', [1,2,3,4,5], 'yTick',0:2:10,'box','off');
    ylabel(morphParaNames{i+2});
    title(gca,morphParaNames{i+2});
    hold off;

    [p_tem,h_tem,stat_tem]=ranksum(Y(:,1),Y(:,2));
    [ave,sd,se] = imean(Y,1);
    VENL_ave(i) = ave(1);
    VENS_ave(i) = ave(2);
    VENL_sd(i) = sd(1);
    VENS_sd(i) = sd(2);
    p(i) = p_tem;
    ransumv(i) = stat_tem.ranksum;
end

paraNames = morphParaNames(3:end);
ave_Table = [VENL_ave',VENS_ave'];
sd_Table = [VENL_sd',VENS_sd'];
dataTable = [VENL_ave',VENL_sd',VENS_ave',VENS_sd',p',ransumv'];


% %% ==============helper=============================
% %% ===============calculate the statistic values===============
% paraName = 'nBasalBranch';
% paraIdx = find(strcmp(morphParaNames,paraName));
% data_tem1 = morphData_VENL(:,paraIdx);%d
% data_tem2 = morphData_VENS(:,paraIdx);%
% [ave1,sd1,se1] = imean(data_tem1);
% [ave2,sd2,se2] = imean(data_tem2);
% [p,h]=ranksum(data_tem1,data_tem2);
% fprintf('VEN-L: %6.2f±%6.2f.\n',ave1,sd1);
% fprintf('VEN-S: %6.2f±%6.2f.\n',ave2,sd2);
% fprintf('p = : %6.5f.\n',p);