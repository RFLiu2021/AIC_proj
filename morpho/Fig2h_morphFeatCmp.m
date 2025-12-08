
clear all;
close all;

% define the cell groups.
cluster_L23IT = {1,15,5};
cluster_IT_rorb = {2,0,7,13,19};
cluster_L5IT = {3,6,9,11,14,18};
cluster_VEN_L5ET =  {16};
cluster_VEN_L56CT = {8};
cluster_L6PC = {4,12,17};
colorTable = {'#008000','#98A14E','#4666A6','#E50012','#00CED1','#6FBC1E'};

% load the morph feature data
dataFile = './dataMorph/morphFeatMatrix.xlsx';
[dataMat,txt] = xlsread(dataFile);
cellNames = txt(2:end,1);
expertCellType = txt(2:end,3);
predictClusterIDs = dataMat(:,3);
morphParaNames = txt(1,5:end);
dataMat = dataMat(:,4:end);

% initialize the cell group data
morphData_L23IT = [];  cellNames_L23IT = [];
morphData_ITrorb = [];  cellNames_ITrorb = [];
morphData_L5IT = [];  cellNames_L5IT = [];
morphData_VEN = [];  cellNames_VEN = [];
morphData_L5ET = [];  cellNames_L5ET = [];
morphData_L6PC = [];  cellNames_L6PC = [];
% morphData_L56CT = [];  cellNames_L56CT = [];

for i = 1:length(predictClusterIDs)
    switch predictClusterIDs(i)
        case cluster_L23IT
            morphData_L23IT(end+1,:) = dataMat(i,:);
            cellNames_L23IT{end+1} = cellNames{i};
        case cluster_IT_rorb
            morphData_ITrorb(end+1,:) = dataMat(i,:);
            cellNames_ITrorb{end+1} = cellNames{i};
        case cluster_L5IT
            morphData_L5IT(end+1,:) = dataMat(i,:);
            cellNames_L5IT{end+1} = cellNames{i};
        case cluster_VEN_L5ET
            if contains(expertCellType{i},'VEN') | contains(expertCellType{i},'ven')
                morphData_VEN(end+1,:) = dataMat(i,:);
                cellNames_VEN{end+1} = cellNames{i};
            elseif contains(expertCellType{i},'PC') | contains(expertCellType{i},'pc')
                morphData_L5ET(end+1,:) = dataMat(i,:);
                cellNames_L5ET{end+1} = cellNames{i};
            end
        case cluster_VEN_L56CT
            if contains(expertCellType{i},'VEN') | contains(expertCellType{i},'ven')
                morphData_VEN(end+1,:) = dataMat(i,:);
                cellNames_VEN{end+1} = cellNames{i};
            end
        case cluster_L6PC
            morphData_L6PC(end+1,:) = dataMat(i,:);
            cellNames_L6PC{end+1} = cellNames{i};
    end
end 



examPara= {'apicalLen','totalVolume','lengthFractionAboveSoma_basal',...
    'centroid_apicalY', 'centroid_allBasalDenY', 'somaSize'};

% plot the high variable parameters (example parameters)..
figure('position',[500,300,300,500]);
markerSize = 30;
s(1)=subplot(2,1,1);hold on;
col_idx = find(strcmp(morphParaNames, examPara{4}));
nL23IT = size(morphData_L23IT(:,col_idx),1);
nITrorb = size(morphData_ITrorb(:,col_idx),1);
nL5IT = size(morphData_L5IT(:,col_idx),1);
nVEN = size(morphData_VEN(:,col_idx),1);
nL5ET = size(morphData_L5ET(:,col_idx),1);
nL6PC = size(morphData_L6PC(:,col_idx),1);
nMax = max([nL23IT,nL5IT,nVEN,nL6PC,nL5ET]);
Y = nan(nMax,6);
Y(1:nL23IT,1) = reshape(morphData_L23IT(:,col_idx),[nL23IT,1]);
Y(1:nITrorb,2) = reshape(morphData_ITrorb(:,col_idx),[nITrorb,1]);
Y(1:nL5IT,3) = reshape(morphData_L5IT(:,col_idx),[nL5IT,1]);
Y(1:nVEN,4) = reshape(morphData_VEN(:,col_idx),[nVEN,1]);
Y(1:nL5ET,5) = reshape(morphData_L5ET(:,col_idx),[nL5ET,1]);
Y(1:nL6PC,6) = reshape(morphData_L6PC(:,col_idx),[nL6PC,1]);
scatter(ones(size(morphData_L23IT,1),1),morphData_L23IT(:,col_idx),markerSize,'marker','.','jitter','on','jitterAmount',0.1,'markerEdgeColor',colorTable{1});
scatter(2*ones(size(morphData_ITrorb,1),1),morphData_ITrorb(:,col_idx),markerSize,'marker','.','jitter','on','jitterAmount',0.1,'markerEdgeColor',colorTable{2});
scatter(3*ones(size(morphData_L5IT,1),1),morphData_L5IT(:,col_idx),markerSize,'marker','.','jitter','on','jitterAmount',0.1,'markerEdgeColor',colorTable{3});
scatter(4*ones(size(morphData_VEN,1),1),morphData_VEN(:,col_idx),markerSize,'marker','.','jitter','on','jitterAmount',0.1,'markerEdgeColor',colorTable{4});
scatter(5*ones(size(morphData_L5ET,1),1),morphData_L5ET(:,col_idx),markerSize,'marker','.','jitter','on','jitterAmount',0.1,'markerEdgeColor',colorTable{5});
scatter(6*ones(size(morphData_L6PC,1),1),morphData_L6PC(:,col_idx),markerSize,'marker','.','jitter','on','jitterAmount',0.1,'markerEdgeColor',colorTable{6});
h=boxplot(Y,'Notch','off','Labels',{'L23IT','IT_RORB','L5IT','VEN','L5ET','L6PC'},'whisker',2,'symbol','','colors','k');
set(h,{'linew'},{0.5});
set(gca,'xlim',[0,6.9],'xTick', 1:7,'ylim',[-200,1500], 'yTick',0:500:1600,'box','off')%,'yscale','log');
ylabel('centroid_apicalY (um)');
title(s(1),'centroid_apicalY');
hold off;


s(2)=subplot(2,1,2);hold on;
col_idx = find(strcmp(morphParaNames, examPara{5}));
nL23IT = size(morphData_L23IT(:,col_idx),1);
nITrorb = size(morphData_ITrorb(:,col_idx),1);
nL5IT = size(morphData_L5IT(:,col_idx),1);
nVEN = size(morphData_VEN(:,col_idx),1);
nL5ET = size(morphData_L5ET(:,col_idx),1);
nL6PC = size(morphData_L6PC(:,col_idx),1);
nMax = max([nL23IT,nL5IT,nVEN,nL6PC,nL5ET]);
Y = nan(nMax,6);
Y(1:nL23IT,1) = reshape(morphData_L23IT(:,col_idx),[nL23IT,1]);
Y(1:nITrorb,2) = reshape(morphData_ITrorb(:,col_idx),[nITrorb,1]);
Y(1:nL5IT,3) = reshape(morphData_L5IT(:,col_idx),[nL5IT,1]);
Y(1:nVEN,4) = reshape(morphData_VEN(:,col_idx),[nVEN,1]);
Y(1:nL5ET,5) = reshape(morphData_L5ET(:,col_idx),[nL5ET,1]);
Y(1:nL6PC,6) = reshape(morphData_L6PC(:,col_idx),[nL6PC,1]);
scatter(ones(size(morphData_L23IT,1),1),morphData_L23IT(:,col_idx),markerSize,'marker','.','jitter','on','jitterAmount',0.1,'markerEdgeColor',colorTable{1});
scatter(2*ones(size(morphData_ITrorb,1),1),morphData_ITrorb(:,col_idx),markerSize,'marker','.','jitter','on','jitterAmount',0.1,'markerEdgeColor',colorTable{2});
scatter(3*ones(size(morphData_L5IT,1),1),morphData_L5IT(:,col_idx),markerSize,'marker','.','jitter','on','jitterAmount',0.1,'markerEdgeColor',colorTable{3});
scatter(4*ones(size(morphData_VEN,1),1),morphData_VEN(:,col_idx),markerSize,'marker','.','jitter','on','jitterAmount',0.1,'markerEdgeColor',colorTable{4});
scatter(5*ones(size(morphData_L5ET,1),1),morphData_L5ET(:,col_idx),markerSize,'marker','.','jitter','on','jitterAmount',0.1,'markerEdgeColor',colorTable{5});
scatter(6*ones(size(morphData_L6PC,1),1),morphData_L6PC(:,col_idx),markerSize,'marker','.','jitter','on','jitterAmount',0.1,'markerEdgeColor',colorTable{6});
h=boxplot(Y,'Notch','off','Labels',{'L23IT','IT_RORB','L5IT','VEN','L5ET','L6PC'},'whisker',2,'symbol','','colors','k');
set(h,{'linew'},{0.5});
set(gca,'xlim',[0,6.9],'xTick', 1:10,'ylim',[-250,200], 'yTick',-300:100:400,'box','off');
ylabel('centroid_BasalY (um)');
title(s(2),'centroid_BasalY');
hold off;







% suppFig13_plotAllmorphFeatCmp;




