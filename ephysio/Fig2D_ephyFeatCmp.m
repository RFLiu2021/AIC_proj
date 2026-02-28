clear all;
close all;

[numInfo, txtInfo, table_all] = xlsread('./data_ephy/Ins_patchseq_ephs.xlsx','Sheet1');

predictClass = numInfo(:,1);
cellNames = txtInfo(2:end,1);
parameterNames = txtInfo(1,4:end);
expertCellType = txtInfo(2:end,2);

ephsData = numInfo(:,2:end);

% remove the cells without seuencing or failed quality control.
ephsData = ephsData(~isnan(predictClass),:);
predictClass = predictClass(~isnan(predictClass));

cluster_L23IT = {1,15,5};
cluster_IT_rorb = {2,0,7,13,19};
cluster_L5IT = {3,6,9,11,14,18};
cluster_VEN_L5ET =  {16};
cluster_VEN_L56CT = {8};
cluster_L6PC = {4,12,17};

% ephsData_L23IT = numInfo(any(cluster_L23IT,predictClass),:);
% ephsData_L5ITT = numInfo(any(cluster_L5IT,predictClass),:);
% ephsData_ET = numInfo(any(cluster_LET,predictClass),:);
% ephsData_L5ET = numInfo(any(cluster_L5ET,predictClass),:);
% ephsData_L6PC = numInfo(any(cluster_L6PC,predictClass),:);

ephsData_L23IT = [];L23IT_names ={};
ephsData_ITrorb = [];ITrorb_names ={};
ephsData_L5IT = []; L5IT_names={};
ephsData_VEN = [];  VEN_names={};
ephsData_L5ET = []; L5ET_names = {};
ephsData_L6PC = []; L6PC_names = {};
for i = 1:length(predictClass)
    switch predictClass(i)
        case cluster_L23IT
            ephsData_L23IT(end+1,:) = ephsData(i,:);
            L23IT_names{end+1} = cellNames{i};
        case cluster_IT_rorb
            ephsData_ITrorb(end+1,:) = ephsData(i,:);
            ITrorb_names{end+1} = cellNames{i};    
        case cluster_L5IT
            ephsData_L5IT(end+1,:) = ephsData(i,:);
            L5IT_names{end+1} = cellNames{i};
        case cluster_VEN_L5ET
            if contains(expertCellType{i},'VEN') | contains(expertCellType{i},'ven')
                ephsData_VEN(end+1,:) = ephsData(i,:);
                VEN_names{end+1} = cellNames{i};
            elseif contains(expertCellType{i},'PC') | contains(expertCellType{i},'pc')
                ephsData_L5ET(end+1,:) = ephsData(i,:);
                L5ET_names{end+1} = cellNames{i};
            end
        case cluster_VEN_L56CT
            if contains(expertCellType{i},'VEN') | contains(expertCellType{i},'ven')
                ephsData_VEN(end+1,:) = ephsData(i,:);
                VEN_names{end+1} = cellNames{i};
            end
        case cluster_L6PC
            ephsData_L6PC(end+1,:) = ephsData(i,:);
            L6PC_names{end+1} = cellNames{i};
    end
end 

colorTable = {'#008000','#98A14E','#4666A6','#E50012','#00CED1','#6FBC1E'};

% plot the high variable parameters (example parameters)..
figure('Position',[700,500,1000,300]);
markerSize = 50;
% -------------------------------Sag ratio--------------------------
s(1)=subplot(1,3,1);hold on;
data_tem1 = ephsData_L23IT(:,1);data_tem1 = data_tem1(data_tem1>1); % exclude the abnoral traces.
data_tem2 = ephsData_ITrorb(:,1);data_tem2 = data_tem2(data_tem2>1);
data_tem3 = ephsData_L5IT(:,1);data_tem3 = data_tem3(data_tem3>1);
data_tem4 = ephsData_VEN(:,1);data_tem4 = data_tem4(data_tem4>1);
data_tem5 = ephsData_L5ET(:,1);data_tem5 = data_tem5(data_tem5>1);
data_tem6 = ephsData_L6PC(:,1);data_tem6 = data_tem6(data_tem6>1);
nL23IT = size(data_tem1,1);
nITrorb = size(data_tem2,1);
nL5IT = size(data_tem3,1);
nVEN = size(data_tem4,1);
nL5ET = size(data_tem5,1);
nL6PC = size(data_tem6,1);
nMax = max([nL23IT,nITrorb,nL5IT,nVEN,nL5ET,nL6PC]);
Y = nan(nMax,6);
Y(1:nL23IT,1) = reshape(data_tem1,[nL23IT,1]);
Y(1:nITrorb,2) = reshape(data_tem2,[nITrorb,1]);
Y(1:nL5IT,3) = reshape(data_tem3,[nL5IT,1]);
Y(1:nVEN,4) = reshape(data_tem4,[nVEN,1]);
Y(1:nL5ET,5) = reshape(data_tem5,[nL5ET,1]);
Y(1:nL6PC,6) = reshape(data_tem6,[nL6PC,1]);
scatter(ones(length(data_tem1),1),data_tem1,markerSize,'marker','.','jitter','on','jitterAmount',0.1,'markerEdgeColor',colorTable{1});
scatter(2*ones(length(data_tem2),1),data_tem2,markerSize,'marker','.','jitter','on','jitterAmount',0.1,'markerEdgeColor',colorTable{2});
scatter(3*ones(length(data_tem3),1),data_tem3,markerSize,'marker','.','jitter','on','jitterAmount',0.1,'markerEdgeColor',colorTable{3});
scatter(4*ones(length(data_tem4),1),data_tem4,markerSize,'marker','.','jitter','on','jitterAmount',0.1,'markerEdgeColor',colorTable{4});
scatter(5*ones(length(data_tem5),1),data_tem5,markerSize,'marker','.','jitter','on','jitterAmount',0.1,'markerEdgeColor',colorTable{5});
scatter(6*ones(length(data_tem6),1),data_tem6,markerSize,'marker','.','jitter','on','jitterAmount',0.1,'markerEdgeColor',colorTable{6});
h=boxplot(Y,'Notch','off','Labels',{'L23IT','ITrorb','L5IT','VEN','L5ET','L6PC'},'whisker',1,'symbol','','colors','k');
set(h,{'linew'},{0.5});
set(gca,'xlim',[0,7],'ylim',[0,120],'xTick', [1,2,3,4,5,6], 'yTick',0:40:120,'box','off');
ylabel('Tau (ms)');
title(s(1),'Membrane time constant');
hold off;
% ------------------ Input resistance --------------------
s(2)=subplot(1,3,2);hold on;
data_tem1 = ephsData_L23IT(:,3);
data_tem2 = ephsData_ITrorb(:,3);
data_tem3 = ephsData_L5IT(:,3);
data_tem4 = ephsData_VEN(:,3);
data_tem5 = ephsData_L5ET(:,3);
data_tem6 = ephsData_L6PC(:,3);
nL23IT = size(data_tem1,1);
nITrorb = size(data_tem2,1);
nL5IT = size(data_tem3,1);
nVEN = size(data_tem4,1);
nL5ET = size(data_tem5,1);
nL6PC = size(data_tem6,1);
nMax = max([nL23IT,nITrorb,nL5IT,nVEN,nL5ET,nL6PC]);
Y = nan(nMax,6);
Y(1:nL23IT,1) = reshape(data_tem1,[nL23IT,1]);
Y(1:nITrorb,2) = reshape(data_tem2,[nITrorb,1]);
Y(1:nL5IT,3) = reshape(data_tem3,[nL5IT,1]);
Y(1:nVEN,4) = reshape(data_tem4,[nVEN,1]);
Y(1:nL5ET,5) = reshape(data_tem5,[nL5ET,1]);
Y(1:nL6PC,6) = reshape(data_tem6,[nL6PC,1]);
scatter(ones(length(data_tem1),1),data_tem1,markerSize,'marker','.','jitter','on','jitterAmount',0.1,'markerEdgeColor',colorTable{1});
scatter(2*ones(length(data_tem2),1),data_tem2,markerSize,'marker','.','jitter','on','jitterAmount',0.1,'markerEdgeColor',colorTable{2});
scatter(3*ones(length(data_tem3),1),data_tem3,markerSize,'marker','.','jitter','on','jitterAmount',0.1,'markerEdgeColor',colorTable{3});
scatter(4*ones(length(data_tem4),1),data_tem4,markerSize,'marker','.','jitter','on','jitterAmount',0.1,'markerEdgeColor',colorTable{4});
scatter(5*ones(length(data_tem5),1),data_tem5,markerSize,'marker','.','jitter','on','jitterAmount',0.1,'markerEdgeColor',colorTable{5});
scatter(6*ones(length(data_tem6),1),data_tem6,markerSize,'marker','.','jitter','on','jitterAmount',0.1,'markerEdgeColor',colorTable{6});
h=boxplot(Y,'Notch','off','Labels',{'L23IT','ITrorb','L5IT','VEN','L5ET','L6PC'},'whisker',1,'symbol','','colors','k');
set(h,{'linew'},{0.5});
set(gca,'xlim',[0,7],'ylim',[0,1000],'xTick', [1,2,3,4,5,6], 'yTick',0:250:1000,'box','off');
ylabel('Rm (MOhm)');
title(s(2),'Input resistance');
hold off;

% ------------------ Threshold --------------------
s(3)=subplot(1,3,3);
hold on;
data_tem1 = ephsData_L23IT(:,12); data_tem1 = data_tem1(data_tem1<5); data_tem1 = data_tem1(data_tem1>1); % exclude the unresonable.
data_tem2 = ephsData_ITrorb(:,12); data_tem2 = data_tem2(data_tem2<5); data_tem2 = data_tem2(data_tem2>1);
data_tem3 = ephsData_L5IT(:,12); data_tem3 = data_tem3(data_tem3<5); data_tem3 = data_tem3(data_tem3>1);
data_tem4 = ephsData_VEN(:,12); data_tem4 = data_tem4(data_tem4<5); data_tem4 = data_tem4(data_tem4>1);
data_tem5 = ephsData_L5ET(:,12); data_tem5 = data_tem5(data_tem5<5); data_tem5 = data_tem5(data_tem5>1);
data_tem6 = ephsData_L6PC(:,12); data_tem6 = data_tem6(data_tem6<5); data_tem6 = data_tem6(data_tem6>1);
nL23IT = size(data_tem1,1);
nITrorb = size(data_tem2,1);
nL5IT = size(data_tem3,1);
nVEN = size(data_tem4,1);
nL5ET = size(data_tem5,1);
nL6PC = size(data_tem6,1);
nMax = max([nL23IT,nITrorb,nL5IT,nVEN,nL5ET,nL6PC]);
Y = nan(nMax,6);
Y(1:nL23IT,1) = reshape(data_tem1,[nL23IT,1]);
Y(1:nITrorb,2) = reshape(data_tem2,[nITrorb,1]);
Y(1:nL5IT,3) = reshape(data_tem3,[nL5IT,1]);
Y(1:nVEN,4) = reshape(data_tem4,[nVEN,1]);
Y(1:nL5ET,5) = reshape(data_tem5,[nL5ET,1]);
Y(1:nL6PC,6) = reshape(data_tem6,[nL6PC,1]);
scatter(ones(length(data_tem1),1),data_tem1,markerSize,'marker','.','jitter','on','jitterAmount',0.1,'markerEdgeColor',colorTable{1});
scatter(2*ones(length(data_tem2),1),data_tem2,markerSize,'marker','.','jitter','on','jitterAmount',0.1,'markerEdgeColor',colorTable{2});
scatter(3*ones(length(data_tem3),1),data_tem3,markerSize,'marker','.','jitter','on','jitterAmount',0.1,'markerEdgeColor',colorTable{3});
scatter(4*ones(length(data_tem4),1),data_tem4,markerSize,'marker','.','jitter','on','jitterAmount',0.1,'markerEdgeColor',colorTable{4});
scatter(5*ones(length(data_tem5),1),data_tem5,markerSize,'marker','.','jitter','on','jitterAmount',0.1,'markerEdgeColor',colorTable{5});
scatter(6*ones(length(data_tem6),1),data_tem6,markerSize,'marker','.','jitter','on','jitterAmount',0.1,'markerEdgeColor',colorTable{6});
h=boxplot(Y,'Notch','off','Labels',{'L23IT','ITrorb','L5IT','VEN','L5ET','L6PC'},'whisker',1,'symbol','','colors','k');

set(h,{'linew'},{0.5});
set(gca,'xlim',[0,7],'ylim',[0,8],'xTick', [1,2,3,4,5,6], 'yTick',0:2:6,'box','off');
ylabel('AP width (ms)');
title(s(3),'AP width');
hold off;


% ========================plot all ephysio features====================
suppFig10a_plot_all_ephysFeatures();










