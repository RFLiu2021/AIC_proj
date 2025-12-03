clear all;
close all;

[numInfo, txtInfo, table_all] = xlsread('./data_ephy/Ins_patchseq_ephs.xlsx','Sheet1');

predictClass = numInfo(:,1);
cellNames = txtInfo(2:end,1);
% parameterNames = txtInfo(1,4:end);
expertCellType = txtInfo(2:end,2);
% ephsData = numInfo(:,2:end);

% remove the cells without seuencing or failed quality control.
% ephsData = ephsData(~isnan(predictClass),:);
predictClass = predictClass(~isnan(predictClass));

cluster_L23IT = {1,15,5};
cluster_IT_rorb = {2,0,7,13,19};
cluster_L5IT = {3,6,9,11,14,18};
cluster_VEN_L5ET =  {16};
cluster_VEN_L56CT = {8};
cluster_L6PC = {4,12,17};

ephsData_L23IT = [];L23IT_names ={};
ephsData_ITrorb = [];ITrorb_names ={};
ephsData_L5IT = []; L5IT_names={};
ephsData_VEN = [];  VEN_names={};
ephsData_L5ET = []; L5ET_names = {};
ephsData_L6PC = []; L6PC_names = {};
for i = 1:length(predictClass)
    switch predictClass(i)
        case cluster_L23IT
            % ephsData_L23IT(end+1,:) = ephsData(i,:);
            L23IT_names{end+1} = cellNames{i};
        case cluster_IT_rorb
            % ephsData_ITrorb(end+1,:) = ephsData(i,:);
            ITrorb_names{end+1} = cellNames{i};    
        case cluster_L5IT
            % ephsData_L5IT(end+1,:) = ephsData(i,:);
            L5IT_names{end+1} = cellNames{i};
        case cluster_VEN_L5ET
            if contains(expertCellType{i},'VEN') | contains(expertCellType{i},'ven')
                % ephsData_VEN(end+1,:) = ephsData(i,:);
                VEN_names{end+1} = cellNames{i};
            elseif contains(expertCellType{i},'PC') | contains(expertCellType{i},'pc')
                % ephsData_L5ET(end+1,:) = ephsData(i,:);
                L5ET_names{end+1} = cellNames{i};
            end
        case cluster_VEN_L56CT
            if contains(expertCellType{i},'VEN') | contains(expertCellType{i},'ven')
                % ephsData_VEN(end+1,:) = ephsData(i,:);
                VEN_names{end+1} = cellNames{i};
            end
        case cluster_L6PC
            % ephsData_L6PC(end+1,:) = ephsData(i,:);
            L6PC_names{end+1} = cellNames{i};
    end
end 

cellTypes = {"L23IT","IT Rorb","L5IT","VEN","L5ET","L6PC"};
colorTable = {'#008000','#98A14E','#4666A6','#E50012','#00CED1','#6FBC1E'};










%% ==================load cell data from .mat files================
PC_names = [L23IT_names,ITrorb_names,L5ET_names,L5IT_names,L6PC_names];
VEN_names;

nPC = length(PC_names);
nVEN = length(VEN_names);

nSpikeMat_ven = nan(nVEN,50);
for i = 1:nVEN
    VEN_fileName_t = fullfile('./data_ephy/seqCellData',VEN_names{i});
    load(VEN_fileName_t);
    curr0_idx = m_FP.curr_index_0;
    counter=1;
    for j = curr0_idx:length(Attr(1).APproperties)
        if (~isempty(Attr(1).APproperties(j).numSpikes) && counter <=50)
            nSpikeMat_ven(i,counter) = Attr(1).APproperties(j).numSpikes;
        end
        counter = counter+1;
    end
end
% % quality control.
% maxN = max(nSpikeMat_ven,[],2);
% idx=maxN>=5;
% nSpikeMat_ven = nSpikeMat_ven(idx,:);

nL23IT = length(L23IT_names);
nSpikeMat_L23IT = nan(nL23IT,50);
for i = 1:nL23IT
    L23IT_fileName_t = fullfile('./data_ephy/seqCellData',L23IT_names{i});
    load(L23IT_fileName_t);
    curr0_idx = m_FP.curr_index_0;
    counter=1;
    for j = curr0_idx:length(Attr(1).APproperties)
        if (~isempty(Attr(1).APproperties(j).numSpikes) && counter <=50)
            nSpikeMat_L23IT(i,counter) = Attr(1).APproperties(j).numSpikes;
        end
        counter = counter+1;
    end
end
% % quality control.
% maxN = max(nSpikeMat_L23IT,[],2);
% idx=find(maxN>=5);
% nSpikeMat_L23IT = nSpikeMat_L23IT(idx,:);


nITrorb = length(ITrorb_names);
nSpikeMat_ITrorb = nan(nITrorb,50);
for i = 1:nITrorb
    ITrorb_fileName_t = fullfile('./data_ephy/seqCellData',ITrorb_names{i});
    load(ITrorb_fileName_t);
    curr0_idx = m_FP.curr_index_0;
    counter=1;
    for j = curr0_idx:length(Attr(1).APproperties)
        if (~isempty(Attr(1).APproperties(j).numSpikes) && counter <=50)
            nSpikeMat_ITrorb(i,counter) = Attr(1).APproperties(j).numSpikes;
        end
        counter = counter+1;
    end
end
% % quality control.
% maxN = max(nSpikeMat_ITrorb,[],2);
% idx=find(maxN>=5);
% nSpikeMat_ITrorb = nSpikeMat_ITrorb(idx,:);


nL5IT = length(L5IT_names);
nSpikeMat_L5IT = nan(nL5IT,50);
for i = 1:nL5IT
    L5IT_fileName_t = fullfile('./data_ephy/seqCellData',L5IT_names{i});
    load(L5IT_fileName_t);
    curr0_idx = m_FP.curr_index_0;
    counter=1;
    for j = curr0_idx:length(Attr(1).APproperties)
        if (~isempty(Attr(1).APproperties(j).numSpikes) && counter <=50)
            nSpikeMat_L5IT(i,counter) = Attr(1).APproperties(j).numSpikes;
        end
        counter = counter+1;
    end
end
% % quality control.
% maxN = max(nSpikeMat_L5IT,[],2);
% idx=find(maxN>=5);
% nSpikeMat_L5IT = nSpikeMat_L5IT(idx,:);


nL5ET = length(L5ET_names);
nSpikeMat_L5ET = nan(nL5ET,50);
for i = 1:nL5ET
    L5ET_fileName_t = fullfile('./data_ephy/seqCellData',L5ET_names{i});
    load(L5ET_fileName_t);
    curr0_idx = m_FP.curr_index_0;
    counter=1;
    for j = curr0_idx:length(Attr(1).APproperties)
        if (~isempty(Attr(1).APproperties(j).numSpikes) && counter <=50)
            nSpikeMat_L5ET(i,counter) = Attr(1).APproperties(j).numSpikes;
        end
        counter = counter+1;
    end
end
% % quality control.
% maxN = max(nSpikeMat_L5ET,[],2);
% idx=find(maxN>=5);
% nSpikeMat_L5ET = nSpikeMat_L5ET(idx,:);

nL6PC = length(L6PC_names);
nSpikeMat_L6PC = nan(nL6PC,50);
for i = 1:nL6PC
    L6PC_fileName_t = fullfile('./data_ephy/seqCellData',L6PC_names{i});
    load(L6PC_fileName_t);
    curr0_idx = m_FP.curr_index_0;
    counter=1;
    for j = curr0_idx:length(Attr(1).APproperties)
        if (~isempty(Attr(1).APproperties(j).numSpikes) && counter <=50)
            nSpikeMat_L6PC(i,counter) = Attr(1).APproperties(j).numSpikes;
        end
        counter = counter+1;
    end
end
% % quality control.
% maxN = max(nSpikeMat_L6PC,[],2);
% idx=find(maxN>=5);
% nSpikeMat_L6PC = nSpikeMat_L6PC(idx,:);


%% ============== plot the firing rates vs input current intensity===========
[FRs_ven,~,se_ven] = imean(nSpikeMat_ven,1);  % mean number of spikes.
[FRs_L23IT,~,se_L23IT] = imean(nSpikeMat_L23IT,1);
[FRs_ITrorb,~,se_ITrorb] = imean(nSpikeMat_ITrorb,1);
[FRs_L5IT,~,se_L5IT] = imean(nSpikeMat_L5IT,1);
[FRs_L5ET,~,se_L5ET] = imean(nSpikeMat_L5ET,1);
[FRs_L6PC,~,se_L6PC] = imean(nSpikeMat_L6PC,1);

markersize = 10;
figure;hold on;
I_in = [0:49]*20;
errorbar(I_in,FRs_L23IT/0.6,0,se_L23IT/0.6,'Color',colorTable{1},"LineWidth",0.5,"LineStyle","-");   % calculate the firing rates. (number of spikes in 600 ms).
plot(I_in,FRs_L23IT/0.6,'Color',colorTable{1},"LineWidth",0.5,"LineStyle","-","Marker",".","MarkerSize",markersize);
errorbar(I_in,FRs_ITrorb/0.6,0,se_ITrorb/0.6,'Color',colorTable{2});
plot(I_in,FRs_ITrorb/0.6,'Color',colorTable{2},"LineWidth",0.5,"LineStyle","-","Marker",".","MarkerSize",markersize);
errorbar(I_in,FRs_L5IT/0.6,0,se_L5IT/0.6,'Color',colorTable{3});
plot(I_in,FRs_L5IT/0.6,'Color',colorTable{3},"LineWidth",0.5,"LineStyle","-","Marker",".","MarkerSize",markersize);
errorbar(I_in,FRs_ven/0.6,0,se_ven/0.6,'Color',colorTable{4});
plot(I_in,FRs_ven/0.6,'Color',colorTable{4},"LineWidth",0.5,"LineStyle","-","Marker",".","MarkerSize",markersize);
errorbar(I_in,FRs_L5ET/0.6,0,se_L5ET/0.6,'Color',colorTable{5});
plot(I_in,FRs_L5ET/0.6,'Color',colorTable{5},"LineWidth",0.5,"LineStyle","-","Marker",".","MarkerSize",markersize);
errorbar(I_in,FRs_L6PC/0.6,0,se_L6PC/0.6,'Color',colorTable{6});
plot(I_in,FRs_L6PC/0.6,'Color',colorTable{6},"LineWidth",0.5,"LineStyle","-","Marker",".","MarkerSize",markersize);
set(gca,'xlim',[-40,400],'ylim',[0,15]);
set(gca,"XTick",0:100:600,"YTick",0:5:15);
xlabel("Input current intensity (pA)");
ylabel("Firing rate (spikes/s)");
legend(cellTypes);


%% ==================two-way ANOVA analysis followed by Bonferroni Correction ============
data = [FRs_L23IT(1:21)/0.6, FRs_ITrorb(1:21)/0.6, FRs_L5IT(1:21)/0.6, FRs_ven(1:21)/0.6, FRs_L5ET(1:21)/0.6,FRs_L6PC(1:21)/0.6];
[p,tbl,stats] = anova2(data);
if p(1)<0.05
    [c, m, h] = multcompare(stats, "Estimate","column", 'CriticalValueType', 'bonferroni');
    % 比较结果矩阵
    tbl2 = array2table(c,"VariableNames", ...
          ["Cell type A","Cell type B","Lower Limit","A-B","Upper Limit","Adjusted p-value"]);
    disp(tbl2);  
end





