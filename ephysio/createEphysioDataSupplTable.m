clear all;
close all;

[numInfo, txtInfo, table_all] = xlsread('./data_ephy/Ins_patchseq_ephs.xlsx','Sheet1');

predictClass = numInfo(:,1);
cellNames = txtInfo(2:end,1);
parameterNames = txtInfo(1,4:end);
expertCellType = txtInfo(2:end,2);
nFeature = length(parameterNames);

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




%% =============sort data and output to excel file==============
% --------------------create L23IT data table----------------
[mean,std,~] = imean(ephsData_L23IT,1);
% create ephysio data table.
dtTable = table();
dtTable.parameterNames=parameterNames';
dtTable.mean = mean;
dtTable.std = std;

% output the adjusted p-values table to an excel file.
writetable(dtTable, './ephysioFeatData_byCellType.xlsx','Sheet','ephysFeat_L23IT','WriteRowNames', true);


% --------------------create ITrorb data table----------------
[mean,std,~] = imean(ephsData_ITrorb,1);
% create ephysio data table.
dtTable = table();
dtTable.parameterNames=parameterNames';
dtTable.mean = mean;
dtTable.std = std;

% output the adjusted p-values table to an excel file.
writetable(dtTable, './ephysioFeatData_byCellType.xlsx','Sheet','ephysFeat_ITrorb','WriteRowNames', true);


% --------------------create L5IT data table----------------
[mean,std,~] = imean(ephsData_L5IT,1);
% create ephysio data table.
dtTable = table();
dtTable.parameterNames=parameterNames';
dtTable.mean = mean;
dtTable.std = std;

% output the adjusted p-values table to an excel file.
writetable(dtTable, './ephysioFeatData_byCellType.xlsx','Sheet','ephysFeat_L5IT','WriteRowNames', true);


% --------------------create VEN data table----------------
[mean,std,~] = imean(ephsData_VEN,1);
% create ephysio data table.
dtTable = table();
dtTable.parameterNames=parameterNames';
dtTable.mean = mean;
dtTable.std = std;

% output the adjusted p-values table to an excel file.
writetable(dtTable, './ephysioFeatData_byCellType.xlsx','Sheet','ephysFeat_VEN','WriteRowNames', true);


% --------------------create L5ET data table----------------
[mean,std,~] = imean(ephsData_L5ET,1);
% create ephysio data table.
dtTable = table();
dtTable.parameterNames=parameterNames';
dtTable.mean = mean;
dtTable.std = std;

% output the adjusted p-values table to an excel file.
writetable(dtTable, './ephysioFeatData_byCellType.xlsx','Sheet','ephysFeat_L5ET','WriteRowNames', true);


% --------------------create L6PC data table----------------
[mean,std,~] = imean(ephsData_L6PC,1);
% create ephysio data table.
dtTable = table();
dtTable.parameterNames=parameterNames';
dtTable.mean = mean;
dtTable.std = std;

% output the adjusted p-values table to an excel file.
writetable(dtTable, './ephysioFeatData_byCellType.xlsx','Sheet','ephysFeat_L6PC','WriteRowNames', true);


