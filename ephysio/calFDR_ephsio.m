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



% ===============create FDR pvalues list for L23IT vs others================
% prepare data.
ephsData_other=[ephsData_ITrorb;ephsData_L5IT; ephsData_VEN; ephsData_L5ET;ephsData_L6PC];
for i = 1:nFeature
    [raw_p(i),h(i)] = ranksum( ephsData_L23IT(:,i),ephsData_other(:,i) );
end
%calculate fdr p-values.
[fdr_qvalues] = mafdr(raw_p, 'BHFDR', true);
% create FDR p-values table.
resTable = table();
resTable.parameterNames=parameterNames';
resTable.FDR_pval = fdr_qvalues';
% sort the table.
[~,sorted_idx] = sort(resTable.FDR_pval);
sorted_resTable = resTable(sorted_idx,:);

% output the adjusted p-values table to an excel file.
writetable(sorted_resTable, './ephysioFeat_FDR_results.xlsx','Sheet','FDR_L23IT_vs_other','WriteRowNames', true);



% ===============create FDR pvalues list for ITrorb vs others================
% prepare data.
ephsData_other=[ephsData_L23IT;ephsData_L5IT; ephsData_VEN; ephsData_L5ET;ephsData_L6PC];
for i = 1:nFeature
    [raw_p(i),h(i)] = ranksum( ephsData_ITrorb(:,i),ephsData_other(:,i) );
end
%calculate fdr p-values.
[fdr_qvalues] = mafdr(raw_p, 'BHFDR', true);
% create FDR p-values table.
resTable = table();
resTable.parameterNames=parameterNames';
resTable.FDR_pval = fdr_qvalues';
% sort the table.
[~,sorted_idx] = sort(resTable.FDR_pval);
sorted_resTable = resTable(sorted_idx,:);

% output the adjusted p-values table to an excel file.
writetable(sorted_resTable, './ephysioFeat_FDR_results.xlsx','Sheet','FDR_ITrorb_vs_other','WriteRowNames', true);



% ===============create FDR pvalues list for L5IT vs others================
% prepare data.
ephsData_other=[ephsData_L23IT;ephsData_ITrorb; ephsData_VEN; ephsData_L5ET;ephsData_L6PC];
for i = 1:nFeature
    [raw_p(i),h(i)] = ranksum( ephsData_L5IT(:,i),ephsData_other(:,i) );
end
%calculate fdr p-values.
[fdr_qvalues] = mafdr(raw_p, 'BHFDR', true);
% create FDR p-values table.
resTable = table();
resTable.parameterNames=parameterNames';
resTable.FDR_pval = fdr_qvalues';
% sort the table.
[~,sorted_idx] = sort(resTable.FDR_pval);
sorted_resTable = resTable(sorted_idx,:);

% output the adjusted p-values table to an excel file.
writetable(sorted_resTable, './ephysioFeat_FDR_results.xlsx','Sheet','FDR_L5IT_vs_other','WriteRowNames', true);



% ===============create FDR pvalues list for VEN vs others================
% prepare data.
ephsData_other=[ephsData_L23IT;ephsData_ITrorb; ephsData_L5IT; ephsData_L5ET;ephsData_L6PC];
for i = 1:nFeature
    [raw_p(i),h(i)] = ranksum( ephsData_VEN(:,i),ephsData_other(:,i) );
end
%calculate fdr p-values.
[fdr_qvalues] = mafdr(raw_p, 'BHFDR', true);
% create FDR p-values table.
resTable = table();
resTable.parameterNames=parameterNames';
resTable.FDR_pval = fdr_qvalues';
% sort the table.
[~,sorted_idx] = sort(resTable.FDR_pval);
sorted_resTable = resTable(sorted_idx,:);

% output the adjusted p-values table to an excel file.
writetable(sorted_resTable, './ephysioFeat_FDR_results.xlsx','Sheet','FDR_VEN_vs_other','WriteRowNames', true);


% ===============create FDR pvalues list for L5ET vs others================
% prepare data.
ephsData_other=[ephsData_L23IT;ephsData_ITrorb; ephsData_L5IT; ephsData_VEN;ephsData_L6PC];
for i = 1:nFeature
    [raw_p(i),h(i)] = ranksum( ephsData_L5ET(:,i),ephsData_other(:,i) );
end
%calculate fdr p-values.
[fdr_qvalues] = mafdr(raw_p, 'BHFDR', true);
% create FDR p-values table.
resTable = table();
resTable.parameterNames=parameterNames';
resTable.FDR_pval = fdr_qvalues';
% sort the table.
[~,sorted_idx] = sort(resTable.FDR_pval);
sorted_resTable = resTable(sorted_idx,:);

% output the adjusted p-values table to an excel file.
writetable(sorted_resTable, './ephysioFeat_FDR_results.xlsx','Sheet','FDR_L5ET_vs_other','WriteRowNames', true);



% ===============create FDR pvalues list for L6PC vs others================
% prepare data.
ephsData_other=[ephsData_L23IT;ephsData_ITrorb; ephsData_L5IT; ephsData_VEN; ephsData_L5ET];
for i = 1:nFeature
    [raw_p(i),h(i)] = ranksum( ephsData_L6PC(:,i),ephsData_other(:,i) );
end
%calculate fdr p-values.
[fdr_qvalues] = mafdr(raw_p, 'BHFDR', true);
% create FDR p-values table.
resTable = table();
resTable.parameterNames=parameterNames';
resTable.FDR_pval = fdr_qvalues';
% sort the table.
[~,sorted_idx] = sort(resTable.FDR_pval);
sorted_resTable = resTable(sorted_idx,:);

% output the adjusted p-values table to an excel file.
writetable(sorted_resTable, './ephysioFeat_FDR_results.xlsx','Sheet','FDR_L6PC_vs_other','WriteRowNames', true);



%% ==================================================================================
% ==================calculate the FDR for pair-wise comparisons (L23IT)============

% --------------------create FDR pvalues list for L23IT vs ITrorb----------------
for i = 1:nFeature   % get p-values of Wilcoxon ranksum test.
    [raw_p(i),h(i)] = ranksum( ephsData_L23IT(:,i),ephsData_ITrorb(:,i) );
end
%calculate fdr p-values.
[fdr_qvalues] = mafdr(raw_p, 'BHFDR', true);
% create FDR p-values table.
resTable = table();
resTable.parameterNames=parameterNames';
resTable.FDR_pval = fdr_qvalues';
% sort the table.
[~,sorted_idx] = sort(resTable.FDR_pval);
sorted_resTable = resTable(sorted_idx,:);

% output the adjusted p-values table to an excel file.
writetable(sorted_resTable, './ephysioFeat_FDR_results.xlsx','Sheet','FDR_L23IT_vs_ITrorb','WriteRowNames', true);

% --------------------create FDR pvalues list for L23IT vs L5IT----------------
for i = 1:nFeature   % get p-values of Wilcoxon ranksum test.
    [raw_p(i),h(i)] = ranksum( ephsData_L23IT(:,i),ephsData_L5IT(:,i) );
end
%calculate fdr p-values.
[fdr_qvalues] = mafdr(raw_p, 'BHFDR', true);
% create FDR p-values table.
resTable = table();
resTable.parameterNames=parameterNames';
resTable.FDR_pval = fdr_qvalues';
% sort the table.
[~,sorted_idx] = sort(resTable.FDR_pval);
sorted_resTable = resTable(sorted_idx,:);

% output the adjusted p-values table to an excel file.
writetable(sorted_resTable, './ephysioFeat_FDR_results.xlsx','Sheet','FDR_L23IT_vs_L5IT','WriteRowNames', true);

% --------------------create FDR pvalues list for L23IT vs VEN----------------
for i = 1:nFeature   % get p-values of Wilcoxon ranksum test.
    [raw_p(i),h(i)] = ranksum( ephsData_L23IT(:,i),ephsData_VEN(:,i) );
end
%calculate fdr p-values.
[fdr_qvalues] = mafdr(raw_p, 'BHFDR', true);
% create FDR p-values table.
resTable = table();
resTable.parameterNames=parameterNames';
resTable.FDR_pval = fdr_qvalues';
% sort the table.
[~,sorted_idx] = sort(resTable.FDR_pval);
sorted_resTable = resTable(sorted_idx,:);

% output the adjusted p-values table to an excel file.
writetable(sorted_resTable, './ephysioFeat_FDR_results.xlsx','Sheet','FDR_L23IT_vs_VEN','WriteRowNames', true);

% --------------------create FDR pvalues list for L23IT vs L5ET----------------
for i = 1:nFeature   % get p-values of Wilcoxon ranksum test.
    [raw_p(i),h(i)] = ranksum( ephsData_L23IT(:,i),ephsData_L5ET(:,i) );
end
%calculate fdr p-values.
[fdr_qvalues] = mafdr(raw_p, 'BHFDR', true);
% create FDR p-values table.
resTable = table();
resTable.parameterNames=parameterNames';
resTable.FDR_pval = fdr_qvalues';
% sort the table.
[~,sorted_idx] = sort(resTable.FDR_pval);
sorted_resTable = resTable(sorted_idx,:);

% output the adjusted p-values table to an excel file.
writetable(sorted_resTable, './ephysioFeat_FDR_results.xlsx','Sheet','FDR_L23IT_vs_L5ET','WriteRowNames', true);

% --------------------create FDR pvalues list for L23IT vs L6PC----------------
for i = 1:nFeature   % get p-values of Wilcoxon ranksum test.
    [raw_p(i),h(i)] = ranksum( ephsData_L23IT(:,i),ephsData_L6PC(:,i) );
end
%calculate fdr p-values.
[fdr_qvalues] = mafdr(raw_p, 'BHFDR', true);
% create FDR p-values table.
resTable = table();
resTable.parameterNames=parameterNames';
resTable.FDR_pval = fdr_qvalues';
% sort the table.
[~,sorted_idx] = sort(resTable.FDR_pval);
sorted_resTable = resTable(sorted_idx,:);

% output the adjusted p-values table to an excel file.
writetable(sorted_resTable, './ephysioFeat_FDR_results.xlsx','Sheet','FDR_L23IT_vs_L6PC','WriteRowNames', true);

% ==================calculate the FDR for pair-wise comparisons (ITrorb)============
% --------------------create FDR pvalues list for ITrorb vs L5IT ----------------
for i = 1:nFeature   % get p-values of Wilcoxon ranksum test.
    [raw_p(i),h(i)] = ranksum( ephsData_ITrorb(:,i),ephsData_L5IT(:,i) );
end
%calculate fdr p-values.
[fdr_qvalues] = mafdr(raw_p, 'BHFDR', true);
% create FDR p-values table.
resTable = table();
resTable.parameterNames=parameterNames';
resTable.FDR_pval = fdr_qvalues';
% sort the table.
[~,sorted_idx] = sort(resTable.FDR_pval);
sorted_resTable = resTable(sorted_idx,:);

% output the adjusted p-values table to an excel file.
writetable(sorted_resTable, './ephysioFeat_FDR_results.xlsx','Sheet','FDR_ITrorb_vs_L5IT','WriteRowNames', true);


% --------------------create FDR pvalues list for ITrorb vs VEN----------------
for i = 1:nFeature   % get p-values of Wilcoxon ranksum test.
    [raw_p(i),h(i)] = ranksum( ephsData_ITrorb(:,i),ephsData_VEN(:,i) );
end
%calculate fdr p-values.
[fdr_qvalues] = mafdr(raw_p, 'BHFDR', true);
% create FDR p-values table.
resTable = table();
resTable.parameterNames=parameterNames';
resTable.FDR_pval = fdr_qvalues';
% sort the table.
[~,sorted_idx] = sort(resTable.FDR_pval);
sorted_resTable = resTable(sorted_idx,:);

% output the adjusted p-values table to an excel file.
writetable(sorted_resTable, './ephysioFeat_FDR_results.xlsx','Sheet','FDR_ITrorb_vs_VEN','WriteRowNames', true);

% --------------------create FDR pvalues list for ITrorb vs L5ET----------------
for i = 1:nFeature   % get p-values of Wilcoxon ranksum test.
    [raw_p(i),h(i)] = ranksum( ephsData_ITrorb(:,i),ephsData_L5ET(:,i) );
end
%calculate fdr p-values.
[fdr_qvalues] = mafdr(raw_p, 'BHFDR', true);
% create FDR p-values table.
resTable = table();
resTable.parameterNames=parameterNames';
resTable.FDR_pval = fdr_qvalues';
% sort the table.
[~,sorted_idx] = sort(resTable.FDR_pval);
sorted_resTable = resTable(sorted_idx,:);

% output the adjusted p-values table to an excel file.
writetable(sorted_resTable, './ephysioFeat_FDR_results.xlsx','Sheet','FDR_ITrorb_vs_L5ET','WriteRowNames', true);

% --------------------create FDR pvalues list for ITrorb vs L6PC----------------
for i = 1:nFeature   % get p-values of Wilcoxon ranksum test.
    [raw_p(i),h(i)] = ranksum( ephsData_ITrorb(:,i),ephsData_L6PC(:,i) );
end
%calculate fdr p-values.
[fdr_qvalues] = mafdr(raw_p, 'BHFDR', true);
% create FDR p-values table.
resTable = table();
resTable.parameterNames=parameterNames';
resTable.FDR_pval = fdr_qvalues';
% sort the table.
[~,sorted_idx] = sort(resTable.FDR_pval);
sorted_resTable = resTable(sorted_idx,:);

% output the adjusted p-values table to an excel file.
writetable(sorted_resTable, './ephysioFeat_FDR_results.xlsx','Sheet','FDR_ITrorb_vs_L6PC','WriteRowNames', true);


% ==================calculate the FDR for pair-wise comparisons (L5IT)============
% --------------------create FDR pvalues list for L5IT vs VEN ----------------
for i = 1:nFeature   % get p-values of Wilcoxon ranksum test.
    [raw_p(i),h(i)] = ranksum( ephsData_L5IT(:,i),ephsData_VEN(:,i) );
end
%calculate fdr p-values.
[fdr_qvalues] = mafdr(raw_p, 'BHFDR', true);
% create FDR p-values table.
resTable = table();
resTable.parameterNames=parameterNames';
resTable.FDR_pval = fdr_qvalues';
% sort the table.
[~,sorted_idx] = sort(resTable.FDR_pval);
sorted_resTable = resTable(sorted_idx,:);

% output the adjusted p-values table to an excel file.
writetable(sorted_resTable, './ephysioFeat_FDR_results.xlsx','Sheet','FDR_L5IT_vs_VEN','WriteRowNames', true);

% --------------------create FDR pvalues list for L5IT vs L5ET ----------------
for i = 1:nFeature   % get p-values of Wilcoxon ranksum test.
    [raw_p(i),h(i)] = ranksum( ephsData_L5IT(:,i),ephsData_L5ET(:,i) );
end
%calculate fdr p-values.
[fdr_qvalues] = mafdr(raw_p, 'BHFDR', true);
% create FDR p-values table.
resTable = table();
resTable.parameterNames=parameterNames';
resTable.FDR_pval = fdr_qvalues';
% sort the table.
[~,sorted_idx] = sort(resTable.FDR_pval);
sorted_resTable = resTable(sorted_idx,:);

% output the adjusted p-values table to an excel file.
writetable(sorted_resTable, './ephysioFeat_FDR_results.xlsx','Sheet','FDR_L5IT_vs_L5ET','WriteRowNames', true);

% --------------------create FDR pvalues list for L5IT vs L6PC ----------------
for i = 1:nFeature   % get p-values of Wilcoxon ranksum test.
    [raw_p(i),h(i)] = ranksum( ephsData_L5IT(:,i),ephsData_L6PC(:,i) );
end
%calculate fdr p-values.
[fdr_qvalues] = mafdr(raw_p, 'BHFDR', true);
% create FDR p-values table.
resTable = table();
resTable.parameterNames=parameterNames';
resTable.FDR_pval = fdr_qvalues';
% sort the table.
[~,sorted_idx] = sort(resTable.FDR_pval);
sorted_resTable = resTable(sorted_idx,:);

% output the adjusted p-values table to an excel file.
writetable(sorted_resTable, './ephysioFeat_FDR_results.xlsx','Sheet','FDR_L5IT_vs_L6PC','WriteRowNames', true);


% ==================calculate the FDR for pair-wise comparisons (VEN)============
% --------------------create FDR pvalues list for VEN vs L5ET ----------------
for i = 1:nFeature   % get p-values of Wilcoxon ranksum test.
    [raw_p(i),h(i)] = ranksum( ephsData_VEN(:,i),ephsData_L5ET(:,i) );
end
%calculate fdr p-values.
[fdr_qvalues] = mafdr(raw_p, 'BHFDR', true);
% create FDR p-values table.
resTable = table();
resTable.parameterNames=parameterNames';
resTable.FDR_pval = fdr_qvalues';
% sort the table.
[~,sorted_idx] = sort(resTable.FDR_pval);
sorted_resTable = resTable(sorted_idx,:);

% output the adjusted p-values table to an excel file.
writetable(sorted_resTable, './ephysioFeat_FDR_results.xlsx','Sheet','FDR_VEN_vs_L5ET','WriteRowNames', true);

% --------------------create FDR pvalues list for VEN vs L6PC ----------------
for i = 1:nFeature   % get p-values of Wilcoxon ranksum test.
    [raw_p(i),h(i)] = ranksum( ephsData_VEN(:,i),ephsData_L6PC(:,i) );
end
%calculate fdr p-values.
[fdr_qvalues] = mafdr(raw_p, 'BHFDR', true);
% create FDR p-values table.
resTable = table();
resTable.parameterNames=parameterNames';
resTable.FDR_pval = fdr_qvalues';
% sort the table.
[~,sorted_idx] = sort(resTable.FDR_pval);
sorted_resTable = resTable(sorted_idx,:);

% output the adjusted p-values table to an excel file.
writetable(sorted_resTable, './ephysioFeat_FDR_results.xlsx','Sheet','FDR_VEN_vs_L6PC','WriteRowNames', true);


% ==================calculate the FDR for pair-wise comparisons (L5ET)============
% --------------------create FDR pvalues list for L5ET vs L6PC ----------------
for i = 1:nFeature   % get p-values of Wilcoxon ranksum test.
    [raw_p(i),h(i)] = ranksum( ephsData_L5ET(:,i),ephsData_L6PC(:,i) );
end
%calculate fdr p-values.
[fdr_qvalues] = mafdr(raw_p, 'BHFDR', true);
% create FDR p-values table.
resTable = table();
resTable.parameterNames=parameterNames';
resTable.FDR_pval = fdr_qvalues';
% sort the table.
[~,sorted_idx] = sort(resTable.FDR_pval);
sorted_resTable = resTable(sorted_idx,:);

% output the adjusted p-values table to an excel file.
writetable(sorted_resTable, './ephysioFeat_FDR_results.xlsx','Sheet','FDR_L5ET_vs_L6PC','WriteRowNames', true);

