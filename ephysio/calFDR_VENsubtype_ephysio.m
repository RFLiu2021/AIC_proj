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

cluster_L23IT = {1,15,5,2,0};
cluster_L5IT = {3,6,7,9,11,13,14,18};
cluster_VEN_L5ET =  {16};
cluster_VEN_L56CT = {8};
cluster_L6PC = {4,12,17};


ephsData_VENL= [];  VENL_names={};
ephsData_VENS = [];  VENS_names={};
ephsData_L5ET = []; L5ET_names = {};
ephsData_L56CT = []; L56CT_names = {};
for i = 1:length(predictClass)
    switch predictClass(i)
        case cluster_VEN_L5ET
            if contains(expertCellType{i},'VEN') || contains(expertCellType{i},'ven')
                ephsData_VENL(end+1,:) = ephsData(i,:);
                VENL_names{end+1} = cellNames{i};
            elseif contains(expertCellType{i},'PC') || contains(expertCellType{i},'pc')
                ephsData_L5ET(end+1,:) = ephsData(i,:);
                L5ET_names{end+1} = cellNames{i};
            end
        case cluster_VEN_L56CT
            if contains(expertCellType{i},'VEN') || contains(expertCellType{i},'ven')
                ephsData_VENS(end+1,:) = ephsData(i,:);
                VENS_names{end+1} = cellNames{i};
            elseif contains(expertCellType{i},'PC') || contains(expertCellType{i},'pc')
                ephsData_L56CT(end+1,:) = ephsData(i,:);
                L56CT_names{end+1} = cellNames{i};
            end
    end
end 

colorDict = {[255,56,56]/255,[255,41,251]/255,[21,21,211]/255};


% ==================calculate the FDR for VEN_L vs VEN_S across all ephysio features ============
[ave_venl, sd_venl] = imean(ephsData_VENL,1);
[ave_vens, sd_vens] = imean(ephsData_VENS,1);

for i = 1:nFeature   % get p-values of Wilcoxon ranksum test.
    [raw_p(i),h(i)] = ranksum( ephsData_VENL(:,i),ephsData_VENS(:,i) );
end
%calculate fdr p-values.
[fdr_qvalues] = mafdr(raw_p, 'BHFDR', true);
% create FDR p-values table.
resTable = table();
resTable.parameterNames=parameterNames';
resTable.VENL_ave=ave_venl;
resTable.VENL_sd=sd_venl;
resTable.VENS_ave=ave_vens;
resTable.VENS_sd= sd_vens;
resTable.FDR_pval = fdr_qvalues';
% sort the table.
[~,sorted_idx] = sort(resTable.FDR_pval);
sorted_resTable = resTable(sorted_idx,:);

% output the adjusted p-values table to an excel file.
writetable(sorted_resTable, './FDR_VENL_vs_VENS_results.xlsx','Sheet','FDR_ephysio','WriteRowNames', true);










