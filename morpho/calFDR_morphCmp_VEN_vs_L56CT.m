clear all;
close all;

[numInfo1, txtInfo1, table_all1] = xlsread('../ephysio/data_ephy/Ins_patchseq_ephs.xlsx','Sheet1');
[numInfo2, txtInfo2, table_all2] = xlsread('./dataMorph/morphFeatMatrix.xlsx','morphFeatMatrix');

cellNames1 = txtInfo1(2:end,1);
predictClass = numInfo1(:,1);
cellNames2 = txtInfo2(2:end,1);
morphData = numInfo2(:,4:end);
morphParaNames = txtInfo2(1,5:end);
expertCellTypeID = numInfo2(:,1);

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



nFeature = size(morphData_L56CT,2);
% ----------------------------FDR calculation VENL vs VENS---------------------------
raw_p = nan(1,nFeature);
for i = 1:nFeature    %get raw p values.
    [raw_p(i)] = ranksum( morphData_VENL(:,i),morphData_VENS(:,i) );
end
%calculate fdr p-values.
[fdr_qvalues] = mafdr(raw_p, 'BHFDR', true);

% create FDR p-values table.
resTable = table();
resTable.morphParaNames=morphParaNames';
resTable.FDR_pval = fdr_qvalues';
% sort the table.
[~,sorted_idx] = sort(resTable.FDR_pval);
sorted_resTable = resTable(sorted_idx,:);

% output the adjusted p-values table to an excel file.
writetable(sorted_resTable, './FDR_morph_VEN_vs_L56CT_results.xlsx','Sheet','FDR_VENL_vs_VENS','WriteRowNames', true);

% ----------------------------FDR calculation  VENL vs L56CT---------------------------
raw_p = nan(1,nFeature);
for i = 1:nFeature    %get raw p values.
    [raw_p(i)] = ranksum( morphData_VENL(:,i),morphData_L56CT(:,i) );
end
%calculate fdr p-values.
[fdr_qvalues] = mafdr(raw_p, 'BHFDR', true);

% create FDR p-values table.
resTable = table();
resTable.morphParaNames=morphParaNames';
resTable.FDR_pval = fdr_qvalues';
% sort the table.
[~,sorted_idx] = sort(resTable.FDR_pval);
sorted_resTable = resTable(sorted_idx,:);

% output the adjusted p-values table to an excel file.
writetable(sorted_resTable, './FDR_morph_VEN_vs_L56CT_results.xlsx','Sheet','FDR_VENL_vs_L56CT','WriteRowNames', true);


% ----------------------------FDR calculation  VENS vs L56CT---------------------------
raw_p = nan(1,nFeature);
for i = 1:nFeature    %get raw p values.
    [raw_p(i)] = ranksum( morphData_VENS(:,i),morphData_L56CT(:,i) );
end
%calculate fdr p-values.
[fdr_qvalues] = mafdr(raw_p, 'BHFDR', true);

% create FDR p-values table.
resTable = table();
resTable.morphParaNames=morphParaNames';
resTable.FDR_pval = fdr_qvalues';
% sort the table.
[~,sorted_idx] = sort(resTable.FDR_pval);
sorted_resTable = resTable(sorted_idx,:);

% output the adjusted p-values table to an excel file.
writetable(sorted_resTable, './FDR_morph_VEN_vs_L56CT_results.xlsx','Sheet','FDR_VENS_vs_L56CT','WriteRowNames', true);