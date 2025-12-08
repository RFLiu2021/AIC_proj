clear all;
% close all;

ven_Long_files = {};
ven_Long_files{end+1} = 'dataMorph/asc&mat/s173.mat';
ven_Long_files{end+1} = 'dataMorph/asc&mat/s174.mat';
ven_Long_files{end+1} = 'dataMorph/asc&mat/s427.mat';
ven_Long_files{end+1} = 'dataMorph/asc&mat/s475.mat';
ven_Long_files{end+1} = 'dataMorph/asc&mat/s476.mat';
ven_Long_files{end+1} = 'dataMorph/asc&mat/s647.mat';
ven_Long_files{end+1} = 'dataMorph/asc&mat/s649.mat';
ven_Long_files{end+1} = 'dataMorph/asc&mat/s673.mat';
ven_Long_files{end+1} = 'dataMorph/asc&mat/B553.mat';
ven_Long_files{end+1} = 'dataMorph/asc&mat/B721.mat';
ven_Long_files{end+1} = 'dataMorph/unSeqCell_morph/VENL/unSM1114.mat';
ven_Long_files{end+1} = 'dataMorph/unSeqCell_morph/VENL/unSM1115.mat';
ven_Long_files{end+1} = 'dataMorph/unSeqCell_morph/VENL/unSM1116.mat';
ven_Long_files{end+1} = 'dataMorph/unSeqCell_morph/VENL/unSM1117.mat';
ven_Long_files{end+1} = 'dataMorph/unSeqCell_morph/VENL/unSM1118.mat';
ven_Long_files{end+1} = 'dataMorph/unSeqCell_morph/VENL/unSM1119.mat';
ven_Long_files{end+1} = 'dataMorph/unSeqCell_morph/VENL/unSM1120.mat';
ven_Long_files{end+1} = 'dataMorph/unSeqCell_morph/VENL/unSM1121.mat';
ven_Long_files{end+1} = 'dataMorph/unSeqCell_morph/VENL/unSM1122.mat';
ven_Long_files{end+1} = 'dataMorph/unSeqCell_morph/VENL/unSM1123.mat';
ven_Long_files{end+1} = 'dataMorph/unSeqCell_morph/VENL/unSM1124.mat';
ven_Long_files{end+1} = 'dataMorph/unSeqCell_morph/VENL/unSM1125.mat';
ven_Long_files{end+1} = 'dataMorph/unSeqCell_morph/VENL/unSM1126.mat';
ven_Long_files{end+1} = 'dataMorph/unSeqCell_morph/VENL/unSM1127.mat';
ven_Long_files{end+1} = 'dataMorph/unSeqCell_morph/VENL/unSM1128.mat';
ven_Long_files{end+1} = 'dataMorph/unSeqCell_morph/VENL/unSM1129.mat';
ven_Long_files{end+1} = 'dataMorph/unSeqCell_morph/VENL/unSM1130.mat';
ven_Long_files{end+1} = 'dataMorph/unSeqCell_morph/VENL/unSM1131.mat';
ven_Long_files{end+1} = 'dataMorph/unSeqCell_morph/VENL/unSM1132.mat';
ven_Long_files{end+1} = 'dataMorph/unSeqCell_morph/VENL/unSM1133.mat';
ven_Long_files{end+1} = 'dataMorph/unSeqCell_morph/VENL/unSM1134.mat';
ven_Long_files{end+1} = 'dataMorph/unSeqCell_morph/VENL/unSM1135.mat';
ven_Long_files{end+1} = 'dataMorph/unSeqCell_morph/VENL/unSM1136.mat';
ven_Long_files{end+1} = 'dataMorph/unSeqCell_morph/VENL/unSM1137.mat';
ven_Long_files{end+1} = 'dataMorph/unSeqCell_morph/VENL/unSM1138.mat';
ven_Long_files{end+1} = 'dataMorph/unSeqCell_morph/VENL/unSM1139.mat';
ven_Long_files{end+1} = 'dataMorph/unSeqCell_morph/VENL/unSM1140.mat';
ven_Long_files{end+1} = 'dataMorph/unSeqCell_morph/VENL/unSM1141.mat';
% 
ven_Short_files = {};
ven_Short_files{end+1} = 'dataMorph/asc&mat/s471.mat';
ven_Short_files{end+1} = 'dataMorph/asc&mat/s495.mat';
ven_Short_files{end+1} = 'dataMorph/asc&mat/s538.mat';
ven_Short_files{end+1} = 'dataMorph/asc&mat/n71.mat';
ven_Short_files{end+1} = 'dataMorph/asc&mat/n143.mat';
ven_Short_files{end+1} = 'dataMorph/asc&mat/n154.mat';
ven_Short_files{end+1} = 'dataMorph/asc&mat/n155.mat';
ven_Short_files{end+1} = 'dataMorph/asc&mat/n176.mat';
ven_Short_files{end+1} = 'dataMorph/asc&mat/n179.mat';
ven_Short_files{end+1} = 'dataMorph/asc&mat/s466.mat';
ven_Short_files{end+1} = 'dataMorph/unSeqCell_morph/VENS/unSM1142.mat';
ven_Short_files{end+1} = 'dataMorph/unSeqCell_morph/VENS/unSM1143.mat';
ven_Short_files{end+1} = 'dataMorph/unSeqCell_morph/VENS/unSM1144.mat';
ven_Short_files{end+1} = 'dataMorph/unSeqCell_morph/VENS/unSM1145.mat';
ven_Short_files{end+1} = 'dataMorph/unSeqCell_morph/VENS/unSM1146.mat';
ven_Short_files{end+1} = 'dataMorph/unSeqCell_morph/VENS/unSM1147.mat';
ven_Short_files{end+1} = 'dataMorph/unSeqCell_morph/VENS/unSM1148.mat';
ven_Short_files{end+1} = 'dataMorph/unSeqCell_morph/VENS/unSM1149.mat';
ven_Short_files{end+1} = 'dataMorph/unSeqCell_morph/VENS/unSM1150.mat';
ven_Short_files{end+1} = 'dataMorph/unSeqCell_morph/VENS/unSM1151.mat';
ven_Short_files{end+1} = 'dataMorph/unSeqCell_morph/VENS/unSM1152.mat';
ven_Short_files{end+1} = 'dataMorph/unSeqCell_morph/VENS/unSM1153.mat';
ven_Short_files{end+1} = 'dataMorph/unSeqCell_morph/VENS/unSM1154.mat';
ven_Short_files{end+1} = 'dataMorph/unSeqCell_morph/VENS/unSM1155.mat';

ven_Short_files{end+1} = 'dataMorph/unSeqCell_morph/VENS/unSM1156.mat';
ven_Short_files{end+1} = 'dataMorph/unSeqCell_morph/VENS/unSM1157.mat';
ven_Short_files{end+1} = 'dataMorph/unSeqCell_morph/VENS/unSM1158.mat';
ven_Short_files{end+1} = 'dataMorph/unSeqCell_morph/VENS/unSM1159.mat';
ven_Short_files{end+1} = 'dataMorph/unSeqCell_morph/VENS/unSM1160.mat';
ven_Short_files{end+1} = 'dataMorph/unSeqCell_morph/VENS/unSM1161.mat';
ven_Short_files{end+1} = 'dataMorph/unSeqCell_morph/VENS/unSM1162.mat';
ven_Short_files{end+1} = 'dataMorph/unSeqCell_morph/VENS/unSM1163.mat';
ven_Short_files{end+1} = 'dataMorph/unSeqCell_morph/VENS/unSM1164.mat';
ven_Short_files{end+1} = 'dataMorph/unSeqCell_morph/VENS/unSM1165.mat';



PC_files = {};
PC_files{end+1} = 'dataMorph/asc&mat/s485.mat';
PC_files{end+1} = 'dataMorph/asc&mat/s486.mat';
PC_files{end+1} = 'dataMorph/asc&mat/B557.mat';
PC_files{end+1} = 'dataMorph/asc&mat/B725.mat';
PC_files{end+1} = 'dataMorph/asc&mat/B735.mat';
PC_files{end+1} = 'dataMorph/asc&mat/B333.mat';
PC_files{end+1} = 'dataMorph/asc&mat/n15.mat';
PC_files{end+1} = 'dataMorph/asc&mat/n47.mat';
PC_files{end+1} = 'dataMorph/asc&mat/n84.mat';
PC_files{end+1} = 'dataMorph/asc&mat/n190.mat';
PC_files{end+1} = 'dataMorph/asc&mat/n147.mat';
PC_files{end+1} = 'dataMorph/asc&mat/s623.mat';
PC_files{end+1} = 'dataMorph/unSeqCell_morph/PC/unSM1101.mat';
PC_files{end+1} = 'dataMorph/unSeqCell_morph/PC/unSM1102.mat';
PC_files{end+1} = 'dataMorph/unSeqCell_morph/PC/unSM1103.mat';
PC_files{end+1} = 'dataMorph/unSeqCell_morph/PC/unSM1104.mat';
PC_files{end+1} = 'dataMorph/unSeqCell_morph/PC/unSM1105.mat';
PC_files{end+1} = 'dataMorph/unSeqCell_morph/PC/unSM1106.mat';
PC_files{end+1} = 'dataMorph/unSeqCell_morph/PC/unSM1107.mat';
PC_files{end+1} = 'dataMorph/unSeqCell_morph/PC/unSM1108.mat';
PC_files{end+1} = 'dataMorph/unSeqCell_morph/PC/unSM1109.mat';
PC_files{end+1} = 'dataMorph/unSeqCell_morph/PC/unSM1110.mat';
PC_files{end+1} = 'dataMorph/unSeqCell_morph/PC/unSM1111.mat';
PC_files{end+1} = 'dataMorph/unSeqCell_morph/PC/unSM1112.mat';
PC_files{end+1} = 'dataMorph/unSeqCell_morph/PC/unSM1113.mat';

% % =================================calculate sholl of L5ITs ===========================
PC_files{end+1} = 'dataMorph/asc&mat/s175.mat';
PC_files{end+1} = 'dataMorph/asc&mat/s222.mat';
PC_files{end+1} = 'dataMorph/asc&mat/s225.mat';
PC_files{end+1} = 'dataMorph/asc&mat/s438.mat';
PC_files{end+1} = 'dataMorph/asc&mat/s444.mat';
PC_files{end+1} = 'dataMorph/asc&mat/s481.mat';
PC_files{end+1} = 'dataMorph/asc&mat/s485.mat';
PC_files{end+1} = 'dataMorph/asc&mat/s486.mat';
PC_files{end+1} = 'dataMorph/asc&mat/s516.mat';
PC_files{end+1} = 'dataMorph/asc&mat/s533.mat';
PC_files{end+1} = 'dataMorph/asc&mat/s535.mat';
PC_files{end+1} = 'dataMorph/asc&mat/s536.mat';
PC_files{end+1} = 'dataMorph/asc&mat/s543.mat';
PC_files{end+1} = 'dataMorph/asc&mat/s549.mat';
PC_files{end+1} = 'dataMorph/asc&mat/s422.mat';
PC_files{end+1} = 'dataMorph/asc&mat/s467.mat';
PC_files{end+1} = 'dataMorph/asc&mat/s506.mat';
PC_files{end+1} = 'dataMorph/asc&mat/s534.mat';
PC_files{end+1} = 'dataMorph/asc&mat/s551.mat';
PC_files{end+1} = 'dataMorph/asc&mat/s620.mat';
PC_files{end+1} = 'dataMorph/asc&mat/s628.mat';
PC_files{end+1} = 'dataMorph/asc&mat/s629.mat';
PC_files{end+1} = 'dataMorph/asc&mat/s631.mat';
PC_files{end+1} = 'dataMorph/asc&mat/s655.mat';
PC_files{end+1} = 'dataMorph/asc&mat/s657.mat';
PC_files{end+1} = 'dataMorph/asc&mat/s680.mat';
PC_files{end+1} = 'dataMorph/asc&mat/s690.mat';
PC_files{end+1} = 'dataMorph/asc&mat/s696.mat';
PC_files{end+1} = 'dataMorph/asc&mat/s720.mat';
PC_files{end+1} = 'dataMorph/asc&mat/s722.mat';
PC_files{end+1} = 'dataMorph/asc&mat/s725.mat';
PC_files{end+1} = 'dataMorph/asc&mat/s727.mat';
PC_files{end+1} = 'dataMorph/asc&mat/s732.mat';
PC_files{end+1} = 'dataMorph/asc&mat/s734.mat';
PC_files{end+1} = 'dataMorph/asc&mat/B555.mat';
PC_files{end+1} = 'dataMorph/asc&mat/B566.mat';
PC_files{end+1} = 'dataMorph/asc&mat/B569.mat';
PC_files{end+1} = 'dataMorph/asc&mat/B570.mat';
PC_files{end+1} = 'dataMorph/asc&mat/B723.mat';
PC_files{end+1} = 'dataMorph/asc&mat/n12.mat';
PC_files{end+1} = 'dataMorph/asc&mat/n19.mat';
PC_files{end+1} = 'dataMorph/asc&mat/n25.mat';
PC_files{end+1} = 'dataMorph/asc&mat/n34.mat';
PC_files{end+1} = 'dataMorph/asc&mat/n52.mat';
PC_files{end+1} = 'dataMorph/asc&mat/n53.mat';
PC_files{end+1} = 'dataMorph/asc&mat/n64.mat';
PC_files{end+1} = 'dataMorph/asc&mat/n65.mat';
PC_files{end+1} = 'dataMorph/asc&mat/n87.mat';
PC_files{end+1} = 'dataMorph/asc&mat/n105.mat';
PC_files{end+1} = 'dataMorph/asc&mat/n112.mat';
PC_files{end+1} = 'dataMorph/asc&mat/n122.mat';
PC_files{end+1} = 'dataMorph/asc&mat/n126.mat';
PC_files{end+1} = 'dataMorph/asc&mat/n149.mat';
PC_files{end+1} = 'dataMorph/asc&mat/n209.mat';
PC_files{end+1} = 'dataMorph/asc&mat/s170.mat';
PC_files{end+1} = 'dataMorph/asc&mat/s177.mat';
PC_files{end+1} = 'dataMorph/asc&mat/s214.mat';
PC_files{end+1} = 'dataMorph/asc&mat/s218.mat';
PC_files{end+1} = 'dataMorph/asc&mat/s219.mat';
PC_files{end+1} = 'dataMorph/asc&mat/s221.mat';


nVENl = length(ven_Long_files);
basalDenArea_venl = [];
basalDenHeight_venl=[];
axon_dist_venl = [];
for i = 1:nVENl
    load(ven_Long_files{i});
    % [ ASC ] = calDendriteDensityMap_retina( ASC );
    % [ ASC ] = calDendriteField( ASC );
    basalDenArea_venl(end+1) = ASC.basalDendriticFieldArea;
    height_t = 0;
    for j = 1:length(ASC.Dendrites)
        if height_t < abs(min( ASC.Dendrites(j).data(:,4) ))
            height_t = abs(min( ASC.Dendrites(j).data(:,4) ));
        end
    end
    ASC.basalDenHeight = height_t;
    basalDenHeight_venl(end+1) = height_t;
    if ~isempty(ASC.axonDist)
        axon_dist_venl(end+1) = ASC.axonDist;
    else
        axon_dist_venl(end+1) = nan;
    end

    % save(ven_Long_files{i},'ASC');
    
    close all;
    clear ASC height_t;
end

nVENs = length(ven_Short_files);
basalDenArea_vens = [];
basalDenHeight_vens = [];
axon_dist_vens = [];
for i = 1:nVENs
    load(ven_Short_files{i});
    % [ ASC ] = calDendriteDensityMap_retina( ASC );
    % [ ASC ] = calDendriteField( ASC );
    basalDenArea_vens(end+1) = ASC.basalDendriticFieldArea;
    height_t = 0;
    for j = 1:length(ASC.Dendrites)
        if height_t < abs(min( ASC.Dendrites(j).data(:,4) ))
            height_t = abs(min( ASC.Dendrites(j).data(:,4) ));
        end
    end
    ASC.basalDenHeight = height_t;
    basalDenHeight_vens(end+1) = height_t;
    if ~isempty(ASC.axonDist)
        axon_dist_vens(end+1) = ASC.axonDist;
    else
        axon_dist_vens(end+1) = nan;
    end

    % save(ven_Short_files{i},'ASC');
    
    close all;
    clear ASC height_t;
end

nPC = length(PC_files);
basalDenArea_pc = [];
basalDenHeight_pc = [];
axon_dist_pc = [];
for i = 1:nPC
    load(PC_files{i});
    % [ ASC ] = calDendriteDensityMap_retina( ASC );
    % [ ASC ] = calDendriteField( ASC );
    basalDenArea_pc(end+1) = ASC.basalDendriticFieldArea;
    height_t = 0;
    for j = 1:length(ASC.Dendrites)
        if height_t < abs(min( ASC.Dendrites(j).data(:,4) ))
            height_t = abs(min( ASC.Dendrites(j).data(:,4) ));
        end
    end
    ASC.basalDenHeight = height_t;
    basalDenHeight_pc(end+1) = height_t;
    if ~isempty(ASC.axonDist)
        axon_dist_pc(end+1) = ASC.axonDist;
    else
        axon_dist_pc(end+1) = nan;
    end

    % save(PC_files{i},'ASC');
    
    close all;
    clear ASC height_t;
end






%%%% ============ plot the figures ========================
colorDict = {[255,56,56]/255,[255,166,0]/255,[21,21,211]/255};
markerSize = 100;

figure('position',[600,500,900,350]);
subplot(1,3,1); hold on;
nMax = max([nVENl,nVENs,nPC]);
Y = nan(nMax,3);
Y(1:nVENl,1) = reshape(basalDenHeight_venl,[nVENl,1]);
Y(1:nVENs,2) = reshape(basalDenHeight_vens,[nVENs,1]);
Y(1:nPC,3) = reshape(basalDenHeight_pc,[nPC,1]);
scatter(ones(nVENl,1),basalDenHeight_venl,markerSize,'marker','.','jitter','on','jitterAmount',0.1,'markerEdgeColor',colorDict{1});
scatter(2*ones(nVENs,1),basalDenHeight_vens,markerSize,'marker','.','jitter','on','jitterAmount',0.1,'markerEdgeColor',colorDict{2});
scatter(3*ones(nPC,1),basalDenHeight_pc,markerSize,'marker','.','jitter','on','jitterAmount',0.1,'markerEdgeColor',colorDict{3});
h=boxplot(Y,'Notch','off','Labels',{'VENL','VENS','L5 PC'},'whisker',1.5,'symbol','','color','k');
set(h,{'linew'},{1.0});
set(gca,'xlim',[0,4],'ylim',[0,1200],'xTick', [1,2,3], 'yTick',0:300:1200,'box','off');
ylabel('Basal dendritic height (um)');
title('Basal dendritic height');
hold off;

subplot(1,3,2); hold on;
nMax = max([nVENl,nVENs,nPC]);
Y = nan(nMax,3);
Y(1:nVENl,1) = reshape(basalDenArea_venl,[nVENl,1]);
Y(1:nVENs,2) = reshape(basalDenArea_vens,[nVENs,1]);
Y(1:nPC,3) = reshape(basalDenArea_pc,[nPC,1]);
scatter(ones(nVENl,1),basalDenArea_venl,markerSize,'marker','.','jitter','on','jitterAmount',0.1,'markerEdgeColor',colorDict{1});
scatter(2*ones(nVENs,1),basalDenArea_vens,markerSize,'marker','.','jitter','on','jitterAmount',0.1,'markerEdgeColor',colorDict{2});
scatter(3*ones(nPC,1),basalDenArea_pc,markerSize,'marker','.','jitter','on','jitterAmount',0.1,'markerEdgeColor',colorDict{3});
h=boxplot(Y,'Notch','off','Labels',{'VENL','VENS','L5 PC'},'whisker',1.5,'symbol','','color','k');
set(h,{'linew'},{1.0});
set(gca,'xlim',[0,4],'ylim',[0,400000],'xTick', [1,2,3], 'yTick',[0:10:40]*10e+3,'box','off');
ylabel('Basal dendritic area (um2)');
title('Basal dendritic area');
hold off;

subplot(1,3,3); hold on;
nMax = max([nVENl,nVENs,nPC]);
Y = nan(nMax,3);
Y(1:nVENl,1) = reshape(axon_dist_venl,[nVENl,1]);
Y(1:nVENs,2) = reshape(axon_dist_vens,[nVENs,1]);
Y(1:nPC,3) = reshape(axon_dist_pc,[nPC,1]);
scatter(ones(nVENl,1),axon_dist_venl,markerSize,'marker','.','jitter','on','jitterAmount',0.1,'markerEdgeColor',colorDict{1});
scatter(2*ones(nVENs,1),axon_dist_vens,markerSize,'marker','.','jitter','on','jitterAmount',0.1,'markerEdgeColor',colorDict{2});
scatter(3*ones(nPC,1),axon_dist_pc,markerSize,'marker','.','jitter','on','jitterAmount',0.1,'markerEdgeColor',colorDict{3});
h=boxplot(Y,'Notch','off','Labels',{'VENL','VENS','L5 PC'},'whisker',1.5,'symbol','','color','k');
set(h,{'linew'},{1.0});
set(gca,'xlim',[0,4],'ylim',[-5,200],'xTick', [1,2,3], 'yTick',0:50:200,'box','off');
ylabel('Axon distance (um)');
title('Axon distance to soma');
hold off;




% % ----------------------------FDR calculation and statistics---------------------------
% [ave_venl, sd_venl] = imean(morphData_VENL,1);
% [ave_vens, sd_vens] = imean(morphData_VENS,1);
% 
% for i = 1:nFeature    %get raw p values.
%     [raw_p(i),h(i)] = ranksum( morphData_VENL(:,i),morphData_VENS(:,i) );
% end
% %calculate fdr p-values.
% [fdr_qvalues] = mafdr(raw_p, 'BHFDR', true);
% 
% % create FDR p-values table.
% resTable = table();
% resTable.morphParaNames=morphParaNames';
% resTable.VENL_ave=ave_venl;
% resTable.VENL_sd=sd_venl;
% resTable.VENS_ave=ave_vens;
% resTable.VENS_sd= sd_vens;
% resTable.FDR_pval = fdr_qvalues';
% % sort the table.
% [~,sorted_idx] = sort(resTable.FDR_pval);
% sorted_resTable = resTable(sorted_idx,:);
% 
% % output the adjusted p-values table to an excel file.
% writetable(sorted_resTable, './FDR_VENL_vs_VENS_results.xlsx','Sheet','FDR_morpho','WriteRowNames', true);


