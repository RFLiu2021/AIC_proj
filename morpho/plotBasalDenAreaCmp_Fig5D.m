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
ven_Long_files{end+1} = 'dataMorph/unSeqVENs_morph/VENL/180913C1.mat';
ven_Long_files{end+1} = 'dataMorph/unSeqVENs_morph/VENL/180913C3.mat';
ven_Long_files{end+1} = 'dataMorph/unSeqVENs_morph/VENL/180913C4.mat';
ven_Long_files{end+1} = 'dataMorph/unSeqVENs_morph/VENL/180921C1.mat';
ven_Long_files{end+1} = 'dataMorph/unSeqVENs_morph/VENL/180921C3.mat';
ven_Long_files{end+1} = 'dataMorph/unSeqVENs_morph/VENL/181018C2.mat';
ven_Long_files{end+1} = 'dataMorph/unSeqVENs_morph/VENL/181018C5.mat';
ven_Long_files{end+1} = 'dataMorph/unSeqVENs_morph/VENL/181206C1.mat';
ven_Long_files{end+1} = 'dataMorph/unSeqVENs_morph/VENL/190828HS11C2.mat';
ven_Long_files{end+1} = 'dataMorph/unSeqVENs_morph/VENL/20200505VEN1.mat';
ven_Long_files{end+1} = 'dataMorph/unSeqVENs_morph/VENL/Golgi-2R-VEN.mat';
ven_Long_files{end+1} = 'dataMorph/unSeqVENs_morph/VENL/20200505VEN2.mat';
ven_Long_files{end+1} = 'dataMorph/unSeqVENs_morph/VENL/20200503C8.mat';

ven_Long_files{end+1} = 'dataMorph/unSeqVENs_morph/VENL/20221008_VEN_1.mat';
ven_Long_files{end+1} = 'dataMorph/unSeqVENs_morph/VENL/20221008_VEN_3.mat';
ven_Long_files{end+1} = 'dataMorph/unSeqVENs_morph/VENL/20221008_VEN_9.mat';
ven_Long_files{end+1} = 'dataMorph/unSeqVENs_morph/VENL/20221008_VEN_10.mat';
ven_Long_files{end+1} = 'dataMorph/unSeqVENs_morph/VENL/20221008_VEN_13.mat';
ven_Long_files{end+1} = 'dataMorph/unSeqVENs_morph/VENL/20221008_VEN_16.mat';
ven_Long_files{end+1} = 'dataMorph/unSeqVENs_morph/VENL/20191122AS2.mat';
ven_Long_files{end+1} = 'dataMorph/unSeqVENs_morph/VENL/20220525_VEN_1.mat';
ven_Long_files{end+1} = 'dataMorph/unSeqVENs_morph/VENL/20220525_VEN_2.mat';
ven_Long_files{end+1} = 'dataMorph/unSeqVENs_morph/VENL/20220525_VEN_4.mat';
ven_Long_files{end+1} = 'dataMorph/unSeqVENs_morph/VENL/20220527_VEN_3.mat';
ven_Long_files{end+1} = 'dataMorph/unSeqVENs_morph/VENL/20220527_VEN_4.mat';
ven_Long_files{end+1} = 'dataMorph/unSeqVENs_morph/VENL/20220527_VEN_6.mat';
ven_Long_files{end+1} = 'dataMorph/unSeqVENs_morph/VENL/20220605_VEN_3.mat';
ven_Long_files{end+1} = 'dataMorph/unSeqVENs_morph/VENL/20230202VEN_1.mat';
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
ven_Short_files{end+1} = 'dataMorph/unSeqVENs_morph/VENS/180730C4.mat';
ven_Short_files{end+1} = 'dataMorph/unSeqVENs_morph/VENS/180802C8.mat';
ven_Short_files{end+1} = 'dataMorph/unSeqVENs_morph/VENS/181204C1.mat';
ven_Short_files{end+1} = 'dataMorph/unSeqVENs_morph/VENS/181204C2.mat';
ven_Short_files{end+1} = 'dataMorph/unSeqVENs_morph/VENS/181204C8.mat';
ven_Short_files{end+1} = 'dataMorph/unSeqVENs_morph/VENS/190828C2.mat';
ven_Short_files{end+1} = 'dataMorph/unSeqVENs_morph/VENS/190828C4.mat';
ven_Short_files{end+1} = 'dataMorph/unSeqVENs_morph/VENS/20200429C3.mat';
ven_Short_files{end+1} = 'dataMorph/unSeqVENs_morph/VENS/20200429C5.mat';
ven_Short_files{end+1} = 'dataMorph/unSeqVENs_morph/VENS/20191122C3.mat';
ven_Short_files{end+1} = 'dataMorph/unSeqVENs_morph/VENS/20191122C1.mat';
ven_Short_files{end+1} = 'dataMorph/unSeqVENs_morph/VENS/20200807C9.mat';
ven_Short_files{end+1} = 'dataMorph/unSeqVENs_morph/VENS/20200530C2.mat';
ven_Short_files{end+1} = 'dataMorph/unSeqVENs_morph/VENS/190829C1.mat';

ven_Short_files{end+1} = 'dataMorph/unSeqVENs_morph/VENS/20221008_VEN_11.mat';
ven_Short_files{end+1} = 'dataMorph/unSeqVENs_morph/VENS/180802S1C10.mat';
ven_Short_files{end+1} = 'dataMorph/unSeqVENs_morph/VENS/180802S5C4.mat';
ven_Short_files{end+1} = 'dataMorph/unSeqVENs_morph/VENS/181204s12c20.mat';
ven_Short_files{end+1} = 'dataMorph/unSeqVENs_morph/VENS/181204s12c25.mat';
ven_Short_files{end+1} = 'dataMorph/unSeqVENs_morph/VENS/181206S2C2.mat';
ven_Short_files{end+1} = 'dataMorph/unSeqVENs_morph/VENS/181206S7C3.mat';
ven_Short_files{end+1} = 'dataMorph/unSeqVENs_morph/VENS/181206S8C1.mat';
ven_Short_files{end+1} = 'dataMorph/unSeqVENs_morph/VENS/181206S8C2.mat';
ven_Short_files{end+1} = 'dataMorph/unSeqVENs_morph/VENS/181210S18C1.mat';

 
PC_files = {};
PC_files{end+1} = 'dataMorph/asc&mat/s485.mat';
PC_files{end+1} = 'dataMorph/asc&mat/s486.mat';
PC_files{end+1} = 'dataMorph/asc&mat/B557.mat';
PC_files{end+1} = 'dataMorph/asc&mat/B735.mat';
PC_files{end+1} = 'dataMorph/asc&mat/n46.mat';
PC_files{end+1} = 'dataMorph/asc&mat/n15.mat';
PC_files{end+1} = 'dataMorph/asc&mat/n47.mat';
PC_files{end+1} = 'dataMorph/asc&mat/n84.mat';
PC_files{end+1} = 'dataMorph/asc&mat/n190.mat';
PC_files{end+1} = 'dataMorph/asc&mat/s175.mat';
PC_files{end+1} = 'dataMorph/asc&mat/s218.mat';
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
PC_files{end+1} = 'dataMorph/asc&mat/n105.mat';
PC_files{end+1} = 'dataMorph/asc&mat/n112.mat';
PC_files{end+1} = 'dataMorph/asc&mat/n122.mat';
PC_files{end+1} = 'dataMorph/asc&mat/n126.mat';
PC_files{end+1} = 'dataMorph/asc&mat/n147.mat';
PC_files{end+1} = 'dataMorph/asc&mat/n149.mat';
PC_files{end+1} = 'dataMorph/asc&mat/n151.mat';
PC_files{end+1} = 'dataMorph/asc&mat/n209.mat';
PC_files{end+1} = 'dataMorph/asc&mat/s170.mat';
PC_files{end+1} = 'dataMorph/asc&mat/s177.mat';
PC_files{end+1} = 'dataMorph/asc&mat/s214.mat';
PC_files{end+1} = 'dataMorph/asc&mat/s218.mat';
PC_files{end+1} = 'dataMorph/asc&mat/s219.mat';
PC_files{end+1} = 'dataMorph/asc&mat/s221.mat';
PC_files{end+1} = 'dataMorph/unSeqVENs_morph/PC/181204s12c21.mat';
PC_files{end+1} = 'dataMorph/unSeqVENs_morph/PC/181204s12c23.mat';
PC_files{end+1} = 'dataMorph/unSeqVENs_morph/PC/190314s3c6.mat';
PC_files{end+1} = 'dataMorph/unSeqVENs_morph/PC/190315s7c4.mat';
PC_files{end+1} = 'dataMorph/unSeqVENs_morph/PC/190315s10c3.mat';
PC_files{end+1} = 'dataMorph/unSeqVENs_morph/PC/190315s10c4.mat';
PC_files{end+1} = 'dataMorph/unSeqVENs_morph/PC/20200429C4.mat';
PC_files{end+1} = 'dataMorph/unSeqVENs_morph/PC/20200429C2.mat';
PC_files{end+1} = 'dataMorph/unSeqVENs_morph/PC/20200429C3.mat';
PC_files{end+1} = 'dataMorph/unSeqVENs_morph/PC/20200429C6.mat';
PC_files{end+1} = 'dataMorph/unSeqVENs_morph/PC/20220525C2.mat';
PC_files{end+1} = 'dataMorph/unSeqVENs_morph/PC/20220525C4.mat';
PC_files{end+1} = 'dataMorph/unSeqVENs_morph/PC/20230202C1.mat';


nVENl = length(ven_Long_files);
basalDenArea_venl = [];
basalDenHeight_venl=[];
for i = 1:nVENl
    load(ven_Long_files{i});
%     [ ASC ] = calDendriteDensityMap_retina( ASC );
%     [ ASC ] = calDendriteField( ASC );
    basalDenArea_venl(end+1) = ASC.basalDendriticFieldArea;
    height_t = 0;
    for j = 1:length(ASC.Dendrites)
        if height_t < abs(min( ASC.Dendrites(j).data(:,4) ))
            height_t = abs(min( ASC.Dendrites(j).data(:,4) ));
        end
    end
    ASC.basalDenHright = height_t;
    basalDenHeight_venl(end+1) = height_t;
    
    save(ven_Long_files{i},'ASC');
    
    close all;
    clear ASC height_t;
end

nVENs = length(ven_Short_files);
basalDenArea_vens = [];
basalDenHeight_vens = [];
for i = 1:nVENs
    load(ven_Short_files{i});
%     [ ASC ] = calDendriteDensityMap_retina( ASC );
%     [ ASC ] = calDendriteField( ASC );
    basalDenArea_vens(end+1) = ASC.basalDendriticFieldArea;
    height_t = 0;
    for j = 1:length(ASC.Dendrites)
        if height_t < abs(min( ASC.Dendrites(j).data(:,4) ))
            height_t = abs(min( ASC.Dendrites(j).data(:,4) ));
        end
    end
    ASC.basalDenHright = height_t;
    basalDenHeight_vens(end+1) = height_t;
    
    save(ven_Short_files{i},'ASC');
    
    close all;
    clear ASC height_t;
end

nPC = length(PC_files);
basalDenArea_pc = [];
basalDenHeight_pc = [];
for i = 1:nPC
    load(PC_files{i});
%     [ ASC ] = calDendriteDensityMap_retina( ASC );
%     [ ASC ] = calDendriteField( ASC );
    basalDenArea_pc(end+1) = ASC.basalDendriticFieldArea;
    height_t = 0;
    for j = 1:length(ASC.Dendrites)
        if height_t < abs(min( ASC.Dendrites(j).data(:,4) ))
            height_t = abs(min( ASC.Dendrites(j).data(:,4) ));
        end
    end
    ASC.basalDenHright = height_t;
    basalDenHeight_pc(end+1) = height_t;
    
    save(PC_files{i},'ASC');
    
    close all;
    clear ASC height_t;
end

%%%% ============ plot the figures ========================
colorDict = {[255,56,56]/255,[255,41,251]/255,[21,21,211]/255};
markerSize = 100;

figure('position',[600,500,600,350]);
subplot(1,2,1); hold on;
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

subplot(1,2,2); hold on;
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

% % ----------------------------statistics---------------------------
nVENl = length(basalDenArea_venl);
nVENs = length(basalDenArea_vens);
nPC = length(basalDenArea_pc);

nMax = max([nVENl,nVENs,nPC]);
Y = nan(nMax,3);
Y(1:nVENl,1) = reshape(basalDenArea_venl,[nVENl,1]);
Y(1:nVENs,2) = reshape(basalDenArea_vens,[nVENs,1]);
Y(1:nPC,3) = reshape(basalDenArea_pc,[nPC,1]);
[p,tbl,stats] = kruskalwallis(Y);
c=multcompare(stats)

nMax = max([nVENl,nVENs,nPC]);
Y = nan(nMax,3);
Y(1:nVENl,1) = reshape(log(basalDenHeight_venl),[nVENl,1]);
Y(1:nVENs,2) = reshape(log(basalDenHeight_vens),[nVENs,1]);
Y(1:nPC,3) = reshape(log(basalDenHeight_pc),[nPC,1]);
[p,tbl,stats] = kruskalwallis(Y);
c=multcompare(stats)


% %%%%=====================results of calculation====================
% basalDenArea_pc = [0.3704    0.5157    0.6619    1.0682    0.4420    0.6792    0.6845    0.9511    0.5642    ...
%      0.6050    0.4867    0.6757    0.1958    0.2657    0.7138    0.9661    0.3704    0.5157 ...
%      0.3121    0.1436    0.8009    0.5021    0.4364    0.9189 ]*1e+5;
% basalDenArea_venl = [0.8944    0.9745    0.9452    0.6851    0.8862    0.7227    0.2215    0.8324    1.0589 ...
%     2.5341    2.3657    1.3252    2.4302    1.5980    1.2660    1.9320    1.8142    1.5697   ...
%     0.5428    1.1752    2.2352    0.7441]*1e+5;
% basalDenArea_vens = [ 0.7656    0.3929    0.7239    0.5466    0.2693    0.7507    0.8490    0.6942    0.7152  ...
%      0.9491    1.1649    1.6567    0.7220    1.0606    0.5341    0.4611    1.1809    0.9928 ...
%      0.7943    0.6578    0.6865    1.1386    0.4319]*1e+5;


% basalDenHeight_pc = [  170.0600  188.6200  166.2900  317.4700  194.0600  222.8200  147.1100  406.8600  164.0800  ...
%      204.9800  208.6800  180.9500   86.0500  143.5300  265.2200  495.6800  170.0600  188.6200   ...
%       176.0200  109.9800  126.9100  126.3300  158.3300  303.7500];
% basalDenHeight_vens = [529.6100  242.6100  177.9900  225.7000  214.2484  256.3100  245.1900  293.7700  172.0100 ...
%     293.2200  396.8600  601.6500  316.0600  328.1200  195.5100  225.0600  307.9800  284.1400 ...
%     257.5100  219.6200  206.4900  289.3400  180.4000];
% basalDenHeight_venl = [491.0500  596.1500  449.2100  397.8600  496.2000  265.2800  448.7500  643.9800  495.9600   ...
%      889.7600  800.2300  308.2000  476.2100  698.6600  482.1400  488.0100  572.2600  501.9500   ...
%      322.2800  593.4000  884.3500  459.7400];



% %% ================================plot the example neurons =====================================
% % plot the example VEN-L
% load('dataMorph/unSeqVENs_morph/180913C1.mat');
% [ ASC ] = calDendriteDensityMap_retina( ASC );
% [ ASC ] = calDendriteField( ASC );
% set(gca,'xlim',[-500,500], 'ylim',[-900,1100]);
% axis off;
% 
% load('dataMorph/unSeqVENs_morph/20200104C4.mat');
% [ ASC ] = calDendriteDensityMap_retina( ASC );
% [ ASC ] = calDendriteField( ASC );
% set(gca,'xlim',[-500,500], 'ylim',[-600,1400]);
% axis off;
% 
% load( 'dataMorph/asc&mat/s486.mat' );
% [ ASC ] = calDendriteDensityMap_retina( ASC );
% [ ASC ] = calDendriteField( ASC );
% set(gca,'xlim',[-500,500], 'ylim',[-500,1500]);
% axis off;



