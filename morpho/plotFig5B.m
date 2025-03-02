clear all; 
% close all;

% plot Fig 7B left pannel.
%---------- plot the example cells...-----------------
ven_long_example = 'dataMorph/asc&mat/B721.mat';
% ven_long_example = 'dataMorph/asc&mat/s427.mat';
ven_short_example = 'dataMorph/asc&mat/s495.mat';


load(ven_long_example);
lenVSthick_long = ASC.lenVSthick;
load(ven_short_example);
lenVSthick_short = ASC.lenVSthick;

figure();
plot(lenVSthick_short(:,1),lenVSthick_short(:,2)./max(lenVSthick_short(:,2)),'Color','#FFA500','lineWidth',1.0);
hold on;
plot(lenVSthick_long(:,1),lenVSthick_long(:,2)./max(lenVSthick_long(:,2)),'Color','#FF3838','lineWidth',1.0);
set(gca,'xLim',[0,800],'yLim',[0,1],'box','off','yTick',0:0.25:1);


% plot Fig 7B right
ven_Long_files = {};
ven_Long_files{end+1} = 'dataMorph/asc&mat/s174.mat';
ven_Long_files{end+1} = 'dataMorph/asc&mat/s427.mat';
ven_Long_files{end+1} = 'dataMorph/asc&mat/s475.mat';
ven_Long_files{end+1} = 'dataMorph/asc&mat/s476.mat';
ven_Long_files{end+1} = 'dataMorph/asc&mat/s647.mat';
ven_Long_files{end+1} = 'dataMorph/asc&mat/s649.mat';
ven_Long_files{end+1} = 'dataMorph/asc&mat/s673.mat';
ven_Long_files{end+1} = 'dataMorph/asc&mat/B553.mat';
ven_Long_files{end+1} = 'dataMorph/asc&mat/B721.mat';
ven_Short_files = {};
ven_Short_files{end+1} = 'dataMorph/asc&mat/s173.mat';
ven_Short_files{end+1} = 'dataMorph/asc&mat/s471.mat';
ven_Short_files{end+1} = 'dataMorph/asc&mat/s495.mat';
ven_Short_files{end+1} = 'dataMorph/asc&mat/s538.mat';
ven_Short_files{end+1} = 'dataMorph/asc&mat/n71.mat';
ven_Short_files{end+1} = 'dataMorph/asc&mat/n143.mat';
ven_Short_files{end+1} = 'dataMorph/asc&mat/n155.mat';
% ven_Short_files{end+1} = 'dataMorph/asc&mat/n201.mat';


% =======================VEN-L analysis=====================
% get the average values of lenVSthick and arrange the data in a matrix.
nVlong = length(ven_Long_files);
for i = 1:nVlong
    load(ven_Long_files{i});
    lenVSthick = ASC.lenVSthick;
    V_long_tem = nan(1000,1);
    
    maxR = max(lenVSthick(:,2));
    lenVSthick(:,2) = lenVSthick(:,2)./maxR;
    
    lenVSthick(:,1) =  round(lenVSthick(:,1));
    [C,ia,ic] = unique(lenVSthick(:,1));
    lenVSthick = lenVSthick(ia,:);
    V_long_tem = nan(1000,1);
    if lenVSthick(end,1) > 1000
        lenVSthick(lenVSthick(:,1)>1000,:) = [];
    end
    
    % fill the value gaps which filled with NaNs.
    V_tem(lenVSthick(:,1))=lenVSthick(:,2);
    zeroIdx = find(V_tem==0);
    valueIdx = find(V_tem~=0);
    for j = 1:length(zeroIdx)
        if zeroIdx(j) < valueIdx(1)
            V_tem(zeroIdx(j)) = V_tem(valueIdx(1));
        else
            idx_tem = find(valueIdx < zeroIdx(j));
            V_tem(zeroIdx(j)) = V_tem(valueIdx( idx_tem(end) ));
        end
    end
    
    V_long_tem(1:length(V_tem)) = V_tem';
    V_long(:,i) = V_long_tem;
    
    clear V_tem
end
[ave_V_long,sd_V_long,se_V_long] = imean(V_long,2);
idx_isnan = ~isnan(ave_V_long);
idx = find(idx_isnan == 1); idx(end) = 690; % delete the last several noise  points.
ave_V_long(idx(end)+1:end) = [];
sd_V_long(idx(end)+1:end) = [];
se_V_long(idx(end)+1:end) = [];
ave_V_long(isnan(ave_V_long)) = ave_V_long(idx(1));
sd_V_long(isnan(sd_V_long)) = sd_V_long(idx(1));
se_V_long(isnan(se_V_long)) = 0;

% =======================VEN-S analysis=====================
% get the average values of lenVSthick and arrange the data in a matrix.
nVShort = length(ven_Short_files);
for i = 1:nVShort
    load(ven_Short_files{i});
    lenVSthick = ASC.lenVSthick;
    V_short_tem = nan(1000,1);
    
    maxR = max(lenVSthick(:,2));
    lenVSthick(:,2) = lenVSthick(:,2)./maxR;
    
    lenVSthick(:,1) =  round(lenVSthick(:,1));
    [C,ia,ic] = unique(lenVSthick(:,1));
    lenVSthick = lenVSthick(ia,:);
    V_short_tem = nan(1000,1);
    if lenVSthick(end,1) > 1000
        lenVSthick(lenVSthick(:,1)>1000,:) = [];
    end
    
    % fill the value gaps which filled with NaNs.
    V_tem(lenVSthick(:,1))=lenVSthick(:,2);
    zeroIdx = find(V_tem==0);
    valueIdx = find(V_tem~=0);
    % aveVolume(1:valueIdx-1) = aveVolume(valueIdx(1));
    for j = 1:length(zeroIdx)
        if zeroIdx(j) < valueIdx(1)
            V_tem(zeroIdx(j)) = V_tem(valueIdx(1));
        else
            idx_tem = find(valueIdx < zeroIdx(j));
            V_tem(zeroIdx(j)) = V_tem(valueIdx( idx_tem(end) ));
        end
    end
    
    V_short_tem(1:length(V_tem)) = V_tem';
    V_short(:,i) = V_short_tem;
    
    clear V_tem
end
[ave_V_short,sd_V_short,se_V_short] = imean(V_short,2);
idx_isnan = ~isnan(ave_V_short);
idx = find(idx_isnan == 1); idx(end) = 690; % delete the last several noise  points.
ave_V_short(idx(end)+1:end) = [];
sd_V_short(idx(end)+1:end) = [];
se_V_short(idx(end)+1:end) = [];
ave_V_short(isnan(ave_V_short)) = ave_V_short(idx(1));
sd_V_short(isnan(sd_V_short)) = sd_V_short(idx(1));
se_V_short(isnan(se_V_short)) = 0;


% ================plot the figure===============
ave_V_long = smooth(ave_V_long,5);
ave_V_short = smooth(ave_V_short,5);

ave_V_long = ave_V_long(1:600);
se_V_long = se_V_long(1:600);
ave_V_short = ave_V_short(1:164);
se_V_short = se_V_short(1:164);

figure;
hold on;
patchX_vens = [1:length(ave_V_short),length(ave_V_short):-1:1];
patchY_vens = [ave_V_short+se_V_short ;flipud(ave_V_short-se_V_short)]';
patch(patchX_vens,patchY_vens,[255,165,0]/255,'EdgeColor','none','FaceAlpha',0.5);
patchX_venl = [1:length(ave_V_long),length(ave_V_long):-1:1];
patchY_venl = [ave_V_long+se_V_long ;flipud(ave_V_long-se_V_long)]';
patch(patchX_venl,patchY_venl,[255,56,56]/255,'EdgeColor','none','FaceAlpha',0.5);
plot(1:length(ave_V_short),ave_V_short,'Color','#FFA500');
plot(1:length(ave_V_long),ave_V_long,'Color','#FF3838');
set(gca,'xlim',[0,700],'ylim',[0,1],'yTick',0:0.25:1);
legend('VENs','VENl');



% =====================Testing codes ================
% figure; hold on;
% plot(1:length(ave_V_short),ave_V_short,'Color',[1,71/255,209/255]);
% plot(1:length(ave_V_long),ave_V_long,'Color',[1,56/255,56/255]);
% 
% 
% figure;
% for i = 1:nVShort
%     load(ven_Short_files{i});
%     lenVSthick = ASC.lenVSthick;
%     plot(lenVSthick(:,1),lenVSthick(:,2));
% end




% load('./dataMorph/asc&mat/s647');
% figure;plotASC(ASC);

% VEN_filename = 'dataMorph/asc&mat/n201.mat';
% load(VEN_filename);
% [ ASC,h ] = basalBigDenThickAna( ASC,1 );
