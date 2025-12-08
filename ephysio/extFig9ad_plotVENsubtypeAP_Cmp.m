clear all;
close all;

colors = {[255,165,0]/255;[255,56,56]/255};
%% ===========plot the example VEN-L and VEN-S =================
VENS_file = './data_ephy/VENsubtypeCmp_extFig8ab/exampleVENS.mat';  
VENL_file = './data_ephy/VENsubtypeCmp_extFig8ab/exampleVENL.mat';


load(VENL_file);
Attr_VENL = Attr;
fsample = m_FP.fsample;
si_venl = 1e+6/m_FP.fsample;
d =  m_FP.alldt{1};
sweep_venl_1 = Attr(1).idx_sweep_rhoebase;
sweep_venl_2 = Attr(1).idx_sweep_maxFR;
d2_venl = m_FP.alldt{1}(6000:35500,sweep_venl_2);
idx = Attr(1).APproperties(sweep_venl_1).ind_threshold(1);
firstAP_data_venl = d(idx+5:idx+1005,sweep_venl_1);
d_venl = firstAP_data_venl;
firstAP_data_venl = firstAP_data_venl - d(idx,sweep_venl_1); % normalize the spike threshold to 0.
t_venl = si_venl*(1:length(firstAP_data_venl))/1000;

load(VENS_file);
Attr_VENS = Attr;
fsample = m_FP.fsample;
si_vens = 1e+6/m_FP.fsample;
d =  m_FP.alldt{1};
sweep_vens_1 = Attr(1).idx_sweep_rhoebase;
sweep_vens_2 = Attr(1).idx_sweep_maxFR;
d2_vens = m_FP.alldt{1}(6000:35500,sweep_vens_2);
idx = Attr(1).APproperties(sweep_vens_1).ind_threshold(1);
firstAP_data_vens = d(idx:idx+1000,sweep_vens_1);
d_vens = firstAP_data_vens;
firstAP_data_vens = filter_1000Lowpass(d_vens,fsample);
firstAP_data_vens = firstAP_data_vens - d(idx,sweep_vens_1); % normalize the spike threshold to 0.
t_vens = si_vens*(1:length(firstAP_data_vens))/1000;


figure('Position',[500,318,300,250]);
hold on; 
plot(t_venl,firstAP_data_venl,'color',colors{2},'lineWidth',2);
plot(t_vens,firstAP_data_vens,'color',colors{1},'lineWidth',2);
plot([7,9;9,9,],[20,20;20,50],'color','k');  % plot the scale bar.
text(11,24,'4 ms');
text(12.5,43,'30 mV');
set(gca,'xlim',[0,10],'xTick',[5,10,15],'ylim',[-20,100],'box','off');
xlabel('Time (ms)');
ylabel('Voltage (mV)');
legend('VEN-L','VEN-S');


%% plot the spike phase plane.
figure('Position',[500,318,700,200]);
subplot(1,2,1);
% figure();
numSpikes = length(Attr_VENS.APproperties(sweep_vens_2).phaseTrace);
% plot the spike phase plane.
for i = 2:numSpikes-1
    hold on;
    mV = Attr_VENS.APproperties(sweep_vens_2).phaseTrace(i).mV;
    dv_dt = Attr_VENS.APproperties(sweep_vens_2).phaseTrace(i).dv_dt;
    plot( mV, smooth(dv_dt,7),'color',colors{1})
end
colorbar();
xlabel('Amplitude (mV)');
ylabel('dV/dt (mV/ms)');
title('Example VEN-S');
set(gca,'xlim',[-50,80],'xTick',-40:40:80,'ylim',[-20,100],'yTick',[-20,0:25:100]);
hold off;

subplot(1,2,2);
% figure();
numSpikes = length(Attr_VENL.APproperties(sweep_venl_2).phaseTrace);
% plot the spike phase plane.
for i = 2:numSpikes-1
    hold on;
    mV = Attr_VENL.APproperties(sweep_venl_2).phaseTrace(i).mV;
    dv_dt = Attr_VENL.APproperties(sweep_venl_2).phaseTrace(i).dv_dt;
    plot( mV, smooth(dv_dt,7),'color',colors{2})
end
colorbar();
xlabel('Amplitude (mV)');
ylabel('dV/dt (mV/ms)');
title('Example VEN-L');
set(gca,'xlim',[-50,80],'xTick',-40:40:80,'ylim',[-20,100],'yTick',[-20,0:25:100]);
hold off;





%% plot the population comparison of AP width and amplitude.

load('./data_ephy/VENsubtypeCmp_extFig9ab/VENsubtype_APs.mat');

% plot the boxplot style.
figure('Position',[500,318,450,350]);
s(1)=subplot(1,2,1);
hold on; 
nMax = max( [ length(VENL.spike_halfWidth),length(VENS.spike_halfWidth) ] );
Y = nan(nMax,2);
Y(1:numel(VENL.spike_halfWidth),1) = reshape(VENL.spike_halfWidth,[numel(VENL.spike_halfWidth),1]);
Y(1:numel(VENS.spike_halfWidth),2) = reshape(VENS.spike_halfWidth,[numel(VENS.spike_halfWidth),1]);
scatter(ones(nMax,1), Y(:,1),50, 'marker','o','markerEdgeColor',[255,56,56]/255,'jitter','on','jitterAmount', 0.25,'LineWidth',1);
scatter(2*ones(nMax,1), Y(:,2),50, 'marker','o','markerEdgeColor',[255,71,209]/255,'jitter','on','jitterAmount', 0.25,'LineWidth',1);
h=boxplot(Y,'Notch','off','Labels',{'VEN-L','VEN-S'},'whisker',2,'symbol','','colors','k');
set(h,{'linew'},{2});
set(gca,'xlim',[0,3],'ylim',[-0,8],'xTick',[1,2],'yTick',0:2:10,'xTickLabel',{'VENL','VENS'},'box','off');
title(s(1),'Half Spike Width: VENL vs VENS');
ylabel('AP width (ms)');
legend('VENL','VENS')
hold off;

s(2)=subplot(1,2,2);
hold on; 
nMax = max( [ length(VENL.AP_Amplitude_mV),length(VENS.AP_Amplitude_mV) ] );
Y = nan(nMax,2);
Y(1:numel(VENL.AP_Amplitude_mV),1) = reshape(VENL.AP_Amplitude_mV,[numel(VENL.AP_Amplitude_mV),1]);
Y(1:numel(VENS.AP_Amplitude_mV),2) = reshape(VENS.AP_Amplitude_mV,[numel(VENS.AP_Amplitude_mV),1]);
scatter(ones(nMax,1), Y(:,1),50, 'marker','o','markerEdgeColor',[255,56,56]/255,'jitter','on','jitterAmount', 0.25,'LineWidth',1);
scatter(2*ones(nMax,1), Y(:,2),50, 'marker','o','markerEdgeColor',[255,71,209]/255,'jitter','on','jitterAmount', 0.25,'LineWidth',1);
h=boxplot(Y,'Notch','off','Labels',{'VEN-L','VEN-S'},'whisker',2,'symbol','','colors','k');
set(h,{'linew'},{2});
set(gca,'xlim',[0,3],'ylim',[-0,150],'xTick',[1,2],'xTickLabel',{'VENL','VENS'},'box','off');
ylabel('AP amp (mV)');
title(s(2),'AP amp: VENL vs VENS');
legend('VENL','VENS')
hold off;








