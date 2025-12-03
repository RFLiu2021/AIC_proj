clear all


load('./data_ephy/NaCurrent_Fig6d.mat');
cellNames = cellNames;
nData = maxNaCurrent;
NaCurr_VEN = exampleCellTraces.NaCurr_VEN;
NaCurr_PC = exampleCellTraces.NaCurr_PC;
 
%% =============plot example traces =================

nTraces = 17;   % test pulses from -80 to 0 mV.
stimStartT = 20; % Unit: ms.
dataRange = 20; % unit ms;
sampleFreq = 50000;
si_pc = 1000/1e+5;
si_ven = 1000/sampleFreq;


% snap the data from -2 to 10 ms relative to the test pulse starting.
dataIdx = int32([40-2,40+dataRange]./si_ven); 
v_VEN = imean(NaCurr_VEN(dataIdx(1):dataIdx(2),1:nTraces,:),3);
dataIdx = int32([stimStartT-2,stimStartT+dataRange]./si_pc); 
v_PC = imean(NaCurr_PC(dataIdx(1):dataIdx(2),1:nTraces,:),3);
t_ven = [-99:size(v_VEN,1)-100].*si_ven;   % time. unit:ms
t_pc = [-199:size(v_PC,1)-200].*si_pc;   % time. unit:ms

% normalize the base voltage and remove the fluctuations.
baseVEN = imean(v_VEN( 1:int32(0.002*sampleFreq),:),1);
for i = 1:length(baseVEN)
    v_VEN(:,i) = v_VEN(:,i)-baseVEN(i);
end
basePC = imean(v_PC( 1:int32(0.002*sampleFreq),:),1);
for i = 1:length(basePC)
    v_PC(:,i) = v_PC(:,i)-basePC(i);
end

figure('color','w','position',[600,500,500,600]);
subplot(2,1,1); hold on;
plot( t_ven, v_VEN(:,1:nTraces),'color',[255,56,56]/255);
line([5,8;8,8],[-30,-30;-30,-10],'color','k'); % scale bar
set(gca,'xlim',[-1,15],'ylim',[-40,10],'box','off');
xlabel('Time (ms)');
ylabel('Na current (pA)');
text(8,5,'VEN');

subplot(2,1,2);
plot( t_pc, v_PC(:,1:nTraces),'color',[21,21,211]/255);
set(gca,'xlim',[-1,15],'ylim',[-40,10]);
xlabel('Time (ms)');
ylabel('Na current (pA)');
text(8,5,'PC');
box('off');



%% =============plot population analysis boxplot =================


d_PC = abs( nData(nData(:,1)==1,2) );
d_VEN = abs( nData(nData(:,1)==2,2) );

[h,p]=ttest2(d_PC,d_VEN)
[p,h,stat]=ranksum(d_PC,d_VEN)

nMax = max( [ length(d_PC),length(d_VEN) ] );
Y = nan(nMax,2);
Y(1:numel(d_VEN),1) = reshape(d_VEN,[numel(d_VEN),1]);
Y(1:numel(d_PC),2) = reshape(d_PC,[numel(d_PC),1]);

% plot the figure;
figure('position',[600,400,300,450]); hold on;
scatter(ones(nMax,1), Y(:,1),150, 'marker','o','markerEdgeColor',[255,56,56]/255,'jitter','on','jitterAmount', 0.1,'LineWidth',1);
scatter(2*ones(nMax,1), Y(:,2),150, 'marker','o','markerEdgeColor',[21,21,211]/255,'jitter','on','jitterAmount', 0.1,'LineWidth',1);
h=boxplot(Y,'Notch','on','Labels',{'VEN','PC'},'whisker',2,'symbol','','colors','k');
set(h,{'linew'},{2});
legend('VEN','PC');
set(gca,'xlim',[0,3],'ylim',[0,60])




