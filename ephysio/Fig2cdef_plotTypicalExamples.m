% clear all;
% close all;


dataPath = './data_ephy/seqCellData';

% % ----------testing code------------
% fileName = fullfile(dataPath, strcat('n202','.mat'));
% plot3Traces( fileName, 1 );
% set(gca,'ylim',[-150,100]);
% % 
% figure;
% load(fileName);
% plot(m_FP.alldt{1}(:,[2,12,17]));
%% plot example fire patterns.
% Example neurons list ...
% exampleNeurons = {'n181','n146','n151','s647','B725','n180',...
%     'n185','n200','n209','s476','s486','n159',...
%     'n173','n202','n12','s475','B735','n145'};%,'n154'};
% cellTypes = {'L23 IT HPCAL1','IT RORB','L5 IT FEZF2','VEN','L5ET','L6PC'};%,'L56CT'};

exampleNeurons = {'n197','n146','n151','s647','B725','n180',...
    'n178','B561','n209','s476','s486','n159',...
    'n173','n202','n12','s475','B735','n145'};%,'n154'};
cellTypes = {'L23 IT HPCAL1','IT RORB','L5 IT FEZF2','VEN','L5ET','L6PC'};%,'L56CT'};

f = figure;
nCellType = length(exampleNeurons);
for i = 1:nCellType
    fileName = fullfile(dataPath, strcat(exampleNeurons{i},'.mat'));
    h_tem = subplot(3,6,i);
    h_tem = plot3Traces_4fig2( fileName, 1, h_tem );
    set(h_tem,'xlim',[0,1],'ylim',[-150,100])
    if i < 7
        title(cellTypes{i});
    end
end



%% plot typical example Tau features.
% largeTauNeuron = 's174';
% smallTauNeuron = 's438';
exampleNeurons = {'s661','s438'};

f=figure('position',[600,500,300,750]);
for i = 1:2
    fileName = fullfile(dataPath, strcat(exampleNeurons{i},'.mat'));
    load(fileName);
    
    si = 1/m_FP.fsample*1e+6;
    d=m_FP.alldt{1}(:,1);
    d_tem = d;
    d_tem = d_tem(5900:15000);  % get the data during stimulation.
    d_tem = smooth(d_tem,7);
    idx_min = find(d_tem==min(d_tem));
    d_tem = d_tem(1:idx_min(1)); % get the data from stimulation onset to the minimal data point in first 1/3 parts.
    t = (1:length(d_tem))*si*1e-3;   % unit is ms. if you want the finally tau unit is sec, here use 1e-6 instead of 1e-3. Because teh unit of si is uSec.

    beta = tauFit(t',d_tem);
    tau = beta(2);
    idx_fitData = [ 5900, 5900+length(d_tem) ];
    fit1.beta1 = beta;
    fit1.idx_data = idx_fitData;
    [ d_fit1 ] = expFunc( beta, t );
    fit1.d_fit = d_fit1;
       
    t = (1:length(d))*si*1e-3;
    h_tem = subplot(2,1,i); hold on
    plot(t,d,'k');  % plot trace.
    plot( (fit1.idx_data(1):fit1.idx_data(2)-1)*si*1e-3, fit1.d_fit, 'r','LineStyle' , '-','LineWidth', 2 );
    plot([250,350;350,350],[-56,-56;-56,-51],'k');
    text(200,-55,strcat('Tau: ',num2str(tau)) );
    set(h_tem,'xlim',[90,400]);     

    hold off;
end

%% plot typical example sag ratio feature.
% largeSagRatioNeuron = 's486';
% smallSagRatioNeuron = 's544';
exampleNeurons = {'s544','s486'};
deSampleFactor = 20;

f=figure('position',[600,500,300,750]);
for i = 1:2
    fileName = fullfile(dataPath, strcat(exampleNeurons{i},'.mat'));
    load(fileName);
    sagRatio = Attr(1).Sag_index;
    
    si = 1/m_FP.fsample*1e+6;
    d=m_FP.alldt{1}(:,1);
    d = filter_50Lowpass(d,m_FP.fsample);
    d_tem = decimate(d,deSampleFactor);
    t = (1:length(d_tem))*si*deSampleFactor*1e-3;   % unit is ms. if you want the finally tau unit is sec, here use 1e-6 instead of 1e-3. Because teh unit of si is uSec.
       
    t = (1:length(d))*si*1e-3;
    h_tem = subplot(2,1,i); hold on
    plot(t,d,'k');  % plot trace.
    plot([400,500;500,500],[-70,-70;-70,-60],'k');
    text(200,-55,strcat('Sag ratio: ',num2str(sagRatio)) );
    set(h_tem,'xlim',[0,600],'ylim',[-100,-50]);     

    hold off;
end

%% plot typical example AP rise-decay ratio features.
exampleNeurons = {'B569', 's212'};

f=figure('position',[600,500,300,750]);
% plot the first neuron.
fileName = fullfile(dataPath, strcat(exampleNeurons{1},'.mat'));
load(fileName);  
sweep1 = Attr.idx_sweep_rhoebase;
d = m_FP.alldt{1}(:,sweep1);
ind_threshold = Attr.APproperties(sweep1).ind_threshold(1);
d1 = d(18205:18254);
d2 = d(18255:18576);
t1 = (1:length(d1) )/length(d1);%*si;    % unit: ms.
t2 = ((1:length(d2)) +length(d1)+10 )/length(d1);%*si;

subplot(2,1,1); hold on;
plot(t1,d1,'k');
plot(t2,d2,'k');
% plot([5000,7000;7000,7000],[-20,-20;-20,20],'k');
set(gca,'xlim',[0,8]);

% plot the second neuron.
fileName = fullfile(dataPath, strcat(exampleNeurons{2},'.mat'));
load(fileName);  
sweep1 = Attr.idx_sweep_rhoebase;
d = m_FP.alldt{1}(:,sweep1);
ind_threshold = Attr.APproperties(sweep1).ind_threshold(1);
d1 = d(17380:17453);
d2 = d(17454:17714);
t1 = (1:length(d1) )/length(d1);%*si;    % unit: ms.
t2 = ((1:length(d2)) +length(d1)+10 )/length(d1);%*si;

subplot(2,1,2); hold on;
plot(t1,d1,'k');
plot(t2,d2,'k');
% plot([5000,7000;7000,7000],[-20,-20;-20,20],'k');
set(gca,'xlim',[0,8]);


%% plot typical example Rm features.
exampleNeurons = {'n173','n181'};%'s438'};

figure('position',[600,500,300,500]);
% plot the first neuron.
fileName = fullfile(dataPath, strcat(exampleNeurons{1},'.mat'));
load(fileName);  
si = 1/m_FP.fsample*1e+3;
d1 = m_FP.alldt{1}(:,1);
d1 = decimate(d1,10);
t1 = (1:length(d1) )*si*10;    % unit: ms.
subplot(2,1,1); hold on;
plot(t1,d1,'k');
set(gca,'xlim',[0,1000],'ylim',[-120,-20]);

fileName = fullfile(dataPath, strcat(exampleNeurons{2},'.mat'));
load(fileName);  
si = 1/m_FP.fsample*1e+3;
d2 = m_FP.alldt{1}(:,1);
d2 = decimate(d2,10);
t2 = (1:length(d2) )*si*10;    % unit: ms.
subplot(2,1,2); hold on;
plot(t2,d2,'k');
plot([750,950;950,950],[-120,-120;-120,-80],'k');
set(gca,'xlim',[0,1000],'ylim',[-160,-60]);


%% plot typical example AP threhsold features
exampleNeurons = {'n192', 'n147'};

f=figure('position',[600,500,300,500]);
% plot the first neuron.
fileName = fullfile(dataPath, strcat(exampleNeurons{1},'.mat'));
load(fileName);  
sweep1 = Attr.idx_sweep_rhoebase;
d1 = m_FP.alldt{1}(:,sweep1);
d1 = d1(5000:5500);
idx1 = ( 1:length(d1) )+500;
ind_threshold1 = Attr.APproperties(sweep1).AP_threshold(1);
fileName = fullfile(dataPath, strcat(exampleNeurons{2},'.mat'));
load(fileName);  
sweep1 = Attr.idx_sweep_rhoebase;
d2 = m_FP.alldt{1}(:,sweep1);
d2 = d2(5900:6400);
idx2 = 1:length(d2);
ind_threshold2 = Attr.APproperties(sweep1).AP_threshold(1);


plot(idx1,d1,'k'); hold on;
plot([500,750],[-35.6889,-35.6889],'r');
plot(idx2,d2,'k');
plot([0,250],[-47.3311,-47.3311],'r');
plot([900,1100;1100,1100],[-30,-30;-30,20],'k');
% set(gca,'xlim',[0,1000],'ylim',[-120,-20]);



%% plot typical example latency changing slope feature
exampleNeurons = {'s212', 's213'};

f=figure;
% plot the first neuron.
fileName = fullfile(dataPath, strcat(exampleNeurons{1},'.mat'));
load(fileName); 
d = m_FP.alldt{1};
d1 = d(1:20000,14);d2 = d(1:15000,15);d3 = d(1:12500,16);d4 = d(1:11500,17);d5 = d(1:10000,18);
subplot(2,1,1); hold on;
plot(d1,'k');plot(d2,'k');plot(d3,'k');plot(d4,'k');plot(d5,'k');
plot([23000,28000;28000,28000],[-60,-60;-60,-10],'k');
set(gca,'xlim',[5000,30000]);
hold off;


fileName = fullfile(dataPath, strcat(exampleNeurons{2},'.mat'));
load(fileName); 
dd = m_FP.alldt{1};
dd1 = dd(1:30000,15);dd2 = dd(1:20000,16);dd3 = dd(1:15000,17);dd4 = dd(1:13000,18);dd5 = dd(1:12500,19);
subplot(2,1,2); hold on;
plot(dd1,'k');plot(dd2,'k');plot(dd3,'k');plot(dd4,'k');plot(dd5,'k');
set(gca,'xlim',[5000,30000]);
hold off;


%% plot typical example firing rate changing slope feature
exampleNeurons = {'n193', 'n209'};

f=figure;
% plot the first neuron.
fileName = fullfile(dataPath, strcat(exampleNeurons{1},'.mat'));
load(fileName); 
clear FR I_in
FR(1) = Attr.APproperties(12).numSpikes/0.6; % unit: spikes/s.
FR(2) = Attr.APproperties(13).numSpikes/0.6;
FR(3) = Attr.APproperties(14).numSpikes/0.6;
FR(4) = Attr.APproperties(15).numSpikes/0.6;
FR(5) = Attr.APproperties(16).numSpikes/0.6;
FR(6) = Attr.APproperties(17).numSpikes/0.6;
I_in = (12:17)*20-220;    % Unit: pA
P_FR = polyfit(I_in,FR,1);   % linear fit.

subplot(2,1,1); hold on;
plot(I_in,FR,'marker','o','markerFace','k','markerSize',10,'LineStyle','none');
plot(I_in,P_FR(1).*I_in+P_FR(2),'k-','lineWidth',1.2);
set(gca,'xlim',[0,140],'ylim',[0,20]);

fileName = fullfile(dataPath, strcat(exampleNeurons{2},'.mat'));
load(fileName); 
clear FR I_in
FR(1) = Attr.APproperties(14).numSpikes/0.6; % unit: spikes/s.
FR(2) = Attr.APproperties(15).numSpikes/0.6;
FR(3) = Attr.APproperties(16).numSpikes/0.6;
FR(4) = Attr.APproperties(17).numSpikes/0.6;
FR(5) = Attr.APproperties(18).numSpikes/0.6;
FR(6) = Attr.APproperties(19).numSpikes/0.6;
I_in = (14:19)*20-220;    % Unit: pA
P_FR = polyfit(I_in,FR,1);   % linear fit.

subplot(2,1,2); hold on;
plot(I_in,FR,'marker','o','markerFace','k','markerSize',10,'LineStyle','none');
plot(I_in,P_FR(1).*I_in+P_FR(2),'k-','lineWidth',1.2);
plot([100,140],[15,15],'k');
set(gca,'xlim',[40,180],'ylim',[0,20]);




