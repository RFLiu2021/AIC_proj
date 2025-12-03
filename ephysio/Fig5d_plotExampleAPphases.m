clear all;
close all;

matFilePath = './data_ephy/seqcellData/';
ExampleCell{1} = 'n181.mat';       %  ExampleCell_L23 HPCAL1, n193.mat
ExampleCell{2} = 'n105.mat';    %  ExampleCell_ITrorb 
ExampleCell{3} = 'n87.mat';      %  ExampleCell_L5IT FEZF2     s219,s533
ExampleCell{4} = 'B721';       %  ExampleCell_VEN   'VEN_M_0801_C10_0000_fp.mat',B571
ExampleCell{5} = 'B333';       %  ExampleCell_L5ET   'L5ET_M_0802_C11_0001_fp.mat'
ExampleCell{6} = 'n180.mat';       %  ExampleCell_L6PC 

neuronType =  {'L23IT','ITrorb','L5IT','VEN','L5ET','L6PC'};
% colors = {'#009900','#949C4E','#4666A6','#FF3838','#00CED1','#6FBC1E'};
colors = {[0,153,0]/255, [148,156,78]/255, [70,102,106]/255,[255,56,56]/255,[0,206,209]/255,[111,188,30]/255};

figure('Position',[600,500,1000,600]);

% plot an example neuron of  L2/3HAPCAL1.
fileName = strcat(matFilePath,ExampleCell{1});
load(fileName);
si = 1e+3/m_FP.fsample;
d = m_FP.alldt{1};  % extract the data. 
subplot(3,4,3); hold on;
d1 = d(:,1);
t = (1:length(d1))*si;
plot(t',smooth(d1,7),'k-');
plot(t',d(:,14),'Color','k');
plot(t',d(:,25),'Color',colors{1});
subplot(3,4,4);
mv = smooth(d(2600:17000,25),5);
dv_dt = diff(mv)/si;
plot(mv(2:end),dv_dt,'Color',colors{1},'LineStyle','-','LineWidth',1.0);
set(gca,'xlim',[-65,65],'ylim',[-50,300]);

% plot an example neuron of IT_rorb.
fileName = strcat(matFilePath,ExampleCell{2});
load(fileName);
si = 1e+3/m_FP.fsample;
d = m_FP.alldt{1};  % extract the data. 
subplot(3,4,5); hold on;
d1 = d(:,1);
t = (1:length(d1))*si;
plot(t',smooth(d1,7),'k-');
plot(t',smooth(d(:,15),7),'Color','k');
plot(t',smooth(d(:,20),7),'Color',colors{2});
subplot(3,4,6);
mv = filter_2000Lowpass(d(2700:17000,20),m_FP.fsample);
dv_dt = diff(mv)/si;
plot(mv(2:end),dv_dt,'Color',colors{2},'LineStyle','-','LineWidth',1.0);
set(gca,'xlim',[-65,65],'ylim',[-50,250]);

% plot an example neuron of L5IT_FEZF2.
fileName = strcat(matFilePath,ExampleCell{3});
load(fileName);
si = 1e+3/m_FP.fsample;
d = m_FP.alldt{1};  % extract the data. 
subplot(3,4,7); hold on;
d1 = d(:,2);
t = (1:length(d1))*si;
plot(t',smooth(d1,7),'k-');
plot(t',smooth(d(:,15),7),'Color','k');
plot(t',smooth(d(:,35),7),'Color',colors{3});
subplot(3,4,8);
% mv = smooth(d(6000:35000,40),7);
mv = filter_2000Lowpass(d(2700:17500,35),m_FP.fsample);
dv_dt = diff(mv)/si;
plot(mv(2:end),dv_dt,'Color',colors{3},'LineStyle','-','LineWidth',1.0);
set(gca,'xlim',[-60,65],'ylim',[-100,300]);

% % plot an example neuron of VEN.
% fileName = strcat(matFilePath,ExampleCell{4});
% load(fileName);
% si = 1e+3/m_FP.fsample;
% d = m_FP.alldt{1};  % extract the data. 
% subplot(3,4,3); hold on;
% d1 = d(:,2);
% t = (1:length(d1))*si;
% plot(t',smooth(d1,7),'k-');
% plot(t',smooth(d(:,13),7),'Color','k');
% plot(t',smooth(d(:,18),7),'Color',colors{4});
% subplot(3,4,4);
% mv = filter_2000Lowpass(d(3100:17500,18),m_FP.fsample);  % select the spikes traces for plotting.
% dv_dt = diff(mv)/si; % calculate the dv/dt.
% plot(mv(2:end),dv_dt,'Color',colors{4},'LineStyle','-','LineWidth',1.0);
% set(gca,'xlim',[-65,45],'ylim',[-50,200]);

% plot an example neuron of VEN.
fileName = strcat(matFilePath,ExampleCell{4});
load(fileName);
si = 1e+3/m_FP.fsample;
d = m_FP.alldt{1};  % extract the data. 
subplot(3,4,1); hold on;
d1 = d(:,1);
t = (1:length(d1))*si;
plot(t',smooth(d1,7),'k-');
plot(t',smooth(d(:,Attr(1).idx_sweep_rhoebase),7),'Color','k');  %15
plot(t',smooth(d(:,Attr(1).idx_sweep_maxFR),7),'Color',colors{4});  %23
subplot(3,4,2);
mv = smooth(d(7000:35000,23),5);  % select the spikes traces for plotting.
dv_dt = diff(mv)/si; % calculate the dv/dt.
plot(mv(2:end),dv_dt,'Color',colors{4},'LineStyle','-','LineWidth',1.0);
set(gca,'xlim',[-60,45],'ylim',[-50,150]);

% plot an example neuron of L5ET.
fileName = strcat(matFilePath,ExampleCell{5});
load(fileName);
si = 1e+3/m_FP.fsample;
d = m_FP.alldt{1};  % extract the data. 
subplot(3,4,9); hold on;
d1 = d(:,2);
t = (1:length(d1))*si;
plot(t',smooth(d1,7),'k-');
plot(t',smooth(d(:,Attr(1).idx_sweep_rhoebase),7),'Color','k');  
plot(t',smooth(d(:,Attr(1).idx_sweep_maxFR),7),'Color',colors{5}); 
% plot(t',smooth(d(:,12),7),'Color','k');
% plot(t',smooth(d(:,39),7),'Color',colors{5});
subplot(3,4,10);
mv = filter_2000Lowpass(d(3000:17500,Attr(1).idx_sweep_maxFR),m_FP.fsample);  % select the spikes traces for plotting.
% mv = smooth(d(3000:17500,Attr(1).idx_sweep_maxFR),5);
dv_dt = diff(mv)/si; % calculate the dv/dt.
plot(mv(2:end),dv_dt,'Color',colors{5},'LineStyle','-','LineWidth',1.0);
set(gca,'xlim',[-65,45],'ylim',[-70,200]);

% plot an example neuron of L6PC.
fileName = strcat(matFilePath,ExampleCell{6});
load(fileName);
si = 1e+3/m_FP.fsample;
d = m_FP.alldt{1};  % extract the data. 
subplot(3,4,11); hold on;
d1 = d(:,2);
t = (1:length(d1))*si;
plot(t',smooth(d1,7),'k-');
plot(t',smooth(d(:,14),7),'Color','k');
plot(t',smooth(d(:,18),7),'Color',colors{6});
subplot(3,4,12);
mv = filter_2000Lowpass(d(3000:17500,18),m_FP.fsample);  % select the spikes traces for plotting.
dv_dt = diff(mv)/si; % calculate the dv/dt.
plot(mv(2:end),dv_dt,'Color',colors{6},'LineStyle','-','LineWidth',1.0);
set(gca,'xlim',[-65,65],'ylim',[-70,250]);






