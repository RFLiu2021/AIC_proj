
clear all
% close all

stimPeriod = [2020,21980];

path = '';
fileName = {};
fileName{end+1} = './data/PC_vs_VEN_spikes/savedData_PC26deg.dat';
fileName{end+1} = './data/PC_vs_VEN_spikes/savedData_VENL_bigDend5-0.5.dat';

nFile = length(fileName);
% extract data from the files.
d = nan(25000,2);
for i = 1:nFile
    fid = fopen( fileName{i},'r');
    d_tem = fscanf(fid,'%f',[1,inf]);
    d_tem = d_tem(3:end)';
    d(:,i) = d_tem;  % remove the first 2 numbers which not the trace data
    fclose(fid);
    [N,peakIdx,peakV] = findSpikes_model(d_tem,-40); % find spikes.
    
    for j = 1:N-1
        interPoints(j) = round( peakIdx(j)+peakIdx(j+1) )/2;
    end
    traceBoarder = int16( [stimPeriod(1),interPoints,stimPeriod(2)] );
    
    for j = 1:N-1
        trace_tem = d_tem(traceBoarder(i):traceBoarder(i+1));
        d3 = diff(trace_tem,3);
        [~,idx] = max(d3); 
        idx = traceBoarder(j)+idx-30;
        file{i}.amp(j) = peakV(j)-d_tem(idx);
        file{i}.threhsholdIdx(j) = idx;
        
    end
    file{i}.peakIdx = peakIdx;
    file{i}.peakV = peakV;
    
    spike1{i}.spikeTrace = d_tem(traceBoarder(1):traceBoarder(2));
    spike2{i}.spikeTrace = d_tem(traceBoarder(2):traceBoarder(3));
    spike3{i}.spikeTrace = d_tem(traceBoarder(3):traceBoarder(4));
    spike4{i}.spikeTrace = d_tem(traceBoarder(4):traceBoarder(5));
    spike5{i}.spikeTrace = d_tem(traceBoarder(5):traceBoarder(6));
    spike6{i}.spikeTrace = d_tem(traceBoarder(6):traceBoarder(7));
end

% ==================plot spike phase plane =========================
nSpike = 6;
figure('position',[500,800,1000,350]);
subplot(1,2,1);
color=colormap(spring(6));
hold on;
for i = 1:6
    eval( strcat('v=spike',num2str(i),'{2}.spikeTrace;') );
%     v = spike5{i}.spikeTrace;
    dvdt = diff(v)*100;   % the data sample interval is 0.01ms, here change the time unit to ms.
    plot(v(2:end),dvdt,'color',color(i,:));
end
plot([20,40],[150,150],'k-','LineWidth',2);      % plot the scale bar.
plot([40,40],[150,200],'k-','LineWidth',2);      % plot the scale bar.
set(gca,'xlim',[-80,60],'ylim',[-40,250]);
subplot(1,2,2);hold on;
color=colormap(summer(nSpike));
for i = 1:nSpike
    eval( strcat('v=spike',num2str(i),'{1}.spikeTrace;') );
%     v = spike5{i}.spikeTrace;
    dvdt = diff(v)*100;   % the data sample interval is 0.01ms, here change the time unit to ms.
    plot(v(2:end),dvdt,'color',color(i,:));
end
% set(gca,'xlim',[],'ylim',[]);


% ==================plot the zoomed-in spike shapes: PC vs VEN =========================

% allign the first and last spikes and compare their amps and widthes.
d1(:,1) = spike1{1}.spikeTrace(1100:1100+1500);
d1(:,2) = spike1{2}.spikeTrace(1350:1350+1500);

d5(:,1) = spike5{1}.spikeTrace(1570:1570+1500);
d5(:,2) = spike5{2}.spikeTrace(1500:1500+1500)-4;

figure('position',[300,500,1500,300]); 
color = {[21/255,21/255,211/255],[1,56/255,56/255]};
subplot(1,3,1); hold on;
x = (1:size(d1,1))./100; % change the data sample points to time (ms). because the sample point interval is 0.01 ms.
for i =1:2
    plot(x,d1(:,i),'color',color{i});
end
set(gca,'xlim',[0,8],'ylim',[-80,60]);
legend({'PC', 'VEN'});

% normalize and plot the spike heights.
peakIdx_1 = file{1}.peakIdx(1)-traceBoarder(1);
peakIdx_2 = file{2}.peakIdx(1)-traceBoarder(1);
data1 = spike1{1}.spikeTrace;
data2 = spike1{2}.spikeTrace;

data_1 = data1(peakIdx_1-300:peakIdx_1);
data_11 = data1(peakIdx_1+1:peakIdx_1+1000);
data_1 = (data_1-min(data_1))./(max(data_1)-min(data_1));
data_11 = (data_11-min(data_11))./(max(data_11)-min(data_11));

data_2 = data2(peakIdx_2-300:peakIdx_2);
data_22 = data2(peakIdx_2+1:peakIdx_2+1000);
data_2 = (data_2-min(data_2))./(max(data_2)-min(data_2));
data_22 = (data_22-min(data_22))./(max(data_22)-min(data_22));

data1 = nan(1310,1); data2 = nan(1310,1);
data1(1:301) = data_1;
data1(311:1310)=data_11;
data2(1:301) = data_2;
data2(311:1310)=data_22;
x1 = (1:size(data1,1))/100; % change the data sample points to time (ms). because the sample point interval is 0.01 ms.
subplot(1,3,2);
plot(x1,data1,'b');
hold on;
plot(x1,data2,'r');
set(gca,'xlim',[0,11],'ylim',[0,1]);

subplot(1,3,3);hold on;
for i =1:2
    plot(x,d1(:,i),'color',color{i},'Linestyle','-');
    plot(x,d5(:,i),'color',color{i},'Linestyle','--');
end
plot([5,7],[0,0],'k-','LineWidth',2);      % plot the scale bar.
plot([7,7],[0,20],'k-','LineWidth',2);      % plot the scale bar.
set(gca,'xlim',[0,8],'ylim',[-80,60]);


% for j = 1:6
%     figure; hold on;
%     for i = 1:5
%         idx = file{i}.threhsholdIdx(j)-traceBoarder(j);
%         d = spike1{i}.spikeTrace(1:end);
%         plot(d,'color',color{i});
%     end
% end










