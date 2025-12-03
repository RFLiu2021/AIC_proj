
clear all
close all

stimPeriod = [2020,21980];

path = '';
fileName = {};
fileName{end+1} = './data/PC_vs_VEN_spikes/savedData_PC26deg.dat';
fileName{end+1} = './data/PC_vs_VEN_spikes/savedData_VENL_bigDend5-0.2.dat';

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


% =================plot the simulated spike trains===================
colors = {[21/255,21/255,211/255],[1,56/255,56/255]};
t = (1:size(d,1) )/100;   % change the data sample points to time (ms). because the sample point interval is 0.01 ms.
figure('position',[500,800,1000,350]);
subplot(1,2,1); hold on;
hold on;
plot(t,d(:,1),'Color',colors{1});
plot(t,d(:,2),'Color',colors{2});
title("Simulated example traces");
xlabel("Time (ms)");
ylabel("Membrane potential (mV)");
legend({"PC","VEN"});



% allign the first and last spikes and compare their amps and widthes.
d1(:,1) = spike1{1}.spikeTrace(1100:1100+1500);
d1(:,2) = spike1{2}.spikeTrace(1517:1517+1500);
 
subplot(1,2,2); hold on;
x = (1:size(d1,1))./100; % change the data sample points to time (ms). because the sample point interval is 0.01 ms.
for i =1:2
    plot(x,d1(:,i),'color',colors{i});
end
set(gca,'xlim',[0,8],'ylim',[-80,60]);
title("Example APs");
legend({'PC', 'VEN'});
xlabel("Time (ms)");
ylabel("Membrane potential (mV)");