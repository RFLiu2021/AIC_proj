
clear all
close all

AIS_L = 20:10:70;
stimPeriod = [2020,21980];

path = '';
fileName = {};
fileName{end+1} = './data/PC_spikes_vsAIS_L/savedData_PC_soma_AIS20.dat';
fileName{end+1} = './data/PC_spikes_vsAIS_L/savedData_PC_soma_AIS30.dat';
fileName{end+1} = './data/PC_spikes_vsAIS_L/savedData_PC_soma_AIS40.dat';
fileName{end+1} = './data/PC_spikes_vsAIS_L/savedData_PC_soma_AIS50.dat';
fileName{end+1} = './data/PC_spikes_vsAIS_L/savedData_PC_soma_AIS60.dat';
fileName{end+1} = './data/PC_spikes_vsAIS_L/savedData_PC_soma_AIS70.dat';

nFile = length(fileName);
% extract data from the files.
d = nan(25000,5);
for i = 1:nFile
    fid = fopen( fileName{i},'r');
    d_tem = fscanf(fid,'%f',[1,inf]);
    d_tem = d_tem(3:end)';
    d(:,i) = d_tem;  % remove the first 2 numbers which not the trace data
    fclose(fid);
    [N,peakIdx,peakV] = findSpikes_model(d_tem,-20); % find spikes.
    
    if N == 1
        interPoints(1) = round(mean(stimPeriod));
        traceBoarder = int16( [stimPeriod(1),interPoints,stimPeriod(2)] );
    elseif N > 1
        for j = 1:N-1
            interPoints(j) = round( peakIdx(j)+peakIdx(j+1) )/2;
        end
        traceBoarder = int16( [stimPeriod(1),interPoints,stimPeriod(2)] );
    end
    

        trace_tem = d_tem(traceBoarder(1):traceBoarder(2));
        d3 = diff(trace_tem,3);
        [~,idx] = max(d3); 
        idx = traceBoarder(1)+idx-30;
        file{i}.amp = peakV-d_tem(idx);
        file{i}.threhsholdIdx = idx;
        

    file{i}.peakIdx = peakIdx;
    file{i}.peakV = peakV;
    
    spike1{i}.spikeTrace = d_tem(traceBoarder(1):traceBoarder(2));

end

color = {[207,229,227]/255,[184,199,209]/255,[149,177,212]/255,[123,161,168]/255,[107,155,184]/255,[88,178,220]/255,[0,166,198]/255, [21,21,211]/255,[36,22,84]/255};
figure('position',[500,300,600,1000]);subplot(2,1,1);
hold on;
for i = 1:nFile
    v = spike1{i}.spikeTrace;
    dvdt = diff(v)*100;   % the data sample interval is 0.01ms, here change the time unit to ms.
    plot(v(2:end),dvdt,'color',color{i});
end
plot([-70,-50],[100,100],'k-','LineWidth',2);  % plot the scale bar.
plot([-50,-50],[100,140],'k-','LineWidth',2);
legend({ 'AIS L: 20 um', 'AIS L: 30 um', 'AIS L: 40 um', 'AIS L: 50 um', 'AIS L: 60 um', 'AIS L: 70 um'})
set(gca,'xlim',[-80,40],'ylim',[-60,160]);


% allign the first and last spikes and compare their amps and widthes.
for i = 1:nFile
    [~,peakIdx] = max(spike1{i}.spikeTrace);
    peakTrace(:,i) = spike1{i}.spikeTrace( peakIdx-180:peakIdx+1180 );
end
aveThreshold = mean(peakTrace(1,:));
for i = 1:nFile
    peakTrace(:,i) = peakTrace(:,i)-(peakTrace(1,i)-aveThreshold);
end

subplot(2,1,2);   %figure(); 
hold on;
x = (1:size(peakTrace,1))./100; % change the data sample points to time (ms). because the sample point interval is 0.01 ms.
for i =1:nFile
    plot(x,peakTrace(:,i),'color',color{i});
end
% plot([0,0],[-50,20],'k');
plot([6,8],[-40,-40],'k-','LineWidth',2);  % plot the scale bar.
plot([8,8],[-40,-10],'k-','LineWidth',2);
set(gca,'xlim',[0,10],'ylim',[-80,40]);
legend({ 'AIS L: 20 um', 'AIS L: 30 um', 'AIS L: 40 um', 'AIS L: 50 um', 'AIS L: 60 um', 'AIS L: 70 um'})

% calculate and plot comparision of the amp. 
for i = 1:nFile
    amp(i) = max( peakTrace(:,i) )-aveThreshold;
end
figure('position',[500,300,600,1000]);subplot(2,1,1);
plot(AIS_L,amp,'Marker','o','MarkerSize',8,'MarkerEdgeColor','b','LineStyle','none')
hold on;
[t]=polyfit(AIS_L,amp,1);
plot(AIS_L,polyval(t,AIS_L),'LineStyle','-','Color','k');
set(gca,'xlim',[0,90],'ylim',[0,120]); box off;
title('PC spike amp relative to AIS L');
% statistic
stat = LinearModel.fit(AIS_L,amp);
p_slope = stat.Coefficients.pValue('x1');
if p_slope < 0.05
    fprintf('The slope of the fitting line for spike amps is significant from horizontal. p = %6.4f. \n',p_slope);
else
    fprintf('the fitting line for spike amps is NOT significant different from horizontal.p = %6.4f. \n', p_slope);
end


% calculate and plot comparision of the spike width. 
for i = 1:nFile
    halfHeight = amp(i)/2 + aveThreshold;
    [~,peakIdx] = max(peakTrace(:,i));
    spike_leftArm = peakTrace(1:peakIdx,i);
    spike_rightArm = peakTrace(peakIdx:end,i);
    [~,leftIdx] = min(abs(spike_leftArm - halfHeight ));
    [~,rightIdx] = min(abs(spike_rightArm - halfHeight ));
    spikeWidth(i) = abs(rightIdx+peakIdx-leftIdx)/100;
        
end
subplot(2,1,2); %figure;
plot(AIS_L,spikeWidth,'Marker','o','MarkerSize',12,'MarkerEdgeColor','b','LineStyle','none')
hold on;
[t]=polyfit(AIS_L,spikeWidth,1);
plot(AIS_L,polyval(t,AIS_L),'LineStyle','-','Color','k');
set(gca,'xlim',[0,90],'ylim',[0,4]); box off;
title('PC spike width relative to AIS L');
set(gca,'xTick',0:20:80,'yTick',0:1:4);
% statistic
stat = LinearModel.fit(AIS_L,spikeWidth);
p_slope = stat.Coefficients.pValue('x1');
if p_slope < 0.05
    fprintf('The slope of the fitting line for spike widthes is significant from horizontal. p = %6.4f. \n',p_slope);
else
    fprintf('the fitting line for spike widthes is NOT significant different from horizontal.p = %6.4f. \n', p_slope);
end










