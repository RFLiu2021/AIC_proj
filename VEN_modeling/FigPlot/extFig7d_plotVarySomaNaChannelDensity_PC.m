

clear all;
close all;

dataPath = './data/PC_varyNaDensity';


fileName = {};
fileName{end+1} = 'savedData_PC_NaDensity_55S.dat';
fileName{end+1} = 'savedData_PC_NaDensity_60S.dat';
fileName{end+1} = 'savedData_PC_NaDensity_65S.dat';
fileName{end+1} = 'savedData_PC_NaDensity_70S.dat';
fileName{end+1} = 'savedData_PC_NaDensity_75S.dat';
fileName{end+1} = 'savedData_PC_NaDensity_80S.dat';


% plot the phase planes of all spikes.
si = 10;
NaDensity = [55,60,65,70,75,80];
stimPeriod = [2020,21980];
nFile = length(fileName);
for i = 1:nFile
    fid = fopen( fullfile(dataPath, fileName{i}),'r');
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
figure('position',[500,200,1000,800]);

% allign the first and last spikes and compare their amps and widthes.
for i = 1:nFile
    [~,peakIdx] = max(spike1{i}.spikeTrace);
    peakTrace(:,i) = spike1{i}.spikeTrace( peakIdx-120:peakIdx+1120 );
end
aveThreshold = mean(peakTrace(1,:));
for i = 1:nFile
    peakTrace(:,i) = peakTrace(:,i)-(peakTrace(1,i)-aveThreshold);
end

subplot(2,2,1);   %figure(); 
hold on;
x = (1:size(peakTrace,1))./100; % change the data sample points to time (ms). because the sample point interval is 0.01 ms.
for i =1:nFile
    plot(x,peakTrace(:,i),'color',color{i});
end
% plot([0,0],[-50,20],'k');
plot([6,8],[-40,-40],'k-','LineWidth',2);  % plot the scale bar.
plot([8,8],[-40,-10],'k-','LineWidth',2);
legend({ 'Na den: 55 S/cm2', 'Na den: 60 S/cm2', 'Na den: 65 S/cm2', 'Na den: 70 S/cm2', 'Na den: 75 S/cm2', 'Na den: 80 S/cm2'})
ax=gca();
title('Spike phase shape depends on somatic Na density');
ax.XLabel.String = 'Time (ms)';
ax.YLabel.String = 'Voltage (mV)';
set(gca,'xlim',[0,10],'ylim',[-90,50]);
set(ax,'xTick',0:2:10,'yTick',-50:50:50);
box off;


subplot(2,2,2);
hold on;
for i = 1:nFile
    v = spike1{i}.spikeTrace;
    dvdt = diff(v)*100;   % the data sample interval is 0.01ms, here change the time unit to ms.
    plot(v(2:end),dvdt,'color',color{i});
end
plot([-70,-50],[100,100],'k-','LineWidth',2);  % plot the scale bar.
plot([-50,-50],[100,140],'k-','LineWidth',2);
ax=gca();
title('Spike phase shape depends on somatic Na density');
ax.XLabel.String = 'Voltage (mV)';
ax.YLabel.String = 'dv/dt (mV/t)';
set(gca,'xlim',[-90,55],'ylim',[-60,150]);
set(ax,'xTick',-80:20:80,'yTick',-50:50:250);
box off;

% ---------------plot and compare the spikes' amps and widthes--------
% calculate and plot comparision of the amp. 
for i = 1:nFile
    amp(i) = max( peakTrace(:,i) )-aveThreshold;
end
subplot(2,2,3);
plot(NaDensity,amp,'Marker','o','MarkerSize',8,'MarkerEdgeColor','b','LineStyle','none')
hold on;
[t]=polyfit(NaDensity,amp,1);
plot(NaDensity,polyval(t,NaDensity),'LineStyle','-','Color','k');
title('PC spike amp relative to AIS D');
ax = gca();
set(gca,'xlim',[50,85],'ylim',[20,140]); box off;
ax.XLabel.String = 'Somatic membrane conductance (S/um2)';
ax.YLabel.String = 'Spike amplitude (mV)';
set(ax,'xTick',50:5:85,'yTick',20:40:140);
% statistic
stat = LinearModel.fit(NaDensity,amp);
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
    spikeWidth(i) = abs(rightIdx+peakIdx-leftIdx)/100;  % change the unit to ms by devided 100. Because of the sample frequency is 100 per ms.
        
end
subplot(2,2,4);  %figure;
plot(NaDensity,spikeWidth,'Marker','o','MarkerSize',8,'MarkerEdgeColor','b','LineStyle','none')
hold on;
[t]=polyfit(NaDensity,spikeWidth,1);
plot(NaDensity,polyval(t,NaDensity),'LineStyle','-','Color','k');
title('PC spike width relative to AIS D');
ax = gca();
set(gca,'xlim',[50,85],'ylim',[0,4]); box off;
ax.XLabel.String = 'Somatic membrane conductance (S/um2)';
ax.YLabel.String = 'Spike width (ms)';
set(gca,'xTick',50:5:85,'yTick',0:1:4);
% statistic
stat = LinearModel.fit(NaDensity,spikeWidth);
p_slope = stat.Coefficients.pValue('x1');
if p_slope < 0.05
    fprintf('The slope of the fitting line for spike widthes is significant from horizontal. p = %6.4f. \n',p_slope);
else
    fprintf('the fitting line for spike widthes is NOT significant different from horizontal.p = %6.4f. \n', p_slope);
end










