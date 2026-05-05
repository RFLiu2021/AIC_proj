
clear all
close all

AIS_D = [10:10:70];
stimPeriod = [2020,21980];

path = '';
fileName = {};
fileName{end+1} = './data/PC_spikes_vs_AISd/savedData_PC_AIS30_Hillock10.dat';
fileName{end+1} = './data/PC_spikes_vs_AISd/savedData_PC_AIS30_Hillock20.dat';
fileName{end+1} = './data/PC_spikes_vs_AISd/savedData_PC_AIS30_Hillock30.dat';
fileName{end+1} = './data/PC_spikes_vs_AISd/savedData_PC_AIS30_Hillock40.dat';
fileName{end+1} = './data/PC_spikes_vs_AISd/savedData_PC_AIS30_Hillock50.dat';
fileName{end+1} = './data/PC_spikes_vs_AISd/savedData_PC_AIS30_Hillock60.dat';
fileName{end+1} = './data/PC_spikes_vs_AISd/savedData_PC_AIS30_Hillock70.dat';
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
% color = {[0.9,0.9,0.9],[0.8,0.8,0.8],[0.7,0.7,0.7],[0.5,0.5,0.5],[0.3,0.3,0.3],[0.1,0.1,0.1]};
figure('position',[500,200,600,1000]);subplot(2,1,1);
hold on;
for i = 1:nFile
    v = spike1{i}.spikeTrace;
    dvdt = diff(v)*100;   % the data sample interval is 0.01ms, here change the time unit to ms.
    plot(v(2:end),dvdt,'color',color{i});
end
plot([-70,-50],[100,100],'k-','LineWidth',2);  % plot the scale bar.
plot([-50,-50],[100,140],'k-','LineWidth',2);
legend({ 'AIS D: 10 um', 'AIS D: 20 um', 'AIS D: 30 um', 'AIS D: 40 um', 'AIS D: 50 um', 'AIS D: 60 um', 'AIS D: 70 um'})
ax=gca();
title('Spike phase shape depends on AIS D');
ax.XLabel.String = 'Voltage (mV)';
ax.YLabel.String = 'dv/dt (mV/t)';
set(gca,'xlim',[-80,40],'ylim',[-60,160]);
set(ax,'xTick',-80:20:80,'yTick',-40:40:160);
box off;


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
legend({ 'AIS D: 10 um', 'AIS D: 20 um', 'AIS D: 30 um', 'AIS D: 40 um', 'AIS D: 50 um', 'AIS D: 60 um', 'AIS D: 70 um'})
ax=gca();
title('Spike phase shape depends on AIS D');
ax.XLabel.String = 'Time (ms)';
ax.YLabel.String = 'Voltage (mV)';
set(gca,'xlim',[0,10],'ylim',[-80,40]);
set(ax,'xTick',0:2:10,'yTick',-80:40:40);
box off;

% --------------plot the dvdt rising slope----------------.
figure('position',[500,200,600,1000]);
subplot(2,1,1);
hold on;
for i = 1:nFile
    v = spike1{i}.spikeTrace;
    dvdt = diff(v)*100;   % the data sample interval is 0.01ms, here change the time unit to ms.
    
    [idx1,~] = find(dvdt>5);
    [idx2,~] = find(v<-20);
    idx = intersect(idx1,idx2);
    v_t = v(idx);dvdt_t = dvdt(idx);
    plot(v_t,dvdt_t,'-','color',color{i});
    
    [~,idx_25] = min(abs(dvdt_t-25));
    v_25 = v_t(idx_25);
    %slopes = diff(dvdt_t); 
    %slope_25(i) = slopes(idx_25);
    slope_25(i) = (dvdt_t(idx_25+5)-dvdt_t(idx_25-5))/(v_t(idx_25+5)-v_t(idx_25-5));
    b = dvdt_t(idx_25)-slope_25(i)*v_25;
    if i == 1 || i == nFile
        if i==1; x = [-35,-30]; end
        if i==nFile; x = [-42,-39]; end
        y = slope_25(i).*x + b;
        plot(x,y,'k-','LineWidth',2);
    end
end
ax=gca();
title('Insertion of spike phase rising slopes');
ax.XLabel.String = 'Voltage (mV)';
ax.YLabel.String = 'dv/dt (mV/t)';
set(ax,'xlim',[-50,-20],'ylim',[0,70]);
set(ax,'xTick',-50:10:-20,'yTick',0:20:60);
box off;

subplot(2,1,2);
hold on;
plot(AIS_D,slope_25,'Marker','o','MarkerSize',8,'MarkerEdgeColor','b','LineStyle','none')
hold on;
[t]=polyfit(AIS_D,slope_25,1);
plot(AIS_D,polyval(t,AIS_D),'LineStyle','-','Color','k');
ax = gca();
set(ax,'xlim',[0,80],'ylim',[0,8]); box off;
title('PC spike phase rising slope relative to AIS D');
ax.XLabel.String = 'AIS D (um)';
ax.YLabel.String = 'Rising slope';
set(ax,'xTick',0:20:80,'yTick',0:2:8);
% statistic
stat = LinearModel.fit(AIS_D,slope_25);
p_slope = stat.Coefficients.pValue('x1');
if p_slope < 0.05
    fprintf('The slope of the fitting line for spike amps is significant from horizontal. p = %6.4f. \n',p_slope);
else
    fprintf('the fitting line for spike amps is NOT significant different from horizontal.p = %6.4f. \n', p_slope);
end

% ---------------plot and compare the spikes' amps and widthes--------
% calculate and plot comparision of the amp. 
for i = 1:nFile
    amp(i) = max( peakTrace(:,i) )-aveThreshold;
end
figure('position',[500,200,600,1000]);subplot(2,1,1);
plot(AIS_D,amp,'Marker','o','MarkerSize',8,'MarkerEdgeColor','b','LineStyle','none')
hold on;
[t]=polyfit(AIS_D,amp,1);
plot(AIS_D,polyval(t,AIS_D),'LineStyle','-','Color','k');
set(gca,'xlim',[-10,90],'ylim',[0,120]); box off;
title('PC spike amp relative to AIS D');
% 这里尚缺少斜率的检验分析。
ax = gca();
set(gca,'xlim',[-10,90],'ylim',[0,120]); box off;
ax.XLabel.String = 'AIS D (um)';
ax.YLabel.String = 'Spike amplitude (mV)';
set(ax,'xTick',0:20:80,'yTick',0:40:120);
% statistic
stat = LinearModel.fit(AIS_D,amp);
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
subplot(2,1,2);  %figure;
plot(AIS_D,spikeWidth,'Marker','o','MarkerSize',8,'MarkerEdgeColor','b','LineStyle','none')
hold on;
[t]=polyfit(AIS_D,spikeWidth,1);
plot(AIS_D,polyval(t,AIS_D),'LineStyle','-','Color','k');
title('PC spike width relative to AIS D');
ax = gca();
set(gca,'xlim',[-10,90],'ylim',[0,4]); box off;
ax.XLabel.String = 'AIS D (um)';
ax.YLabel.String = 'Spike width (ms)';
set(ax,'xTick',0:20:80,'yTick',0:1:4);
% statistic
stat = LinearModel.fit(AIS_D,spikeWidth);
p_slope = stat.Coefficients.pValue('x1');
if p_slope < 0.05
    fprintf('The slope of the fitting line for spike widthes is significant from horizontal. p = %6.4f. \n',p_slope);
else
    fprintf('the fitting line for spike widthes is NOT significant different from horizontal.p = %6.4f. \n', p_slope);
end




