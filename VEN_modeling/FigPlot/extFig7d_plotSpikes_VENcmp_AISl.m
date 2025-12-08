
clear all
close all

AIS_L = 20:10:70;
stimPeriod = [2020,21980];

path = '';
fileName = {};
% fileName{end+1} = './data/VEN_spikes_vs_AIS_L/savedData_VENL__bigDend5_0.5_AIS5.dat';
% fileName{end+1} = './data/VEN_spikes_vs_AIS_L/savedData_VENL__bigDend5_0.5_AIS10.dat';
fileName{end+1} = './data/VEN_spikes_vs_AIS_L/savedData_VENL__bigDend5_0.5_AIS20.dat';
fileName{end+1} = './data/VEN_spikes_vs_AIS_L/savedData_VENL__bigDend5_0.5_AIS30.dat';
fileName{end+1} = './data/VEN_spikes_vs_AIS_L/savedData_VENL__bigDend5_0.5_AIS40.dat';
fileName{end+1} = './data/VEN_spikes_vs_AIS_L/savedData_VENL__bigDend5_0.5_AIS50.dat';
fileName{end+1} = './data/VEN_spikes_vs_AIS_L/savedData_VENL__bigDend5_0.5_AIS60.dat';
fileName{end+1} = './data/VEN_spikes_vs_AIS_L/savedData_VENL__bigDend5_0.5_AIS70.dat';

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
        idx = traceBoarder+idx-30;
        file{i}.amp = peakV-d_tem(idx);
        file{i}.threhsholdIdx = idx;
        

    file{i}.peakIdx = peakIdx;
    file{i}.peakV = peakV;
    
    spike1{i}.spikeTrace = d_tem(traceBoarder(1):traceBoarder(2));

end

color = {[250,218,201]/255,[247,198,173]/255,[241,154,190]/255,[255,71,209]/255,[255,98,115]/255,[218,35,55]/255,[203,27,10]/255, [81,0,27]/255};

% ----------------------------------------------------------------------
figure('position',[800,200,600,1000]);
subplot(2,1,1);
hold on;
for i = 1:nFile
    v = spike1{i}.spikeTrace(1000:5000);
    dvdt = diff(v)*100;   % the data sample interval is 0.01ms, here change the time unit to ms.
    plot(v(2:end),dvdt,'color',color{i});
end
plot([-70,-50],[100,100],'k-','LineWidth',2);  % plot the scale bar.
plot([-50,-50],[100,140],'k-','LineWidth',2);
legend({ 'AIS L: 20 um', 'AIS L: 30 um', 'AIS L: 40 um', 'AIS L: 50 um', 'AIS L: 60 um', 'AIS L: 70 um'})
ax=gca();
title('Spike phase shape relative to AIS L');
ax.XLabel.String = 'Voltage (mV)';
ax.YLabel.String = 'dv/dt (mV/t)';
set(gca,'xlim',[-80,40],'ylim',[-40,160]);
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

subplot(2,1,2); %figure(); 
hold on;
x = (1:size(peakTrace,1))./100; % change the data sample points to time (ms). because the sample point interval is 0.01 ms.
for i =1:nFile
    plot(x,peakTrace(:,i),'color',color{i});
end
% plot([0,0],[-50,20],'k');
plot([6,8],[-40,-40],'k-','LineWidth',2);  % plot the scale bar.
plot([8,8],[-40,-10],'k-','LineWidth',2);
legend({ 'AIS L: 20 um', 'AIS L: 30 um', 'AIS L: 40 um', 'AIS L: 50 um', 'AIS L: 60 um', 'AIS L: 70 um'})
ax=gca();
title('Spike phase shape depends on AIS L');
ax.XLabel.String = 'Time (ms)';
ax.YLabel.String = 'Voltage (mV)';
set(gca,'xlim',[0,10],'ylim',[-70,40]);
set(ax,'xTick',0:2:10,'yTick',-60:20:40);
box off;

% --------------plot the dvdt rising slope----------------.
figure('position',[500,200,600,1000]);
subplot(2,1,1);
hold on;
for i = 1:nFile
    v = spike1{i}.spikeTrace;
    dvdt = diff(v)*100;   % the data sample interval is 0.01ms, here change the time unit to ms.
    
    [idx1,~] = find(dvdt>6);
    [idx2,~] = find(v<-38);
    idx = intersect(idx1,idx2);
    v_t = v(idx);dvdt_t = dvdt(idx);
    plot(v_t,dvdt_t,'-','color',color{i});
    
    dvdt_t_half = (max(dvdt_t)+min(dvdt_t))/2;
    [~,idx_half] = min(abs(dvdt_t-dvdt_t_half));
    v_half = v_t(idx_half);
    slope_half(i) = (dvdt_t(idx_half+5)-dvdt_t(idx_half-5))/(v_t(idx_half+5)-v_t(idx_half-5));
    b = dvdt_t(idx_half)-slope_half(i)*v_half;
    if i == 1 || i == nFile
        if i==1; x = [-43,-40]; end
        if i==nFile; x = [-47,-45]; end
        y = slope_half(i).*x + b;
        plot(x,y,'k-','LineWidth',2);
    end
end
ax=gca();
title('Insertion of spike phase rising slopes');
ax.XLabel.String = 'Voltage (mV)';
ax.YLabel.String = 'dv/dt (mV/t)';
set(ax,'xlim',[-50,-35],'ylim',[5,40]);
set(ax,'xTick',-50:5:-35,'yTick',5:10:35);
box off;

subplot(2,1,2);
hold on;
plot(AIS_L,slope_half,'Marker','o','MarkerSize',8,'MarkerEdgeColor','r','LineStyle','none')
hold on;
[t]=polyfit(AIS_L,slope_half,1);
plot(AIS_L,polyval(t,AIS_L),'LineStyle','-','Color','k');
ax = gca();
set(ax,'xlim',[0,90],'ylim',[0,15]); box off;
title('VEN spike phase rising slope relative to AIS L');
ax.XLabel.String = 'AIS L (um)';
ax.YLabel.String = 'Rising slope';
set(ax,'xTick',0:20:80,'yTick',0:3:15);
% statistic
stat = LinearModel.fit(AIS_L,slope_half);
p_slope = stat.Coefficients.pValue('x1');
if p_slope < 0.05
    fprintf('The slope of the fitting line for spike phase AIS part is significant from horizontal. p = %6.4f. \n',p_slope);
else
    fprintf('the fitting line for spike amps is NOT significant different from horizontal.p = %6.4f. \n', p_slope);
end


% -----------------------------------------------------------------
% calculate and plot comparision of the amp. 
for i = 1:nFile
    amp(i) = max( peakTrace(:,i) )-aveThreshold;
end
figure('position',[600,200,600,1000]);
subplot(2,1,1);
plot(AIS_L,amp,'Marker','o','MarkerSize',8,'MarkerEdgeColor',[255,56,56]/255,'LineStyle','none')
hold on;
[t]=polyfit(AIS_L,amp,1);
plot(AIS_L,polyval(t,AIS_L),'LineStyle','-','Color','k');
% 这里尚缺少斜率的检验分析。
ax = gca();
set(ax,'xlim',[0,90],'ylim',[0,120]); box off;
title('VEN spike amplitude to AIS L');
ax.XLabel.String = 'AIS L (um)';
ax.YLabel.String = 'Spike amp (mV)';
set(ax,'xTick',0:20:90,'yTick',0:40:120);
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
    spikeWidth(i) = abs(rightIdx+peakIdx-leftIdx)/100;  % change the unit to ms, because the sample frequency is 100 per ms.
        
end
subplot(2,1,2); %figure;
plot(AIS_L,spikeWidth,'Marker','o','MarkerSize',8,'MarkerEdgeColor',[255,56,56]/255,'LineStyle','none')
hold on;
[t]=polyfit(AIS_L,spikeWidth,1);
plot(AIS_L,polyval(t,AIS_L),'LineStyle','-','Color','k');
ax = gca();
set(ax,'xlim',[0,90],'ylim',[0,4]); box off;
title('VEN spike width to AIS L');
ax.XLabel.String = 'AIS L (um)';
ax.YLabel.String = 'Spike width (ms)';
set(ax,'xTick',0:20:80,'yTick',0:1:4);
% statistic
stat = LinearModel.fit(AIS_L,spikeWidth);
p_slope = stat.Coefficients.pValue('x1');
if p_slope < 0.05
    fprintf('The slope of the fitting line for spike widthes is significant from horizontal. p = %6.4f. \n',p_slope);
else
    fprintf('the fitting line for spike widthes is NOT significant different from horizontal.p = %6.4f. \n', p_slope);
end

