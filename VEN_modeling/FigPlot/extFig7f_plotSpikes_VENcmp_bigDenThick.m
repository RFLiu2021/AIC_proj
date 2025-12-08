
clear all
% close all

bigDenDiameter = [1,1.5,2,2.5,3,4,5];
stimPeriod = [2020,21980];

path = '';
fileName = {};
fileName{end+1} = './data/VEN_spikes_vs_bigDenThickVary/savedData_VENL_bigDen1.dat';
fileName{end+1} = './data/VEN_spikes_vs_bigDenThickVary/savedData_VENL_bigDen1_5.dat';
fileName{end+1} = './data/VEN_spikes_vs_bigDenThickVary/savedData_VENL_bigDen2.dat';
fileName{end+1} = './data/VEN_spikes_vs_bigDenThickVary/savedData_VENL_bigDen2_5.dat';
fileName{end+1} = './data/VEN_spikes_vs_bigDenThickVary/savedData_VENL_bigDen3.dat';
fileName{end+1} = './data/VEN_spikes_vs_bigDenThickVary/savedData_VENL_bigDen4.dat';
fileName{end+1} = './data/VEN_spikes_vs_bigDenThickVary/savedData_VENL_bigDen5.dat';

nFile = length(fileName);
% extract data from the files.
d = nan(25000,5);
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
        trace_tem = d_tem(traceBoarder(j):traceBoarder(j+1));
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

color = {[250,218,201]/255,[247,198,173]/255,[241,154,190]/255,[255,71,209]/255,[255,98,115]/255,[218,35,55]/255,[203,27,10]/255};

% allign the first and last spikes and compare their amps and widthes.
d1(:,1) = spike1{1}.spikeTrace(910:910+1100);
d1(:,2) = spike1{2}.spikeTrace(932:932+1100);
d1(:,3) = spike1{3}.spikeTrace(994:994+1100);
d1(:,4) = spike1{4}.spikeTrace(1069:1069+1100);
d1(:,5) = spike1{5}.spikeTrace(1152:1152+1100);
d1(:,6) = spike1{6}.spikeTrace(1272:1272+1100);
d1(:,7) = spike1{7}.spikeTrace(1614:1614+1100);

figure('position',[500,200,800,800]); 
ax=subplot(2,2,1);hold on;
x = (1:size(d1,1))./100; % change the data sample points to time (ms). because the sample point interval is 0.01 ms.
for i =1:nFile
    plot(x,d1(:,i),'color',color{i});
end
plot([6,8],[-20,-20],'k-','LineWidth',2);  % plot the scale bar.
plot([8,8],[-20,0],'k-','LineWidth',2);
set(gca,'xlim',[0,10],'ylim',[-80,60]);
ax.XLabel.String = 'Time (ms)';
ax.YLabel.String = 'Membrane potential (mV)';
legend({'Dendrite diameter: 1 um', '1.5 um', '2 um','2.5 um','3 um','4 um','5 um'})
title('Spike shape');

ax=subplot(2,2,2);hold on;
for i = 1:nFile
    v = spike1{i}.spikeTrace;
    dvdt = diff(v)*100;   % the data sample interval is 0.01ms, here change the time unit to ms.
    plot(v(2:end),dvdt,'color',color{i});
end
line([-1,29;29,29],[130,130;130,180],'color','k','LineWidth',2);
set(gca,'xlim',[-70,30],'ylim',[-80,200]);
ax.XLabel.String = 'Membrane potential (mV)';
ax.YLabel.String = 'dv/dt (mV/ms)';
title('AP phase plane');



% calculate the spike widthes and heights.
heights = max(d1)-d1(1,:);
halfHeights = heights(i)/2 + d1(1,:);
spikeWidthes = nan(nFile,1);
for i = 1:nFile
    [~,peakIdx] = max(d1(:,i));
    spike_leftArm = d1(1:peakIdx,i);
    spike_rightArm = d1(peakIdx:end,i);
    [~,leftIdx] = min(abs(spike_leftArm - halfHeights(i) ));
    [~,rightIdx] = min(abs(spike_rightArm - halfHeights(i) ));
    spikeWidthes(i) = abs(rightIdx+peakIdx-leftIdx)/100;  % change the unit to ms, because the sample frequency is 100 per ms.
end

subplot(2,2,3);
plot(bigDenDiameter,heights,'Marker','o','MarkerSize',8,'MarkerEdgeColor',[255,56,56]/255,'MarkerFaceColor','none','LineStyle','none');
hold on;
[t]=polyfit(bigDenDiameter,heights,1);
plot(bigDenDiameter,polyval(t,bigDenDiameter),'LineStyle','-','Color','k');
ax1 = gca();
set(ax1,'xlim',[0,6],'ylim',[0,140]); box off;
title('VEN spike amp to bigDen diameter');
ax1.XLabel.String = 'BigDen diameter (um)';
ax1.YLabel.String = 'Spike amplitudes (mV)';
set(ax1,'xTick',0:1:6,'yTick',0:20:140);
% statistic
stat = LinearModel.fit(bigDenDiameter,heights);
p_slope = stat.Coefficients.pValue('x1');
if p_slope < 0.05
    fprintf('The slope of the fitting line for spike amplitudes is significant from horizontal. p = %6.4f. \n',p_slope);
else
    fprintf('the fitting line for spike amplitudes is NOT significant different from horizontal.p = %6.4f. \n', p_slope);
end

subplot(2,2,4);
plot(bigDenDiameter,spikeWidthes,'Marker','o','MarkerSize',8,'MarkerEdgeColor',[255,56,56]/255,'LineStyle','none');
hold on;
[t]=polyfit(bigDenDiameter,spikeWidthes,1);
plot(bigDenDiameter,polyval(t,bigDenDiameter),'LineStyle','-','Color','k');
ax = gca();
set(ax,'xlim',[0,6],'ylim',[0,5]); box off;
title('VEN spike width to bigDen diameter');
ax.XLabel.String = 'BigDen diameter (um)';
ax.YLabel.String = 'Spike width (ms)';
set(ax,'xTick',0:1:6,'yTick',0:1:5);
% statistic
stat = LinearModel.fit(bigDenDiameter,spikeWidthes);
p_slope = stat.Coefficients.pValue('x1');
if p_slope < 0.05
    fprintf('The slope of the fitting line for spike widthes is significant from horizontal. p = %6.4f. \n',p_slope);
else
    fprintf('the fitting line for spike widthes is NOT significant different from horizontal.p = %6.4f. \n', p_slope);
end








