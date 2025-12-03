
clear all
% close all

AIS_D = [0,50,70,90,110,130,170];
stimPeriod = [2020,21980];

path = '';
fileName = {};
% fileName{end+1} = './mechs/savedData_PC_soma-1-hill-30.dat';
fileName{end+1} = './data/VEN_spikes_vs_AIS_D/savedData_VENL_soma_1.dat';
fileName{end+1} = './data/VEN_spikes_vs_AIS_D/savedData_VENL_bigDend0-1.dat';
fileName{end+1} = './data/VEN_spikes_vs_AIS_D/savedData_VENL_bigDend3-0.5.dat';
fileName{end+1} = './data/VEN_spikes_vs_AIS_D/savedData_VENL_bigDend3-1.dat';
fileName{end+1} = './data/VEN_spikes_vs_AIS_D/savedData_VENL_bigDend5-0.2.dat';
fileName{end+1} = './data/VEN_spikes_vs_AIS_D/savedData_VENL_bigDend5-0.5.dat';
fileName{end+1} = './data/VEN_spikes_vs_AIS_D/savedData_VENL_bigDend5-1.dat';

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
% color = {[0.9,0.9,0.9],[0.8,0.8,0.8],[0.7,0.7,0.7],[0.5,0.5,0.5],[0.3,0.3,0.3],[0.1,0.1,0.1]};
figure;hold on;
for i = 1:nFile
    v = spike5{i}.spikeTrace;
    dvdt = diff(v)*100;   % the data sample interval is 0.01ms, here change the time unit to ms.
    plot(v(2:end),dvdt,'color',color{i});
end
plot([-70,-50],[100,100],'k-','LineWidth',2);  % plot the scale bar.
plot([-50,-50],[100,150],'k-','LineWidth',2);
legend({'AIS D: 0 um', 'AIS D: 50 um', 'AIS D: 70 um','AIS D: 90 um','AIS D: 110 um', 'AIS D: 130 um', 'AIS D: 170 um',})

% the following points values were gotten from the above pahseplot figure by ginput.
points =[
  -34.8291   65.1854
    8.2910  193.0480
  -36.3015   60.8233
    1.6652  167.1483
  -17.7914  128.9804
    0.4032  129.7983
  -25.0482  102.2628
   -0.2279   98.4460
  -34.1981   62.4591
   -3.6985   60.8233
  -38.1946   37.6499
  -11.1656   30.5616
  -38.9308   22.9280
  -20.4207   12.5682];
for i = 1:7
%     phaseHeightIdx(i) = ( points(2*i,2)-points(2*(i-1)+1,2) )/( points(2*(i-1)+1,2)+points(2*i,2) );
    phaseHeightIdx(i) = ( points(2*i,2)-points(2*(i-1)+1,2) )/( points(2*i,2) );
end
figure; hold on;
plot(AIS_D,phaseHeightIdx,'o','markerFaceColor','r','markerEdgeColor','k');
[t]=polyfit(AIS_D,phaseHeightIdx,1);
plot(AIS_D,polyval(t,AIS_D),'LineStyle','-','Color','k');
set(gca,'xlim',[-20,200],'ylim',[-1,1],'xTick',0:50:200,'yTick',-1:0.5:1);
% statistic
stat = LinearModel.fit(AIS_D,phaseHeightIdx);
p_slope = stat.Coefficients.pValue('x1');


% allign the first and last spikes and compare their amps and widthes.
d1(:,1) = spike1{1}.spikeTrace(1310:1310+1500);
d1(:,2) = spike1{2}.spikeTrace(1355:1355+1500);
d1(:,3) = spike1{3}.spikeTrace(1385:1385+1500);
d1(:,4) = spike1{4}.spikeTrace(1420:1420+1500);
d1(:,5) = spike1{5}.spikeTrace(1510:1510+1500);
d1(:,6) = spike1{6}.spikeTrace(1450:1450+1500);
d1(:,7) = spike1{7}.spikeTrace(1635:1635+1500)-2;

d5(:,1) = spike5{1}.spikeTrace(1520:1520+1500);
d5(:,2) = spike5{2}.spikeTrace(1540:1540+1500);
d5(:,3) = spike5{3}.spikeTrace(1552:1552+1500);
d5(:,4) = spike5{4}.spikeTrace(1570:1570+1500);
d5(:,5) = spike5{5}.spikeTrace(1630:1630+1500)-1;
d5(:,6) = spike5{6}.spikeTrace(1455:1455+1500)-2;
d5(:,7) = spike5{7}.spikeTrace(1415:1415+1500)-3;

figure('position',[300,500,1000,400]); 
subplot(1,2,1); 
hold on;
x = (1:size(d1,1))./100; % change the data sample points to time (ms). because the sample point interval is 0.01 ms.
for i =1:7
    plot(x,d1(:,i),'color',color{i});
end
plot([6,8],[-20,-20],'k-','LineWidth',2);  % plot the scale bar.
plot([8,8],[-20,0],'k-','LineWidth',2);
set(gca,'xlim',[0,10],'ylim',[-80,60]);
legend({'AIS D: 0 um', 'AIS D: 50 um', 'AIS D: 70 um', 'AIS D: 90 um', 'AIS D: 110 um', 'AIS D: 130 um', 'AIS D: 170 um'})
% statistic

subplot(1,2,2);
hold on;
for i =1:7
    plot(x,d5(:,i),'color',color{i});
end
set(gca,'xlim',[0,10],'ylim',[-80,60]);
legend({'AIS D: 0 um', 'AIS D: 50 um', 'AIS D: 70 um', 'AIS D: 90 um', 'AIS D: 110 um', 'AIS D: 130 um', 'AIS D: 170 um'})

% calculate the spike widthes and heights.
heights = max(d1)-d1(1,:);
halfHeights = heights(i)/2 + d1(1,:);
for i = 1:nFile
    [~,peakIdx] = max(d1(:,i));
    spike_leftArm = d1(1:peakIdx,i);
    spike_rightArm = d1(peakIdx:end,i);
    [~,leftIdx] = min(abs(spike_leftArm - halfHeights(i) ));
    [~,rightIdx] = min(abs(spike_rightArm - halfHeights(i) ));
    spikeWidthes(i) = abs(rightIdx+peakIdx-leftIdx)/100;  % change the unit to ms, because the sample frequency is 100 per ms.
end

figure;
subplot(1,2,1);
plot(AIS_D,spikeWidthes,'Marker','o','MarkerSize',8,'MarkerEdgeColor',[255,56,56]/255,'LineStyle','none');
hold on;
[t]=polyfit(AIS_D,spikeWidthes,1);
plot(AIS_D,polyval(t,AIS_D),'LineStyle','-','Color','k');
ax = gca();
set(ax,'xlim',[-30,200],'ylim',[0,4]); box off;
title('VEN spike width to AIS D');
ax.XLabel.String = 'AIS D (um)';
ax.YLabel.String = 'Spike width (ms)';
set(ax,'xTick',0:30:180,'yTick',0:1:14);
% statistic
stat = LinearModel.fit(AIS_D,spikeWidthes);
p_slope = stat.Coefficients.pValue('x1');
if p_slope < 0.05
    fprintf('The slope of the fitting line for spike widthes is significant from horizontal. p = %6.4f. \n',p_slope);
else
    fprintf('the fitting line for spike widthes is NOT significant different from horizontal.p = %6.4f. \n', p_slope);
end

subplot(1,2,2);
plot(AIS_D,heights,'Marker','o','MarkerSize',8,'MarkerEdgeColor',[255,56,56]/255,'MarkerFaceColor','none','LineStyle','none');
hold on;
[t]=polyfit(AIS_D,heights,1);
plot(AIS_D,polyval(t,AIS_D),'LineStyle','-','Color','k');
ax = gca();
set(ax,'xlim',[-30,200],'ylim',[20,120]); box off;
title('VEN spike amp to AIS D');
ax.XLabel.String = 'AIS D (um)';
ax.YLabel.String = 'Spike amplitudes (mV)';
set(ax,'xTick',0:30:180,'yTick',0:20:160);
% statistic
stat = LinearModel.fit(AIS_D,heights);
p_slope = stat.Coefficients.pValue('x1');
if p_slope < 0.05
    fprintf('The slope of the fitting line for spike amplitudes is significant from horizontal. p = %6.4f. \n',p_slope);
else
    fprintf('the fitting line for spike amplitudes is NOT significant different from horizontal.p = %6.4f. \n', p_slope);
end







% --------------plot the dvdt rising slope----------------.
figure('position',[500,200,600,1000]);
subplot(2,1,1);
hold on;
for i = 1:nFile
    v = spike1{i}.spikeTrace;
    dvdt = diff(v)*100;   % the data sample interval is 0.01ms, here change the time unit to ms.
    
    [idx1,~] = find(dvdt>5);
    [idx2,~] = find(v<-42);
    idx = intersect(idx1,idx2);
    v_t = v(idx);dvdt_t = dvdt(idx);
    plot(v_t,dvdt_t,'-','color',color{i});
    
    dvdt_t_half = (max(dvdt_t)+min(dvdt_t))/2;
    [~,idx_half] = min(abs(dvdt_t-dvdt_t_half));
    v_half = v_t(idx_half);
    slope_half(i) = (dvdt_t(idx_half+5)-dvdt_t(idx_half-5))/(v_t(idx_half+5)-v_t(idx_half-5));
    b = dvdt_t(idx_half)-slope_half(i)*v_half;
    if i == 1 || i == nFile
        if i==1; x = [-47,-45]; end
        if i==nFile; x = [-45,-43]; end
        y = slope_half(i).*x + b;
        plot(x,y,'k-','LineWidth',2);
    end
end
ax=gca();
title('Insertion of spike phase rising slopes');
ax.XLabel.String = 'Voltage (mV)';
ax.YLabel.String = 'dv/dt (mV/t)';
set(ax,'xlim',[-50,-40],'ylim',[5,40]);
set(ax,'xTick',-50:5:-35,'yTick',5:10:35);
box off;

subplot(2,1,2);
hold on;
plot(AIS_D,slope_half,'Marker','o','MarkerSize',8,'MarkerEdgeColor','r','LineStyle','none')
hold on;
[t]=polyfit(AIS_D,slope_half,1);
plot(AIS_D,polyval(t,AIS_D),'LineStyle','-','Color','k');
ax = gca();
set(ax,'xlim',[0,170],'ylim',[0,15]); box off;
title('VEN spike phase rising slope relative to AIS L');
ax.XLabel.String = 'AIS D (um)';
ax.YLabel.String = 'Rising slope';
set(ax,'xTick',0:20:80,'yTick',0:3:15);
% statistic
stat = LinearModel.fit(AIS_D,slope_half);
p_slope = stat.Coefficients.pValue('x1');
if p_slope < 0.05
    fprintf('The slope of the fitting line for spike amps is significant from horizontal. p = %6.4f. \n',p_slope);
else
    fprintf('the fitting line for spike amps is NOT significant different from horizontal.p = %6.4f. \n', p_slope);
end








