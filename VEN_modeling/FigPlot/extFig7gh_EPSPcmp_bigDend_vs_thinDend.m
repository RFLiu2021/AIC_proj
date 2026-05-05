clear all;
close all;

pathName = './data/thinDend_vs_bigDend_singleSynapticInput/';
fileName ={};
fileName{end+1} = 'savedData_VENL_synapticInputAtThinDend[0](1).dat';
fileName{end+1} = 'savedData_VENL_synapticInputAtThinDend[1](1).dat';
fileName{end+1} = 'savedData_VENL_synapticInputAtThinDend[3](0.15).dat';
fileName{end+1} = 'savedData_VENL_synapticInputAtThinDend[3](0.45).dat';
fileName{end+1} = 'savedData_VENL_synapticInputAtThinDend[3](1).dat';
fileName{end+1} = 'savedData_VENL_synapticInputAtThinDend[5](0.15).dat';
fileName{end+1} = 'savedData_VENL_synapticInputAtThinDend[5](0.55).dat';
fileName{end+1} = 'savedData_VENL_synapticInputAtThinDend[7](0).dat';
fileName{end+1} = 'savedData_VENL_synapticInputAtThinDend[7](0.5).dat';

fileName{end+1} = 'savedData_VENL_synapticInputAtbigDend[0](1).dat';
fileName{end+1} = 'savedData_VENL_synapticInputAtbigDend[1](1).dat';
fileName{end+1} = 'savedData_VENL_synapticInputAtbigDend[3](0.15).dat';
fileName{end+1} = 'savedData_VENL_synapticInputAtbigDend[3](0.45).dat';
fileName{end+1} = 'savedData_VENL_synapticInputAtbigDend[3](1).dat';
fileName{end+1} = 'savedData_VENL_synapticInputAtbigDend[5](0.15).dat';
fileName{end+1} = 'savedData_VENL_synapticInputAtbigDend[5](0.55).dat';
fileName{end+1} = 'savedData_VENL_synapticInputAtbigDend[7](0).dat';
fileName{end+1} = 'savedData_VENL_synapticInputAtbigDend[7](0.5).dat';

nFile = length(fileName);
d = nan(25000,18);
for i = 1:nFile
    fileName_t = strcat(pathName,fileName{i});
    fid = fopen( fileName_t,'r');
    d_tem = fscanf(fid,'%f',[1,inf]);
    d_tem = d_tem(3:end)';
    d(1:length(d_tem),i) = d_tem;  % remove the first 2 numbers which not the trace data
    fclose(fid);
 
    refV = mean( d(950:990,i) );
    [peak,idx] = max(d(:,i));
    amp(i) = peak - refV;
    delay_peak(i) = (idx-1000)/100;   % unit: ms.
    
    halfHeight = (refV+peak)/2;
    [v_t,idx_t] = min( abs(halfHeight-d(1001:idx,i)) );
    slope(i) = ( d_tem(1000+idx_t+2)-d_tem(1000+idx_t-2) )/(10*0.01);
end

% dist = [26.2,52.34,55.68,62.42,72.54,89.39,104.99,136.19,173.02,217.65];% distance of synaptic input location to the soma.
dist = [52.34,55.68,62.42,72.54,89.39,104.99,136.19,173.02,217.65];% distance of synaptic input location to the soma.
for i = 1:length(dist)
    legendStr{i} =  strcat( num2str(dist(i)),'um');
end
d = d(900:7000,:);
figure('Position',[600,600,300,600]);
x = 1:size(d,1);
x = x/100;
subplot(2,1,1);
color = [ [207,229,227]/255;[184,199,209]/255;[149,177,212]/255;[123,161,168]/255;[107,155,184]/255;[88,178,220]/255;[0,166,198]/255; [0,89,119]/255; [21,21,211]/255; [0,0,255]/255; ];%;[36,22,84]/255; ];[60,94,145]/255;
set(groot,'defaultAxesColorOrder',color)
plot(x,d(:,1:9))%,'color',color{i});%[116,54,50]/255);
ax1 = gca();
ax1.XLim = [0,30];
ax1.YLim = [-72,-64];
ax1.XLabel.String = 'Time (ms)';
ax1.YLabel.String = 'Voltage (mV)';
legend(legendStr);

title('EPSP evoked along thin dendrite');
subplot(2,1,2);
color = [[250,218,201]/255;[247,198,173]/255;[235,180,113]/255; [239,154,72]/255;[241,154,190]/255;[255,98,115]/255;[255,71,209]/255;[240,86,84]/255;[255,0,0]/255; [203,27,10]/255; ];
set(groot,'defaultAxesColorOrder',color)
plot(x,d(:,10:end));%,'color',[241,154,190]/255);
ax2 = gca();
ax2.XLim = [0,30];
ax2.YLim = [-72,-64];
title('EPSP evoked along thick dendrite');
ax2.XLabel.String = 'Time (ms)';
ax2.YLabel.String = 'Voltage (mV)';
legend(legendStr);

figure('position',[1000,600,500,450]);
plot(dist, amp(:,1:9),'color',[21,21,211]/255,'LineWidth',2,'marker','o','markerSize',8);hold on;
plot(dist, amp(:,10:end),'color',[255,56,56]/255,'LineWidth',2,'marker','o','markerSize',8);
ax = gca();
ax.XLim = [0,250];
ax.YLim = [0,6];
title( 'EPSP amplitude comparision of BigDend vs thinDend');
ax.XLabel.String = 'Synaptic input distance to Soma (um)';
ax.YLabel.String = 'EPSP amp (mV)';
legend('Thin Dendrite','Thick dendrite');
% paired t-test.
[h,p] = ttest(amp(:,1:9),amp(:,10:end))
% fitting and student's t-test for slopes
stat1 = LinearModel.fit(dist,amp(:,1:9));
p_slopeOFamp_thin = stat1.Coefficients.pValue('x1');
stat2 = LinearModel.fit(dist,amp(:,10:end));
p_slopeOFamp_thick = stat2.Coefficients.pValue('x1');






