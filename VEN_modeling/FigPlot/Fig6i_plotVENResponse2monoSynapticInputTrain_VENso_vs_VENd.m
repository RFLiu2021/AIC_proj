clear all;
close all;

xlsFile = './data/axon2soma_vs_axon2Dend_thickDend_monoSynapticTrain_60random/firingRateMatrix_axonLocationCmp.xlsx';

[num,str,raw] = xlsread(xlsFile);

x = num(1,[3,6,9,12,15,18,21,24,27,30,33,36]);    % distance of synaptic input location to the soma, Unit:um
y = num(4:9,1);        % synaptic input freqency.
z_ven = num(4:9,[3,6,9,12,15,18,21,24,27,30,33,36]);  % data_VENthickDend, firing rate. unit: spikes/s
z_pc = num(4:9,[4,7,10,13,16,19,22,25,28,31,34,37]);  % data_PCthinDend, firing rate. unit: spikes/s.

[X,Y]=meshgrid(x,y);

[Xq,Yq] = meshgrid(0:150,0:1200);
Zq_ven = interp2(X,Y,z_ven,Xq,Yq);
figure;h_surf = surf(Xq,Yq,Zq_ven,'cdata',abs(Zq_ven));
h_surf.EdgeColor = 'none';
h_surf.FaceColor = 'flat';
h_surf.FaceAlpha=0.8;
temp1=caxis;colorbar
h_ax = gca;
h_ax.XLabel.String = 'Distance from synaptic input location to soma (um)';
h_ax.YLabel.String = 'Synaptic input frequency (Hz)';
h_ax.ZLabel.String = 'Firing rate (spikes/s)';
% hold on; plot3(x(6)*ones(1,6),y,z_ven(:,6),'k-','LineWidth',2);
set(h_ax,'xlim',[0,150],'ylim',[0,1200],'zlim',[0,40]);


figure();
Zq_pc = interp2(X,Y,z_pc,Xq,Yq);
h_surf = surf(Xq,Yq,Zq_pc);%,'cdata',abs(Zq-300));
h_surf.EdgeColor = 'none';
h_surf.FaceColor = 'flat';
h_surf.FaceAlpha=0.8;

h_ax = gca;
h_ax.XLabel.String = 'Distance from synaptic input location to soma (um)';
h_ax.YLabel.String = 'Synaptic input frequency (Hz)';
h_ax.ZLabel.String = 'Firing rate (spikes/s)';
colorbar();
% hold on; plot3(x(6)*ones(1,6),y,z_pc(:,6),'k-','LineWidth',2);
set(h_ax,'xlim',[0,150],'ylim',[0,1200],'zlim',[0,40]);
caxis(temp1);


figure;plot3(x(4)*ones(1,6),y,z_ven(:,4),'k-','LineWidth',2);
h_ax = gca;
h_ax.XLabel.String = 'Distance from synaptic input location to soma (um)';
h_ax.YLabel.String = 'Synaptic input frequency (Hz)';
h_ax.ZLabel.String = 'Firing rate (spikes/s)';
set(h_ax,'xlim',[0,150],'ylim',[0,1200],'zlim',[0,40]);


%% plot typical traces for comparision. 
% savedData_VENL_soma_axon2BigDend[3](0.7)_monoTrain_500Hz VS. savedData_VENL_soma_monoSynaptic_ThickDend_axon2Soma_500Hz

trace1 = './data/axon2soma_vs_axon2Dend_thickDend_monoSynapticTrain_60random/synapticLocation_Dend[3](0.7)/savedData_VENL_soma_axon2BigDend[3](0.7)_monoTrain_500Hz.dat';
trace2 = './data/thinDend_vs_bigDend_monoSynapticInputTrain_60random/synapticLocation_Dend[3](0.7)/savedData_VENL_soma_monoSynaptic_ThickDend_axon2Soma_500Hz.dat';

fid_a = fopen(trace1,'r');
d_1 = fscanf(fid_a,'%f',[1,inf]);
d_1 = d_1(3:end);
fclose(fid_a);
fid_b = fopen(trace2,'r');
d_2 = fscanf(fid_b,'%f',[1,inf]);
d_2 = d_2(3:end);
fclose(fid_b);

sampleFreq = 100;
t = ( 1:length(d_1) )/sampleFreq;
figure;
% plot the response of VEN whose axon connect to bigDend[1](1), normally.
subplot(4,1,1);
plot(t,d_1,'color',[255,56,56]/255);
text( 600, 25,'D: um, L 25 um, 500 Hz');
h_ax = gca;
h_ax.XLabel.String = 'Time (ms)';
h_ax.YLabel.String = 'Membrane potentials (mV)';
set(h_ax,'xlim',[0,1050],'ylim',[-100,50]);
% plot the response of VEN whose axon connect to bigDend[3](1)
subplot(4,1,2);
plot(t,d_2,'color',[255,71,209]/255);
hold on; plot([700,900;900,900],[-40,-40;-40,10],'k');   % plot scale bars.
text( 600, 25,'D: um, L 25 um, 500 Hz');
h_ax = gca;
h_ax.XLabel.String = 'Time (ms)';
h_ax.YLabel.String = 'Membrane potentials (mV)';
set(h_ax,'xlim',[0,1050],'ylim',[-100,50]);
% plot synaptic stim trace.
% synapticTrace: savedData_synapInputTrain_RecSynapSiteDend_500Hz, rec from dend[3](1), synaptic input upon BigDend[3](1)
synapticTrace = './data/synapticInputTrain/savedData_synapInputTrain_RecSynapSiteDend_500Hz.dat';
fid = fopen(synapticTrace,'r');
d_synap = fscanf(fid,'%f',[1,inf]);
d_synap = d_synap(3:end);
fclose(fid);
subplot(4,1,3);
t1 = ( 1:length(d_synap) )/400;
plot(t1,d_synap,'g');
h_ax = gca;
set(h_ax,'xlim',[0,1050],'ylim',[-75,-40]);
subplot(4,1,4);
t1 = ( 1:length(d_synap) )/400;
plot(t1,d_synap,'g');
h_ax = gca;
set(h_ax,'xlim',[500,600],'ylim',[-75,-40]);



% 
% 
% synapticTrace = './data/synapticInputTrain/savedData_synapInputTrain_RecSynapSiteDend_900Hz.dat';
% % synapticTrace = './mechs/savedData_VENL_soma.dat';
% fid = fopen(synapticTrace,'r');
% d_synap = fscanf(fid,'%f',[1,inf]);
% d_synap = d_synap(3:end);
% fclose(fid);
% t = ( 1:length(d_synap) )/sampleFreq;
% subplot(3,1,3);
% plot(t,d_synap,'k');
% 










