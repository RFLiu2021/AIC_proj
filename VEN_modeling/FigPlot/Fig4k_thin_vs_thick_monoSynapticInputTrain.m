clear all;
close all;

pathName = './data/VEN_vs_PC_response2SynTrainOnBasal/';
% ----------------------------------------axon2soma-----------------------------
f_axon2Soma_thick_s ={};
f_axon2Soma_thick_s{end+1} = 'savedData_VENL_soma_monoSynaptic_ThickDend_axon2Soma_200Hz.dat';
f_axon2Soma_thick_s{end+1} = 'savedData_VENL_soma_monoSynaptic_ThickDend_axon2Soma_300Hz.dat';
f_axon2Soma_thick_s{end+1} = 'savedData_VENL_soma_monoSynaptic_ThickDend_axon2Soma_500Hz.dat';
f_axon2Soma_thick_s{end+1} = 'savedData_VENL_soma_monoSynaptic_ThickDend_axon2Soma_700Hz.dat';
f_axon2Soma_thick_s{end+1} = 'savedData_VENL_soma_monoSynaptic_ThickDend_axon2Soma_900Hz.dat';
f_axon2Soma_thick_s{end+1} = 'savedData_VENL_soma_monoSynaptic_ThickDend_axon2Soma_1100Hz.dat';
f_axon2Soma_thin_s ={};
f_axon2Soma_thin_s{end+1} = 'savedData_VENL_soma_monoSynaptic_ThinDend_axon2Soma_200Hz.dat';
f_axon2Soma_thin_s{end+1} = 'savedData_VENL_soma_monoSynaptic_ThinDend_axon2Soma_300Hz.dat';
f_axon2Soma_thin_s{end+1} = 'savedData_VENL_soma_monoSynaptic_ThinDend_axon2Soma_500Hz.dat';
f_axon2Soma_thin_s{end+1} = 'savedData_VENL_soma_monoSynaptic_ThinDend_axon2Soma_700Hz.dat';
f_axon2Soma_thin_s{end+1} = 'savedData_VENL_soma_monoSynaptic_ThinDend_axon2Soma_900Hz.dat';
f_axon2Soma_thin_s{end+1} = 'savedData_VENL_soma_monoSynaptic_ThinDend_axon2Soma_1100Hz.dat';

% ----------------------------------------axon2dendrite-----------------------------
f_axon2Dend_thick_s ={};
f_axon2Dend_thick_s{end+1} ='savedData_VENL_soma_mono-synapticTrain_bigDend[3](1)-200Hz.dat';
f_axon2Dend_thick_s{end+1} ='savedData_VENL_soma_mono-synapticTrain_bigDend[3](1)-300Hz.dat';
f_axon2Dend_thick_s{end+1} ='savedData_VENL_soma_mono-synapticTrain_bigDend[3](1)-500Hz.dat';
f_axon2Dend_thick_s{end+1} ='savedData_VENL_soma_mono-synapticTrain_bigDend[3](1)-700Hz.dat';
f_axon2Dend_thick_s{end+1} ='savedData_VENL_soma_mono-synapticTrain_bigDend[3](1)-900Hz.dat';
f_axon2Dend_thick_s{end+1} ='savedData_VENL_soma_mono-synapticTrain_bigDend[3](1)-1100Hz.dat';
f_axon2Dend_thin_s ={};
f_axon2Dend_thin_s{end+1} ='savedData_VENL_soma_mono-synapticTrain_thinDend[3](1)-200Hz.dat';
f_axon2Dend_thin_s{end+1} ='savedData_VENL_soma_mono-synapticTrain_thinDend[3](1)-300Hz.dat';
f_axon2Dend_thin_s{end+1} ='savedData_VENL_soma_mono-synapticTrain_thinDend[3](1)-500Hz.dat';
f_axon2Dend_thin_s{end+1} ='savedData_VENL_soma_mono-synapticTrain_thinDend[3](1)-700Hz.dat';
f_axon2Dend_thin_s{end+1} ='savedData_VENL_soma_mono-synapticTrain_thinDend[3](1)-900Hz.dat';
f_axon2Dend_thin_s{end+1} ='savedData_VENL_soma_mono-synapticTrain_thinDend[3](1)-1100Hz.dat';

% --------------------------analysis for axon2soma -------------------------
freqStr = {'200Hz','300Hz','500Hz','700Hz','900Hz','1100Hz'};
nFile = length(f_axon2Dend_thick_s);
% read and analize f_axon2Soma_thick_s data.
d2 = nan(105001,nFile);
figure('Name','RecSomaV-ThickDend-axon2soma','Position',[900,200,600,1000]);
for i = 1:nFile
    fileName_t = strcat(pathName,f_axon2Soma_thick_s{i});
    fid = fopen( fileName_t,'r');
    d_tem = fscanf(fid,'%f',[1,inf]);
    d_tem = d_tem(3:end)';
    d2(1:length(d_tem),i) = d_tem;  % remove the first 2 numbers which not the trace data
    fclose(fid);

    x = ( 1:size(d2,1) )./100;    % set the x value (time, unit: ms. Because the sample frequency is 100 per ms )
    subplot(nFile,1,i);
    plot(x, d2(:,i), 'k-');
    title(freqStr{i});
    set(gca,'xlim',[0,1050],'ylim',[-100,50])
end
ax_2 = gca();
ax_2.XLabel.String = 'Time (ms)';
ax_2.YLabel.String = 'Voltage (mV)';

% read and analize f_axon2Soma_thin_s data.
d4 = nan(105001,nFile);
figure('Name','RecSomaV-ThinDend-axon2some','Position',[900,200,600,1000]);
for i = 1:nFile
    fileName_t = strcat(pathName,f_axon2Soma_thin_s{i});
    fid = fopen( fileName_t,'r');
    d_tem = fscanf(fid,'%f',[1,inf]);
    d_tem = d_tem(3:end)';
    d4(1:length(d_tem),i) = d_tem;  % remove the first 2 numbers which not the trace data
    fclose(fid);

    x = ( 1:size(d4,1) )./100;    % set the x value (time, unit: ms. Because the sample frequency is 100 per ms )
    subplot(nFile,1,i);
    plot(x,d4(:,i), 'k-');
    title(freqStr{i});
    set(gca,'xlim',[0,1050],'ylim',[-100,50])
end
ax_4 = gca();
ax_4.XLabel.String = 'Time (ms)';
ax_4.YLabel.String = 'Voltage (mV)';


