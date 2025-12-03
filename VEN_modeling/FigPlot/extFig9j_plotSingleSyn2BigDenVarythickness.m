clear all;
close all;

pathName = './data/VENL_WMsingleSynstim_varyBigBasalDenThick_FigS13j/';
fileName ={};
fileName{end+1} = 'savedData_VENL_bigDen_1.dat';
fileName{end+1} = 'savedData_VENL_bigDen_2.dat';
fileName{end+1} = 'savedData_VENL_bigDen_3.dat';
fileName{end+1} = 'savedData_VENL_bigDen_4.dat';
fileName{end+1} = 'savedData_VENL_bigDen_5.dat';
fileName{end+1} = 'savedData_VENL_bigDen_6.dat';

nFile = length(fileName);
d = nan(25000,nFile);
for i = 1:nFile
    fileName_t = strcat(pathName,fileName{i});
    fid = fopen( fileName_t,'r');
    d_tem = fscanf(fid,'%f',[1,inf]);
    d_tem = d_tem(3:end)';
    d(1:length(d_tem),i) = d_tem;  % remove the first 2 numbers which not the trace data
    fclose(fid);


end


% for i = 1:nFile
%     legendStr{i} =  strcat( num2str(dist(i)),'um');
% end

figure('Position',[1000,800,500,400]);
x = 1:size(d,1);
x = x/100;
% color = [ [207,229,227]/255;[184,199,209]/255;[149,177,212]/255;[123,161,168]/255;[107,155,184]/255;[88,178,220]/255;[0,166,198]/255; [0,89,119]/255; [21,21,211]/255; [0,0,255]/255; ];%;[36,22,84]/255; ];[60,94,145]/255;
color = [[203,27,10]/255;[218,35,55]/255;[255,98,115]/255;[255,71,209]/255;[241,154,190]/255;[247,198,173]/255;[250,218,201]/255];
set(groot,'defaultAxesColorOrder',color)
plot(x,d)%,'color',color{i});%[116,54,50]/255);
ax1 = gca();
ax1.XLim = [60,200];
ax1.YLim = [-71,-64];
ax1.XLabel.String = 'Time (ms)';
ax1.YLabel.String = 'Voltage (mV)';
% legend(legendStr);








