
% this script read the spike phase height indices from mat files, then plot
% and fitting with linear function. This figure was embeded in the 
% manuscript Fig2 and compared with patch recording data..

% The saved .mat file was copied to folder './ana_AIS2spikePhase' and used
% to compare with the recorded neuron data. Please to read the code 'plot_PhasePeakIdx_vs_AISd.m'



clear all;

nVEN = 1;
peakIdx = nan(7,nVEN);
AIS_D = [0;    50;    70;    90;   110;   130;   170];
for i = 1:nVEN
    fileStr  = strcat('./data/VENs_for_phasePeakIdxFitting/VEN',num2str(i),'/peakHeightIdx.mat');
    load(fileStr);
    peakIdx(:,i) = phaseHeightInd';
end

[idx,sd,se] = imean(peakIdx,2);
figure;hold on;
plot(AIS_D,idx,'o','markerSize',10,'markerFaceColor',[255,56,56]/255,'markerEdgeColor','none');
errorbar(AIS_D,idx,se,'lineStyle','none','Color','k');
set(gca,'xlim',[-20,200],'ylim',[-0.6,0.6],'xTick',0:50:200,'yTick',-0.6:0.3:0.6)
% linear fit the data points.
P=polyfit(AIS_D,idx,1);
plot(AIS_D,polyval(P,AIS_D),'LineStyle','-','Color','k');
box off;



% save('model_PhaseHeightIdx.mat','peakIdx');




