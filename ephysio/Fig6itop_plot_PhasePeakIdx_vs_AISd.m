% clear all;
% close all;

load('./data_ephy/NeuronPeakIdx_summary.mat');

bins = 0:10:200;
aveIdx=[];stdIdx=[];seIdx=[];
figure;
plot(axonDist,phaseHeightIdx,'ro')
% plot(axonDist,pahseArea,'ro')
hold on; plot([bins;bins],[ones(1,length(bins))*-0.5;ones(1,length(bins))*1],'k');
set(gca,'xlim',[0,200],'ylim',[-0.2,0.5])

for i = 1:length(bins)-1
     idx_tem = find( axonDist < bins(i+1) & axonDist >= bins(i) );
     if ~isempty(idx_tem)
         [aveIdx(i),stdIdx(i),seIdx(i)] = imean( phaseHeightIdx(idx_tem) );
%          [aveIdx(i),stdIdx(i),seIdx(i)] = imean( pahseArea(idx_tem) );
     else
         aveIdx(i)=nan;
         stdIdx(i)=nan;
         seIdx(i)=nan;
     end
end


figure;
x = bins(1:end-1);
errorbar(x,aveIdx,seIdx,'marker','o','LineStyle','none','markerFaceColor','r','Color','k');
set(gca,'xlim',[-20,200],'ylim',[-0.6,0.6])
% linear fit the data points.
idx = ~isnan(aveIdx);
xx = x(idx); yy = aveIdx(idx);
P=polyfit(xx,yy,1);
hold on;plot(xx,polyval(P,xx),'LineStyle','-','Color','k');
box off;
set(gca,"XTick",[0:50:200],"YTick",[-0.6:0.3:0.6]);
ylabel("Peak Index");
xlabel("Axon to Soma (um)");




