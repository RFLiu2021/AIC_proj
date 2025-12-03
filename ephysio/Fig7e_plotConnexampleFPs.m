

close all;
clear all;


cellType = {'VEN','PC','PC','VEN','VEN','PC'};
colors = {[255,56,56]/255; [55,179,74]/255; [0,104,56]/255; [236,41,123]/255; [255,153,204]/255; [27,117,187]/255;}; 

eleFileName= './data_ephy/connExampleFPs_Fig5e.xlsx';
[n,t,r] = xlsread(eleFileName,'Fig5e');

m_FP.fsample = 25000;


figure('position',[1000,700,650,360]);
for i = 1:6
    d = n(:,5*(i-1)+1:5*(i-1)+3);
    t = 10*(1:2500)/m_FP.fsample;
    
    subplot(2,3,i); hold on;
    plot( t, d(:,1),'color',colors{i});
    plot( t, d(:,2),'color',colors{i});
    plot( t, d(:,3),'color',colors{i});
    set(gca,'xlim',[0.01,1],'ylim',[-130,50]);
    if i==1
        plot([0.78,0.98;0.98,0.98],[-40,-40;-40,0],'-k');  % adding the scale bar.
    end
    hold off;
end