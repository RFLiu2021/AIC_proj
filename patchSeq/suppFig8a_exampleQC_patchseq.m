clear all;
close all;

[num,txt,tbl] = xlsread('./data/B563_S1C08_R1.csv');

figure('position',[600,500,1000,350]);
plot(num(:,2),num(:,3),'k-','lineWidth',1);
text(71,7000,'LM','color','red');
text(162,16000,'UM','color','red');
set(gca,'xlim',[50,200],'ylim',[500,20000],'yTick',[4500,8000,12000,16000,20000],'yTickLabel',{'4';'8';'12';'16';'20'});
set(gca,'xTick',[73,92,110,125,135,141,150,156,164],'xTickLabel',{'20';'100';'200';'300';'400';'500';'1k';'2k';'5k'});
set(gca,'box','off');
xlabel('Length (bp)');
ylabel('Relative fluorescent unit');

