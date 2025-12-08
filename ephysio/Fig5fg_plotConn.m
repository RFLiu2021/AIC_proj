clear all;
close all;

%% % =========  import data from the HEKA file and save to a mat file. =====
% connectfile1 = '../data_electro/20190315/20190315s4c2r1.dat';
% d = HEKA_Importer(connectfile1);
% 
% % Selected the test round 2.
% presyn_C2.d_postsyn_C1 = d.RecTable.dataRaw{2}{1};    % Cell 1 recorded from HEKA channel 1.
% presyn_C2.d_postsyn_C1 = resample(presyn_C2.d_postsyn_C1,1,25);
% presyn_C2.d_presyn_C2 = d.RecTable.dataRaw{2}{2};    % Cell 2 recorded from HEKA channel 2.    
% presyn_C2.d_presyn_C2 = resample(presyn_C2.d_presyn_C2,1,25);
% presyn_C2.d_postsyn_C3 = d.RecTable.dataRaw{2}{3};    % Cell 3 recorded from HEKA channel 3.    
% presyn_C2.d_postsyn_C3 = resample(presyn_C2.d_postsyn_C3,1,25);
% presyn_C2.d_postsyn_C4 = d.RecTable.dataRaw{2}{6};    % Cell 4 recorded from HEKA channel 6.    
% presyn_C2.d_postsyn_C4 = resample(presyn_C2.d_postsyn_C4,1,25);
% presyn_C2.d_postsyn_C5 = d.RecTable.dataRaw{2}{7};    % Cell 5 recorded from HEKA channel 7.    
% presyn_C2.d_postsyn_C5 = resample(presyn_C2.d_postsyn_C5,1,25);
% presyn_C2.d_postsyn_C6 = d.RecTable.dataRaw{2}{8};    % Cell 6 recorded from HEKA channel 8.    
% presyn_C2.d_postsyn_C6 = resample(presyn_C2.d_postsyn_C6,1,25);
% 
% connectfile2 = '../data_electro/20190315/20190315s4c6r1.dat';   % Cell 4 recorded from channel 6.
% d = HEKA_Importer(connectfile2);
% 
% % Selected the test round 2.
% presyn_C4.d_postsyn_C1 = d.RecTable.dataRaw{2}{1}(1:15100,:);    % Cell 1 recorded from HEKA channel 1.    
% presyn_C4.d_postsyn_C1 = resample(presyn_C4.d_postsyn_C1,1,25);
% presyn_C4.d_postsyn_C2 = d.RecTable.dataRaw{2}{2}(1:15100,:);    % Cell 2 recorded from HEKA channel 2.    
% presyn_C4.d_postsyn_C2 = resample(presyn_C4.d_postsyn_C2,1,25);
% presyn_C4.d_postsyn_C3 = d.RecTable.dataRaw{2}{3}(1:15100,:);    % Cell 3 recorded from HEKA channel 3.    
% presyn_C4.d_postsyn_C3 = resample(presyn_C4.d_postsyn_C3,1,25);
% presyn_C4.d_presyn_C4 = d.RecTable.dataRaw{2}{6}(1:15100,:);    % Cell 4 recorded from HEKA channel 6.    
% presyn_C4.d_presyn_C4 = resample(presyn_C4.d_presyn_C4,1,25);
% presyn_C4.d_postsyn_C5 = d.RecTable.dataRaw{2}{7}(1:15100,:);    % Cell 5 recorded from HEKA channel 7.    
% presyn_C4.d_postsyn_C5 = resample(presyn_C4.d_postsyn_C5,1,25);
% presyn_C4.d_postsyn_C6 = d.RecTable.dataRaw{2}{8}(1:15100,:);    % Cell 6 recorded from HEKA channel 8.  
% presyn_C4.d_postsyn_C6 = resample(presyn_C4.d_postsyn_C6,1,25);
% 
% sampFreq = 1000;
% save("connect0315S4.mat","presyn_C2","presyn_C4","sampFreq");



%%  plot the figure
load("./data_ephy/connect0315S4.mat");

nPoint = size(presyn_C2.d_postsyn_C1,1);
t = (1:nPoint)'; % Unit: 1 ms
colors = {[255,56,56]/255; [55,178,74]/255; [0,104,56]/255; [235,41,123]/255; [255,152,203]/255; [27,117,186]/255 };

figure('Position',[500,300,500,800]);
subplot(6,1,1);
plot(t, filter_300Lowpass( mean(presyn_C2.d_presyn_C2, 2) ,sampFreq), 'Color', colors{2},'LineWidth', 1.5 );
set(gca,'XLim',[10,600],'YLim',[-0.1,0.05]);
title('Cell 2: PC');
subplot(6,1,2);
plot(t, filter_300Lowpass( mean(presyn_C2.d_postsyn_C1, 2) ,sampFreq), 'Color', colors{1},'LineWidth', 1.5);
set(gca,'XLim',[10,600],'YLim',[-0.077,-0.072]);
title('Cell 1: VEN');
subplot(6,1,3);
plot(t, filter_300Lowpass( mean(presyn_C2.d_postsyn_C3, 2) ,sampFreq), 'Color', colors{3},'LineWidth', 1.5 );
set(gca,'XLim',[10,600],'YLim',[-0.075,-0.065]);
title('Cell 3: PC');
subplot(6,1,4);
plot(t, filter_300Lowpass( mean(presyn_C2.d_postsyn_C4, 2) ,sampFreq), 'Color', colors{4},'LineWidth', 1.5 );
set(gca,'XLim',[10,600],'YLim',[-0.075,-0.068]);
title('Cell 4: VEN');
subplot(6,1,5);
plot(t, filter_300Lowpass( mean(presyn_C2.d_postsyn_C5, 2) ,sampFreq), 'Color', colors{5},'LineWidth', 1.5 );
set(gca,'XLim',[10,600],'YLim',[-0.075,-0.068]);
title('Cell 5: VEN');
subplot(6,1,6);
plot(t, filter_300Lowpass( mean(presyn_C2.d_postsyn_C6, 2) ,sampFreq), 'Color', colors{6},'LineWidth', 1.5 );
set(gca,'XLim',[10,600],'YLim',[-0.071,-0.07]);
title('Cell 6: PC');
xlabel('Time (ms)');


figure('Position',[500,300,500,800]);
subplot(6,1,1);
plot(t, filter_300Lowpass( mean(presyn_C4.d_presyn_C4, 2) ,sampFreq), 'Color', colors{4},'LineWidth', 1.5 );
set(gca,'XLim',[10,600],'YLim',[-0.1,0.05]);
title('Cell 4: VEN');
subplot(6,1,2);
plot(t, filter_300Lowpass( mean(presyn_C4.d_postsyn_C1, 2) ,sampFreq), 'Color', colors{1},'LineWidth', 1.5);
set(gca,'XLim',[10,600],'YLim',[-0.077,-0.072]);
title('Cell 1: VEN');
subplot(6,1,3);
plot(t, filter_300Lowpass( mean(presyn_C4.d_postsyn_C2, 2) ,sampFreq), 'Color', colors{2},'LineWidth', 1.5 );
set(gca,'XLim',[10,600],'YLim',[-0.075,-0.067]);
title('Cell 2: PC');
subplot(6,1,4);
plot(t, filter_300Lowpass( mean(presyn_C4.d_postsyn_C3, 2) ,sampFreq), 'Color', colors{3},'LineWidth', 1.5 );
set(gca,'XLim',[10,600],'YLim',[-0.072,-0.065]);
title('Cell 4: PC');
subplot(6,1,5);
plot(t, filter_300Lowpass( mean(presyn_C4.d_postsyn_C5, 2) ,sampFreq), 'Color', colors{5},'LineWidth', 1.5 );
set(gca,'XLim',[10,600],'YLim',[-0.073,-0.068]);
title('Cell 5: VEN');
subplot(6,1,6);
plot(t, filter_300Lowpass( mean(presyn_C4.d_postsyn_C6, 2) ,sampFreq), 'Color', colors{6},'LineWidth', 1.5 );
set(gca,'XLim',[10,600],'YLim',[-0.072,-0.067]);
title('Cell 6: PC');







