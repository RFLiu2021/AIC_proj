clear all;
close all;

patchSeqFile = '../patchSeqMappingTo10X/data/Ins_patchSeq_new.xlsx'; 
[numInfo, txtInfo, table_all] = xlsread(patchSeqFile);
cellNames = txtInfo(1,3:end);
geneNames = txtInfo(7:end,2);
mappingRates = numInfo(1,3:end);
totalCounts = numInfo(3,3:end);
cellType = txtInfo(2,3:end);
counts = numInfo(4:end,3:end);

% plot the nGenes vs totalCounts.
goodData = counts(:,mappingRates>0.4);
[nGene, nCell] = size(counts);
nDetectGenes = sum(counts~=0,1);
figure('position',[500,600,450,400]);hold on;
scatter(totalCounts(mappingRates>0.4),nDetectGenes(mappingRates>0.4),150,'.','markerEdgeColor','k');
scatter(totalCounts(mappingRates<=0.4),nDetectGenes(mappingRates<=0.4),150,'.','markerEdgeColor','r');
set(gca,'XScale','log','xlim',[0,5e+7],'ylim',[0,12000],'yTick',0:4000:16000);
xlabel('Sequencing depth (total counts)');
ylabel('# of detected genes');

% plot the mappingRates vs totalCounts / nGenes.
figure('position',[500,500,800,300]);
subplot(1,2,1);hold on;
scatter(totalCounts(mappingRates>0.4),mappingRates(mappingRates>0.4),80,'.','markerEdgeColor','k');
scatter(totalCounts(mappingRates<=0.4),mappingRates(mappingRates<=0.4),80,'.','markerEdgeColor','r');
set(gca,'XScale','log','ylim',[0,1],'yTick',0:0.2:1);
xlabel('Sequencing depth (total counts)');
ylabel('Mapping rates');

subplot(1,2,2);hold on;
scatter(nDetectGenes,mappingRates,80,'.','markerEdgeColor','k');
scatter(nDetectGenes(mappingRates>0.4),mappingRates(mappingRates>0.4),80,'.','markerEdgeColor','k');
scatter(nDetectGenes(mappingRates<=0.4),mappingRates(mappingRates<=0.4),80,'.','markerEdgeColor','r');
set(gca,'XScale','log','xlim',[0,50000],'ylim',[0,1],'yTick',0:0.2:1);
xlabel('# of detected genes');
ylabel('Mapping rates');

% plot the Neural and glial expression in the patchseq data set.
counts_snap25 = counts(strcmp(geneNames,'SNAP25'),:);  % 13575: index of gene SNAP25
counts_slc17a7 = counts(strcmp(geneNames,'SLC17A7'),:);  % 13575: index of gene SLC17A7
counts_c1qb = counts(strcmp(geneNames,'C1QB'),:);   % find counts of gene C1QB
counts_gad1 = counts(strcmp(geneNames,'GAD1'),:);  % 13575: index of gene GAD1
randNum = randi(9,[1,length(cellNames)])*0.2-0.8;
markerGenes = {'AQP4'...    % Astrocytes
    ,'C1QA'...     % Microglia
    ,'C1QL1'...    % OPC
    ,'MOG'};       % Oligodendrocyte

figure
for i = 1:length(markerGenes)
    geneIdx = find(strcmp(geneNames,markerGenes{i})==1);
    counts_c1qa = counts(geneIdx,:);  % 13575: index of gene C1QA
    
    subplot(2,2,i);hold on;
    scatter(log2(counts_c1qa(mappingRates>0.4)+1),log2(counts_snap25(mappingRates>0.4)+1)+randNum(1:307),30,'k.','jitter',1);
    scatter(log2(counts_c1qa(mappingRates<=0.4)+1),log2(counts_snap25(mappingRates<=0.4)+1)+randNum(1:45),30,'r.','jitter',1);
    set(gca,'xlim',[-1,17],'ylim',[-1,18],'xTick',0:5:20,'yTick',0:5:20);
    xlabel(strcat(markerGenes{i},', log2(x+1) counts'));
    ylabel('SNAP25, log2(x+1) counts');
    hold off;
end

figure('position',[600,500,1000,400]); 
subplot(1,2,1);hold on;
scatter(log2(counts_c1qb(mappingRates>0.4)+1),log2(counts_snap25(mappingRates>0.4)+1)+randNum(1:307),120,'k.','jitter',1);
scatter(log2(counts_c1qb(mappingRates<=0.4)+1),log2(counts_snap25(mappingRates<=0.4)+1)+randNum(1:45),120,'r.','jitter',1);
set(gca,'xlim',[-1,17],'ylim',[-1,18],'xTick',0:5:20,'yTick',0:5:20);
xlabel('C1QB, log2(x+1) counts');
ylabel('SNAP25, log2(x+1) counts');
subplot(1,2,2);hold on;
scatter(log2(counts_gad1(mappingRates>0.4)+1),log2(counts_slc17a7(mappingRates>0.4)+1)+randNum(1:307),120,'k.','jitter',0.5);
scatter(log2(counts_gad1(mappingRates<=0.4)+1),log2(counts_slc17a7(mappingRates<=0.4)+1)+randNum(1:45),120,'r.','jitter',0.5);
set(gca,'xlim',[-1,17],'ylim',[-1,18],'xTick',0:5:20,'yTick',0:5:20);
xlabel('GAD1, log2(x+1) counts');
ylabel('SLC17A7, log2(x+1) counts');
legend({'PatchSeq, passed QC','PatchSeq, failed QC'});



%% plot the mapping correlation vs total counts vs nGene.
patchSeqmapResFile = '../patchSeqMappingTo10X/data/mappingRes/mappingResults_neurons.csv'; 
[numInfo2, txtInfo2, table_all2] = xlsread(patchSeqmapResFile);
corrs = numInfo2(:,6);
cellNames2 = txtInfo2(2:end,1)';

[cells,ind1,ind2] = intersect(cellNames,cellNames2);
totalCounts2 = totalCounts(ind1);
nGene2 = nDetectGenes(ind1);
corrs2 = corrs(ind2);

figure('position',[500,500,800,300]);
subplot(1,2,1);hold on;
scatter(totalCounts2,corrs2,80,'.','markerEdgeColor','k');
set(gca,'XScale','log','xlim',[0,4e+7],'ylim',[0,1.02],'yTick',0:0.2:2);
xlabel('Sequencing depth (total counts)');
ylabel('Max corr. to 10X clusters');
subplot(1,2,2);hold on;
scatter(nGene2,corrs2,80,'.','markerEdgeColor','k');
set(gca,'xlim',[0,1.5e+4],'ylim',[0,1.02],'yTick',0:0.2:2);
set(gca,'XScale','log');
xlabel('# of detected genes');
ylabel('Max corr. to 10X clusters');






%% testing...





















