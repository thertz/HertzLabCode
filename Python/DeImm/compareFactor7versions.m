
clear all;
close all;

addpath('../General/');
addpath(genpath('../ThreadingScripts/'));


seqFilename = 'Factor7.fasta';
HLAFilename = 'alleleList.mat';
fastaFlag = 1;


Seqs = fastaread(seqFilename);
[bindingEnergies] = ADTpredict(HLAFilename,seqFilename,fastaFlag);

bindingThreshold = 4;
for i=1:length(bindingEnergies)

    epitopeMap(i,:) = sum(bindingEnergies{i} <= bindingThreshold,2)./size(bindingEnergies{i},2);
end


numMaps = size(epitopeMap,1);

   
subplot(2,1,1)
plot(epitopeMap', 'LineWidth',2);  
        
a = gca;
set(a,'FontSize',14);
title('Factor 7 epitope Maps');
xlabel('Epitope Start Position');
ylabel('Percent HLA binding');
legend({'Factor 7', 'Factor 7a'})


subplot(2,1,2)
plot(epitopeMap(1,:) - epitopeMap(2,:), 'LineWidth',2)

a = gca;
set(a,'FontSize',14);
title('% Difference Map')
xlabel('Epitope Start Position')
inds = find(epitopeMap(1,:) > epitopeMap(2,:))


printFigure(gcf,'Factor7_EpitopeMaps','png') 