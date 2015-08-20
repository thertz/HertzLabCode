
function [figInd] = plotIndividualResponsesForSingleAntigensByGroup(respVec,labels,groupNames,minTreshold,figInd,fontSize)

	figure(figInd);

	binVec = [0:2048:2^16];
	yLim   =  0.6;

	uniqLabels = unique(labels);
	numLabels  = length(uniqLabels);
	for i=1:numLabels
		currRespVec = sort(respVec(labels == uniqLabels(i)),'descend');
		subplot(numLabels,1,i);
		[N,BIN] = histc(currRespVec,binVec)
		bar(binVec,N./sum(N),'histc');
		currL(i) = length(currRespVec);
		meanResp(i) = mean(currRespVec);
	end

	for i=1:numLabels
		subplot(numLabels,1,i);
		hold on;
		a = gca;
		set(a,'FontSize',fontSize);
		title(groupNames{i});
		set(a,'Ylim',[0 yLim]);
		set(a,'Xlim',[0 65000]);
		if(minTreshold > 0)
			line([minTreshold minTreshold],[0 yLim],'Color','Red','LineStyle','-.');
		end
	end
	figInd = figInd + 1;



