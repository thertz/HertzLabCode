
function [figInd] = plotSingleAntigenCDFsByGroup(respVec,labels,groupNames,minTreshold,figInd,fontSize)

	figure(figInd);
	hold on;

	uniqLabels = unique(labels);
	numLabels  = length(uniqLabels);
	colors = {'Red','Blue','Black','Green'};
	for i=1:numLabels
		currRespVec = sort(respVec(labels == uniqLabels(i)),'descend');
		h = cdfplot(currRespVec);
		set(h,'Color',colors{i});
		currL(i) = length(currRespVec);
		meanResp(i) = mean(currRespVec);
	end
	a = gca;
	set(a,'FontSize',fontSize);
	legend(groupNames);
	if(minTreshold > 0)
		line([minTreshold minTreshold],[0 1],'Color','Red','LineStyle','-.');
	end
	%line([meanResp(i) meanResp(i)],[ 0 1])

	figInd = figInd + 1;
	

