
% Plot summarized responses by group - where summary can be mean/median of each group. Based on plotAllResponsesByGroup.m
% Specifically tailored for Stacey and Erik's obese mouse dataset to compare pre and post boost vaccination responses in 
% the Obese and WT mice in the different vaccination groups.
% Dec. 28 2014
% 

function [figInd] = plotPreVsPostSummaryResponsesByGroupObese(mPeptResponses,expGroupStruct,method,figInd,yLims,fontSize)

if(~exist('yLims','var'))
  yLims = [-500 60000];
end

if(~exist('fontSize','var'))
  fontSize = 12;
end

if(~exist('method','var'))
  method = 'mean';
end

for i=1:length(expGroupStruct)

  figure(figInd);
  figInd = figInd + 1;

  preIndWT  = strmatch('WT_pre',expGroupStruct(i).groupNames);
  postIndWT = strmatch('WT_post',expGroupStruct(i).groupNames);

  preIndOb  = strmatch('Ob_pre',expGroupStruct(i).groupNames);
  postIndOb = strmatch('Ob_post',expGroupStruct(i).groupNames);

  preIndsWT  = expGroupStruct(i).groupInds{preIndWT};
  postIndsWT = expGroupStruct(i).groupInds{postIndWT};
  
  preIndsOb  = expGroupStruct(i).groupInds{preIndOb};
  postIndsOb = expGroupStruct(i).groupInds{postIndOb};
  
  disp(sprintf('found %d samples from group %s',length(preIndsWT),expGroupStruct(i).groupNames{preIndWT}));
  disp(sprintf('found %d samples from group %s\n',length(postIndsWT),expGroupStruct(i).groupNames{postIndWT}));
  
  disp(sprintf('found %d samples from group %s',length(preIndsOb),expGroupStruct(i).groupNames{preIndOb}));
  disp(sprintf('found %d samples from group %s\n\n',length(postIndsOb),expGroupStruct(i).groupNames{postIndOb}));

  switch(method)
    case 'mean'
      subplot(2,1,1)
      hold on;
      if(length(preIndsWT) == 1)
        plot(mPeptResponses(preIndsWT,:),'r'); 
      else
        plot(mean(mPeptResponses(preIndsWT,:)),'r'); 
      end
      if(length(postIndsWT) == 1)
        plot(mPeptResponses(postIndsWT,:),'k');
      else
        plot(mean(mPeptResponses(postIndsWT,:)),'k');
      end
      subplot(2,1,2)
      hold on;
      if(length(preIndsOb) == 1)
        plot(mPeptResponses(preIndsOb,:),'r');
      else
        plot(mean(mPeptResponses(preIndsOb,:)),'r');
      end
      if(length(postIndsOb) == 1)
        plot(mPeptResponses(postIndsOb,:),'k');
      else
        plot(mean(mPeptResponses(postIndsOb,:)),'k');
      end
    case 'median'

      subplot(2,1,1)
      hold on;
      if(length(preIndsWT) == 1)
        plot(mPeptResponses(preIndsWT,:),'r'); 
      else
        plot(median(mPeptResponses(preIndsWT,:)),'r'); 
      end
      if(length(postIndsWT) == 1)
        plot(mPeptResponses(postIndsWT,:),'k');
      else
        plot(median(mPeptResponses(postIndsWT,:)),'k');
      end
      subplot(2,1,2)
      hold on;
      if(length(preIndsOb) == 1)
        plot(mPeptResponses(preIndsOb,:),'r');
      else
        plot(median(mPeptResponses(preIndsOb,:)),'r');
      end
      if(length(postIndsOb) == 1)
        plot(mPeptResponses(postIndsOb,:),'k');
      else
        plot(median(mPeptResponses(postIndsOb,:)),'k');
      end
  end

  subplot(2,1,1)
  title([sprintf('%s WT pre (n = %d) vs. post (n=%d) responses ',expGroupStruct(i).adjuvant,length(preIndsWT),length(postIndsWT))]);
  xlabel('Antigen #');
  ylabel('MFI');

  a = gca;
  set(a,'FontSize',fontSize);
  set(a,'ylim',yLims);
  set(a,'yTick',[10:16])
  set(a,'yTickLabel',[10:16]);
  legend('Pre','Post')
  grid on;
     

  subplot(2,1,2)
  title([sprintf('%s Obese pre (n = %d) vs. post (n=%d) responses ',expGroupStruct(i).adjuvant,length(preIndsOb),length(postIndsOb))]);
  xlabel('Antigen #');
  ylabel('MFI');

  a = gca;
  set(a,'FontSize',fontSize);
  set(a,'ylim',yLims);
  set(a,'yTick',[10:16])
  set(a,'yTickLabel',[10:16]);
  legend('Pre','Post')
  grid on;

  set(gcf,'Position',[100,200,1800,1000]);

end

