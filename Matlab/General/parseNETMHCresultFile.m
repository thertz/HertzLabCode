
%parse a NETMHC3.0 result file from standalone version.

% Wednesday May 20 2009 14:58


% NetMHC version 3.0. 9mer predictions using Artificial Neural Networks. Allele A1101. 
% Strong binder threshold  50 nM. Weak binder threshold score 500 nM



% ----------------------------------------------------------------------------------------------------
%  pos    peptide      logscore affinity(nM) Bind Level    Protein Name  Allele
% ----------------------------------------------------------------------------------------------------
%    0  MKAILVVML         0.015        42599            gi|227831776|gb|ACP41935|   A1101
%    1  KAILVVMLY         0.676           33         SB gi|227831776|gb|ACP41935|   A1101


%OR OTHER FORMAT....

%NetMHC version 3.0. 9mer predictions using Scoring Matrix. Allele A3303. 
%Strong binder threshold  13.13. Weak binder threshold score   8.67



%----------------------------------------------------------------------------------------------------
% pos    peptide      Score Bind Level    Protein Name  Allele
%----------------------------------------------------------------------------------------------------
%   0  MKAILVVML       -5.8            gi|227831776|gb|ACP41935|   A3303
%   1  KAILVVMLY       -5.6            gi|227831776|gb|ACP41935|   A3303
%   2  AILVVMLYT       -2.4            gi|227831776|gb|ACP41935|   A3303



function [alleleList,bindingEnergies,scoringMatrixIndices] = parseNETMHCresultFile(resFile)

scoringMatrixIndices = [];

fd = fopen(resFile);
i = 1;
while(~feof(fd))
  
  if(i>1)
    sepLine = fgetl(fd);
  end
  
  matrixFormatFlag = 0;
  for j=1:11
    line = fgetl(fd);
    if(strfind(line,'Scoring Matrix'))
      matrixFormatFlag = 1;
      scoringMatrixIndices(i) = 1;
    end
  end
    
  if(matrixFormatFlag) %read scoring matrix NETMHC output
    [C] = textscan(fd,'%d %s %f %s %s %s');

    if(isempty(C{1}))
      break;
    end
    
    %now remove all the SB WB bla bla
    inds1 = strmatch('SB',C{4});
    inds2 = strmatch('WB',C{4});
    inds = [inds1; inds2];
    C{4}(inds) = C{5}(inds);
    C{5}(inds) = C{6}(inds);
    
    alleleList{i} = C{5}{1};
    bindingEnergies(:,i) = C{3};
    
  else
    [C] = textscan(fd,'%d %s %f %f %s %s %s');
    
    if(isempty(C{1}))
      break;
    end
    
    %now remove all the SB WB bla bla
    inds1 = strmatch('SB',C{5});
    inds2 = strmatch('WB',C{5});
    inds = [inds1; inds2];
    C{5}(inds) = C{6}(inds);
    C{6}(inds) = C{7}(inds);
    
    alleleList{i} = C{6}{1};
    bindingEnergies(:,i) = C{4};
  end
  
  i=i+1;
end

fclose(fd);
