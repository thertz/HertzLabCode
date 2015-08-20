
%parse a NETMHC3.0 result file from standalone version.


% HLA-A0101 : Estimated prediction accuracy  0.811 (using nearest neighbor HLA-A0101)

% -----------------------------------------------------------------------------------
%   pos        HLA    peptide        Identity 1-log50k(aff) Affinity(nM) Bind Level
% -----------------------------------------------------------------------------------
%     0 HLA-A*0101  MKAILVVML gi|227831776|gb         0.028     37082.90

function [alleleList,bindingEnergies] = parseNETMHCpanResultFile(resFile)

fd = fopen(resFile);
i = 1;
while(~feof(fd))

  [C] = textscan(fd,'%d %s %s %s %f %f %s %s','commentStyle','#');
  
  while(isempty(C{1}) && ~feof(fd))
    line = fgetl(fd);
    [C] = textscan(fd,'%d %s %s %s %f %f %s %s','commentStyle','#');
  end

  if(isempty(C{1}))
    break;
  end
    
  alleleList{i} = strrep(C{2}{1},'*','_');
  alleleList{i} = strrep(alleleList{i},'HLA-','');
  bindingEnergies(:,i) = C{5};
    
  i=i+1;
end

fclose(fd);