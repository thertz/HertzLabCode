% create a list of peptides for Webby and Ali to test the MF59 induced response.
clear all
close all
 


[arrayPepData] = readH5N1peptideData(); 

inds = [87:89]; 
pepLength = 9;
pepInterval = 3; 

for i=inds  
  currSeq = arrayPepData(i).sequence;
  
  disp(sprintf('Native peptide: %s',currSeq));
  
  str = [ arrayPepData(i).sequence ,'\n'];
  
  pepL = length(currSeq); 
  for begInd = [1:pepInterval:pepL-pepLength];
   currPep =  arrayPepData(i).sequence(begInd:begInd+pepLength);
    str = [str currPep,'\n'];
    for j=1:pepLength
      str = [str currPep(1:j-1),'A',currPep(j+1:end),'\n'];
    end
  end

  sprintf(str)
end
 

for i=inds  
  currSeq = arrayPepData(i).sequence;
  
  disp(sprintf('Native peptide: %s',currSeq));
  
  str = [ arrayPepData(i).sequence ,'\n'];
  
  pepL = length(currSeq); 
  for begInd = [1:pepInterval:pepL-pepLength];
   currPep =  arrayPepData(i).sequence(begInd:begInd+pepLength);
    str = [str currPep,'\n'];
    for j=1:pepLength
      str = [str currPep(1:j-1),'A',currPep(j+1:end),'\n'];
    end
  end

  sprintf(str)
end
 