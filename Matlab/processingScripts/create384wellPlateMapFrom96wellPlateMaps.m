
%created 384 well plates from a set of 96 well plates. Assembly options are
%as folows: 
% Column - Standard option in which we load 2 columns of each 96 well plate
% in a single column of the 384 well plate. 
%
% Rows - this optin is with the rainin liquidator, which loads full 96 well plates starting at A1 and then A2 and
% then B1 and then B2.
%
% Manual-Rows - this option is for assembling 384 manually by row and not by column.
% loading 6 wells in each row one after the other to create blocks. This only
% works for partial plates where it is assumed that there are not more than 6
% columns full in each 96 well plate.
% Blocks are created as follows: (note that this only works for partial
% plates, so even with 4 96 well plates you won't have more than 4 blocks!
% 1 
% 2 
% 3 
% 4 
% Manual-Rows-Full - use 6 channel pipette to fill in rows of 384, each with two rows of the 96 well plate
% format. (new July 2013)

function [plateMap] = create384wellPlateMapFrom96wellPlateMaps(plateMaps96wellFormat,assemblyMode)

%first check that we have no more than 4 96 well plates to begin with.
if(length(plateMaps96wellFormat) > 4)
  disp('Too many 96 well plates provided, exiting empty')
  plateMap = [];
  return;
end

plateMap.title = {plateMaps96wellFormat.title};
plateMap.data = cell(16,24);
[plateMap.data{:}] = deal('');
switch(assemblyMode)
  
  case 'Columns'
   
   colIndex = 1;
   for i=1:length(plateMaps96wellFormat)
     
     for j=1:2:size(plateMaps96wellFormat(i).data,2)
       plateMap.data(1:2:16,colIndex) = plateMaps96wellFormat(i).data(:,j);
       plateMap.data(2:2:16,colIndex) = plateMaps96wellFormat(i).data(:,j+1);
       colIndex = colIndex+1;
     end
   end
 
 case 'Rows'
     
  for i=1:length(plateMaps96wellFormat)
    rowIndex = ceil(i/2);
    colIndex = (1-mod(i,2))+1;
    
    for j=1:size(plateMaps96wellFormat(i).data,1)
    
      plateMap.data(rowIndex,colIndex:2:end) = plateMaps96wellFormat(i).data(j,:);
      rowIndex = rowIndex + 2;
    end
  end
    
 case 'Manual-Rows'
  
  rowIndex = 1;
  
  for i=1:length(plateMaps96wellFormat)
    
    for j=1:2:size(plateMaps96wellFormat(i).data,1)
      %note that this never Assembles data into the right 4 blocks on plate....
      plateMap.data(rowIndex,1:2:12) = plateMaps96wellFormat(i).data(j,1:6);
      plateMap.data(rowIndex,2:2:12) = plateMaps96wellFormat(i).data(j+1,1:6);
      rowIndex = rowIndex + 1;
    end
  end
  
  
 case 'Manual-Rows-Full'
  
  rowIndex = 1;
  
  for i=1:length(plateMaps96wellFormat)
    
    for j=1:2:size(plateMaps96wellFormat(i).data,1)
      %note that this never Assembles data into the right 4 blocks on plate....
      plateMap.data(rowIndex,1:2:12) = plateMaps96wellFormat(i).data(j,1:6);
      plateMap.data(rowIndex,2:2:12) = plateMaps96wellFormat(i).data(j,7:12);
      plateMap.data(rowIndex,13:2:24) = plateMaps96wellFormat(i).data(j+1,1:6);
      plateMap.data(rowIndex,14:2:24) = plateMaps96wellFormat(i).data(j+1,7:12);
      rowIndex = rowIndex + 1;
    end
  end
  
  case '12Channel-Full-Row'
      %use this if using a 12 channel pipette to load samples
      %1,13,2,14,3,15....23,12,24
  
  rowIndex = 1;
  
  for i=1:length(plateMaps96wellFormat)
    
    for j=1:2:size(plateMaps96wellFormat(i).data,1)
      %note that this never Assembles data into the right 4 blocks on plate....
      plateMap.data(rowIndex,1:2:24) = plateMaps96wellFormat(i).data(j,1:12);
      plateMap.data(rowIndex,2:2:24) = plateMaps96wellFormat(i).data(j+1,1:12);
    
      rowIndex = rowIndex + 1;
    end
  end
  
end
 
% fd = fopen('plateMap.txt','w');
% for i=1:16
%   for j=1:24
%     fprintf(fd,'%s\t',['#',num2str(M(i,j))]);
%   end
%   fprintf(fd,'\n');
% end
% fclose(fd);