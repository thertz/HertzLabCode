%2d version of pilug, where we collect stats on pairs on  a grid....
%singletons are considered as a pair with the same input on both columns...
function [table]= pilug2d(data)

temp= unique(data,'rows');
if(size(data,2)==2)
  table(:,1:2)= temp;
else
  table(:,1:2)= temp';
end

for i=1:size(table,1)
    table(i,3)= length(intersect(find(data(:,1)==table(i,1)),find(data(:,2)==table(i,2))));
end
    