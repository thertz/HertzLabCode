function labels=unionfind(edg,n)
% UNIONFIND - implementation of the well-known union-find
%    algorithm. Return labeling file of the connected components,
%    such that a larger component has a smaller label.
%
%    labels = UNIONFIND(edg,n), where edg is the edges list, and n
%    the number of points in the problem.

parent = -1*ones(n,1);
grsize = ones(n,1);
for i=1:size(edg,1)
 a = edg(i,1); b = edg(i,2);
 while(parent(b)~=-1)
   b = parent(b);
 end;
 while(parent(a)~=-1)
   a = parent(a);
 end;
 if (a~=b)
   if grsize(a) > grsize(b)
     c = a; a = b; b = c;
   end;
   grsize(b) = grsize(b)+grsize(a);
   parent(a) = b;
 end;
end;
labels = ones(1,n);
for i=1:n
 a = i;
 while(parent(a)~=-1)
   a = parent(a);
 end;
 labels(i) = a;
end;


[y,x]=hist(labels,1:max(labels));
[mm,mi]=sort(-y);
d(mi) = 1:length(mi);
labels = d(labels);
