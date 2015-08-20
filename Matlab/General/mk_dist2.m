function [LOC]=mk_dist2(dist,k)
%this program goes back to checking each connection twice (once for each neibor)
%mcc -irv make_f !
%===========================================

N=size(dist,1);

LOC=zeros(k,N);
ONES=ones(1,N);
p=inf;

if k==N
  LOC=1:N^2;LOC(1:N+1:N^2)=2;
else


for i=1:N
        dist(i,i)=p;
        for j=1:k
                val=min(dist(i,:));
                sorted=find(dist(i,:)==val);
                LOC(j,i)=sorted(1);
                dist(i,sorted(1))=p;
        end;
        LOC(:,i)=LOC(:,i)+(i-1)*N;
end;

end





