function [aa_indexes,x]=get_aa_indexes(aa_string,justWarn);

if nargin < 2
    justWarn=0;
end

aa_initial = ['A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y'];

aa_num=length(aa_initial);

[n1,n2]=size(aa_string);

temp1 = repmat(reshape(repmat(aa_initial,n2,1),1,n2*aa_num),n1,1);
temp2 = repmat(aa_string,1,aa_num);

temp=temp1-temp2;

[a1,a2]=find(temp==0);

a3=mod(a2,n2);a3(find(a3==0))=n2;
a2 = ceil(a2/n2);

num=length(a1);

for i=1:num    
    aa_indexes(a1(i),a3(i))=a2(i);
end

[x,y] = find(aa_indexes==0);

if ~isempty(x)
     warning(['Terminal error: Couldnt find an index for the peptides on indexes : ' num2str(x')]);
     num=length(x);
     wrong=zeros(1,num);
     for i=1:num
        wrong(i)=aa_string(x(i),y(i)); 
     end
     if x > 1
         [x] = sort(x);x=x([1 ; find(x(2:end)-x(1:end-1))+1]);  
     end
     if justWarn
        warning(['The following found strings (on above indexes) do not represent a valid aminoacid : ' char(wrong) ]);
     else
        error(['The following found strings (on above indexes) do not represent a valid aminoacid : ' char(wrong) ]); 
     end
end
