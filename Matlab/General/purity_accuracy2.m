function [pu , ac]=purity_accuracy2(s,t,nc);


freq=zeros(nc,nc);

for i=1:nc
    inds=find(s==i);
    ms(i)=length(inds);
    if isempty(inds)
        freq(i,:)=0;
        ms(i)=1;        % to avoid divide by zero in acurracy calculation.
    else
        for j=1:nc
            freq(i,j)=length(find(t(inds)==j)); % freq(i,j) is the number of (s=i,t=j)
        end
    end
    inds=find(t==i);
    mt(i)=length(inds);
end

inds=find(mt==0);	% get rid of 0/0 : when mt(i)=0 all freq(i,:)=0 too.
mt(inds)=1;	

pu=(freq./( mt'*ones(1,nc))')';
ac=freq./( ones(nc,1)*ms)';

mt(inds)=0;		% return 0 to where it was

%pu=mean(max(pu'));
%ac=mean(max(ac'));
pu=mt/sum(mt)*max(pu')';
ac=ms/sum(ms)*max(ac')';
