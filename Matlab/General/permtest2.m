function [stats,p,samp1i,samp2i]=permtest2(N,func,samp1,samp2)
%function [stats,p,samp1i,samp2i]=permtest2(N,func,samp1,samp2)
%   Performs a permutation test by permuting members of samp1
%   and samp2, each time computing func(samp1,samp2)
%   returns a 1-sided p-value and optionally the whole dist of func() evals
n1=length(samp1);
n2=length(samp2);

if size(samp1,1) ~= n1
    samp1=samp1';
end
if size(samp2,1) ~= n2
    samp2=samp2';
end

stats=zeros(N,1);

if nargout>2
    samp1i=zeros(n1,N);
    samp2i=zeros(n2,N);
end


allSamples=[samp1;samp2];
parfor i=1:N
    ind=randperm(n1+n2);
    tmp1=ind(1:n1);
    tmp2=ind(n1+1:n1+n2);
    stats(i)=func(allSamples(tmp1),allSamples(tmp2));
    if nargout>2
        samp1i(:,i)=tmp1;
        samp2i(:,i)=tmp2;
    end
end
p=sum(stats>=func(samp1,samp2))/N;
end