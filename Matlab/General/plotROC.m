function [x,y] = plotROC(pos, neg)

epsilon = 1e-250;
% N = 1000; 
% x = linspace(min([pos,neg]), max([pos,neg]), N);
x = sort([pos, neg]);
N = length(x);
for i=1:N, 
    TP(i) = length(find(pos<=x(i)));
    FP(i) = length(find(neg<=x(i)));
    TN(i) = length(find(neg>x(i)));
    FN(i) = length(find(pos>x(i)));
end

x = FP./(FP+TN+epsilon);
y = TP./(TP+FN);
if nargout<2,
    plot(x, y, '+-') ;
    xlabel('FP/(FP+TN)');
    ylabel('TP/(TP+FN)');
end