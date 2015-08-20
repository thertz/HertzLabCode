function JBpvalue = JarqueBera(Variable)
%
% JBpvalue = JarqueBera(Variable)
%
% returns the p-value of the Jarqu-Bera
% test of Normality.
%
% NOTE: a low p-value (e.g. 0.03) => NOT Normal
%       a high p-value (e.g. 0.1) => Normal
%
%
% Copyright(c): PNath@London.edu 25-Jan-2001
%


NoOfDataPoints = prod(size(Variable));
Variable = reshape(Variable, NoOfDataPoints, 1);
VariableSkewness = skewness(Variable);
VariableKurtosis = kurtosis(Variable);

JB = NoOfDataPoints.*(VariableSkewness.^2/6 + (VariableKurtosis-3).^2/24);
JBpvalue = 1-chis_cdf(abs(JB),2);
