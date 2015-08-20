function [skewCoeff,skewStat,skewDf,skewPvalue,kurtCoeff,kurtStat,kurtPvalue] = Mskekur(X,alpha)
%Mardia's multivariate skewness and kurtosis.
%Calculates the Mardia's multivariate skewness and kurtosis coefficients
%as well as their corresponding statistical test. For large sample size 
%the multivariate skewness is asymptotically distributed as a Chi-square 
%random variable; here it is corrected for small sample size. Likewise,
%the multivariate kurtosis it is distributed as a unit-normal.
%
%   Syntax: function [Mskekur] = Mskekur(X,alpha) 
%      
%     Inputs:
%          X - multivariate data matrix [Size of matrix must be n(data)-by-p(variables)]. 
%      alpha - significance level (default = 0.05). 
%     Output:
%          Complete statistical analysis table of both multivariate
%          Mardia's skewness and kurtosis.
%
%    Example: For the example of Pope et al. (1980) given by Stevens (1992, p. 249), 
%             with 12 cases (n = 12) and three variables (p = 3). We are interested
%             to calculate and testing its multivariate skewnees and kurtosis with a
%             significance level = 0.05.
%                      --------------    --------------
%                       x1   x2   x3      x1   x2   x3
%                      --------------    --------------
%                      2.4  2.1  2.4     4.5  4.9  5.7
%                      3.5  1.8  3.9     3.9  4.7  4.7
%                      6.7  3.6  5.9     4.0  3.6  2.9
%                      5.3  3.3  6.1     5.7  5.5  6.2
%                      5.2  4.1  6.4     2.4  2.9  3.2
%                      3.2  2.7  4.0     2.7  2.6  4.1
%                      --------------    --------------
%
%          Total data matrix must be:
%           X=[2.4 2.1 2.4;3.5 1.8 3.9;6.7 3.6 5.9;5.3 3.3 6.1;5.2 4.1 6.4;
%           3.2 2.7 4.0;4.5 4.9 5.7;3.9 4.7 4.7;4.0 3.6 2.9;5.7 5.5 6.2;2.4 2.9 3.2;
%           2.7 2.6 4.1];
%
%             Calling on Matlab the function: 
%                Mskekur(X)
%
%             Answer is:
%
% Analysis of the Mardia's multivariate skewness and kurtosis.
% ---------------------------------------------------------------
% Multivariate       Coefficient      Statistic     df       P
% ---------------------------------------------------------------
% Skewness             1.7660           4.9908      10    0.8918
% Kurtosis             9.6443          -1.6936            0.0903
% ---------------------------------------------------------------
% With a given significance level of: 0.05
% The multivariate skewness results not significative.
% The multivariate kurtosis results not significative.
%

%  Created by A. Trujillo-Ortiz and R. Hernandez-Walls
%         Facultad de Ciencias Marinas        
%         Universidad Autonoma de Baja California 
%         Apdo. Postal 453  
%         Ensenada, Baja California
%         Mexico  
%         atrujo@uabc.mx
%
%  May 22, 2003.
%
%  To cite this file, this would be an appropriate format:
%  Trujillo-Ortiz, A. and R. Hernandez-Walls. (2003). Mskekur: Mardia's multivariate skewness
%    and kurtosis coefficients and its hypotheses testing. A MATLAB file. [WWW document]. URL 
%    http://www.mathworks.com/matlabcentral/fileexchange/loadFile.do?objectId=3519&objectType=FILE
%
%  References:
% 
%  Mardia, K. V. (1970), Measures of multivariate skewnees and kurtosis with
%         applications. Biometrika, 57(3):519-530.
%  Mardia, K. V. (1974), Applications of some measures of multivariate skewness
%         and kurtosis for testing normality and robustness studies. Sankhyâ A,
%         36:115-128
%  Stevens, J. (1992), Applied Multivariate Statistics for Social Sciences. 2nd. ed.
%         New-Jersey:Lawrance Erlbaum Associates Publishers. pp. 247-248.
%

if nargin < 2, 
    alpha = 0.05;  %(default)
end; 

if nargin < 1, 
    error('Requires at least one input arguments.'); 
end; 

[n,p] = size(X);

difT = [];
for	j = 1:p
   eval(['difT=[difT,(X(:,j)-mean(X(:,j)))];']);
end;

S = cov(X);  %variance-covariance matrix
D = difT*inv(S)*difT';  %Mahalanobis' distances matrix
b1p = (sum(sum(D.^3)))/n^2;  %multivariate skewness coefficient
b2p=trace(D.^2)/n;  %multivariate kurtosis coefficient

k = ((p+1)*(n+1)*(n+3))/(n*(((n+1)*(p+1))-6));  %small sample correction
v = (p*(p+1)*(p+2))/6;  %degrees of freedom
g1 = (n*b1p*k)/6;  %skewness test statistic:it approximates to a chi-square distribution
P1 = 1 - chi2cdf(g1,v);  %significance value associated to the skewness

g2 = (b2p-(p*(p+2)))/(sqrt((8*p*(p+2))/n));  %kurtosis test statistic:it approximates to
                                             %a unit-normal distribution
P2 = 2*(1-normcdf(abs(g2)));  %significance value associated to the kurtosis

disp('Analysis of the Mardia''s multivariate skewness and kurtosis.')
fprintf('---------------------------------------------------------------\n');
disp('Multivariate       Coefficient      Statistic     df       P')
fprintf('---------------------------------------------------------------\n');
fprintf('Skewness        %11.4f%17.4f%8i%10.4f\n\n',b1p,g1,v,P1);
fprintf('Kurtosis        %11.4f%17.4f%18.4f\n',b2p,g2,P2);
fprintf('---------------------------------------------------------------\n');
fprintf('With a given significance level of: %.2f\n', alpha);
if P1 >= alpha;
   fprintf('The multivariate skewness results not significative.\n');
else 
   fprintf('The multivariate skewness results significative.\n');
end;
if P2 >= alpha;
   fprintf('The multivariate kurtosis results not significative.\n\n');
else 
   fprintf('The multivariate kurtosis results significative.\n\n');
end;

skewCoeff  = b1p;
skewStat   = g1;
skewDf     = v;
skewPvalue = P1;

kurtCoeff  = b2p;
kurtStat   = g2;
kurtPvalue = P2;