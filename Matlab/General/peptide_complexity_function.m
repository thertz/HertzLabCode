% Complexity for a peptide 9mer
% I need the Bioinformatics Toolbox to run this
function[pepcomplexity] = peptide_complexity_function(Peptide)
    J= length(Peptide);
    for j=1:J
        AA = aacount(Peptide); 
        aa(1) = AA.A;
        aa(2) = AA.R;
        aa(3) = AA.N;
        aa(4) = AA.D;
        aa(5) = AA.C;
        aa(6) = AA.Q;
        aa(7) = AA.E;
        aa(8) = AA.G;
        aa(9) = AA.H;
        aa(10) = AA.I;
        aa(11) = AA.L;
        aa(12) = AA.K;
        aa(13) = AA.M;
        aa(14) = AA.F;
        aa(15) = AA.P;
        aa(16) = AA.S;
        aa(17) = AA.T;
        aa(18) = AA.W;
        aa(19) = AA.Y;
        aa(20) = AA.V;
        % according to this order A R N D C Q E G H I L K M F P S T W Y V
        multiplica = cumprod(factorial(aa'));
        pepcomplexity = (1/J)*log(factorial(J)/multiplica(length(aa)))/log(J);
%        pepcomplexity = (1/J).*log(factorial(J)./multiplica(20))./log(J);
    end
end




% %The histogram of complexities
% X = pepcomplexity(:);
% Y = find(X); %Y has the non-zero indices of X
% W = X(Y); %W has the values of X with Y indices
% %hist(W); % for a count histogram
% %for a frequency histogram instead of just count histogram
% [xout, nbin] = hist(W);
% xout2 = xout/length(W);
% figure;
% bar(nbin, xout2);


% The histogram actually agrees with the Swissprot histograms with window
% size 10 in the paper
% Wootton et al. in Computers Chem. Vol. 17, NO. 2, pp. 149-163. 1993.
% For a window of size 9 we can set a threshold: If complexity <= .55 the
% sequence is a low complex 9-mer according the paper.