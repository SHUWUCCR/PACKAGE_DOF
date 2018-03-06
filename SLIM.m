function [B,G,N] = SLIM(X,SNS)
% input: X(m*SNS,n), m = is the number of modes, n is the length of time series
% SNS is number of seasons
% TAU, lag time
% output: B, G
M = size(X,1)/SNS;

B = zeros(M,M,SNS-1);
G = zeros(M,M,SNS-1);
N = zeros(M,M,SNS-1);

for isns = 1:SNS-1
    
    x1 = X((isns-1)*M+1:isns*M,:);
    x2 = X(isns*M+1:(isns+1)*M,:);
    
    Nyr = size(x1,2);
    
    TTAU = x2*x1'/Nyr;
    T0 = x1*x1'/Nyr;
    
    iG = TTAU*inv(T0);
    iB = real(logm(iG));
    iN = -iB*T0+T0*iB;
    
    B(:,:,isns) = iB;
    G(:,:,isns) = iG;
    N(:,:,isns) = iN;
    
end



