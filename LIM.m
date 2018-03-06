function [B,G,N] = LIM(X,TAU)
% input: X(m,n), m = is the number of modes, n is the length of time series
% TAU, lag time
% output: B, G

x1 = X(:,1:end-TAU);
x2 = X(:,TAU+1:end);

Nyr = size(x1,2);

TTAU = x2*x1'/Nyr;
T0 = x1*x1'/Nyr;

G = TTAU*inv(T0);
B = real(logm(G))/TAU;
N = -B*T0+T0*B;




