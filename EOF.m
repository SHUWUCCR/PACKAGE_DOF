% EOF decomposition
% [contri L T]=EOF(x)
% L is the eigenvector
% T is the time coefficient
% contri is the contribution of the first M eigenvalues
function [contri L] = EOF(x)
M = 10;
onoff = 0;
n = length(x(1,:));
m = length(x(:,1));
sprintf('n = %d    m = %d',n,m)
if m>n
    S = x'*x/n;
    sub_x = x/sqrt(n);
    onoff = 1;
else
    S = x*x'/n;
end
[L,D] = eig(S,'nobalance');
D = diag(D);
D = reshape(D,1,length(D));
D = fliplr(D);
contri(1:M) = D(1:M)/sum(D);
if 1 == onoff
    L = sub_x*L;
    for i = 1:M
         L(:,n-i+1) = L(:,n-i+1)/sqrt(D(i));
    end
end
L = fliplr(L);
L = L(:,1:M);
disp('ok,your programme of EOF has carried out successfully');



