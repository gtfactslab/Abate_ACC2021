function [Ac1 Ac2 Bc Cc] = metaSystem (A1, A2,B,C,c)
% input is a system xdot = Ax + Bu, y = Cx and an order "c"
% output is the (reduced) meta system of meta-level "c"
% (such a system will be used to construct a Lyapunov function of order 2c)

% get system dimension
n = size(A1,1);

% Comput Bc from Wc and Ac
Ac1 = sparse(A1);
Ac2 = sparse(A2);
Bc = B;
Cc = C;
for i = 1:c-1
    % Compute Ac - sparse matrixes are used for computation speed
    Ac1 = kron(speye(n), Ac1) + kron(A1, speye(n^i));
    Ac2 = kron(speye(n), Ac2) + kron(A2, speye(n^i));
    Bc = kron(B,Bc);
    Cc = kron(C,Cc);   
end
Ac1 = full(Ac1);
Ac2 = full(Ac2);
end