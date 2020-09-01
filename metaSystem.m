function [Ac Bc Cc] = metaSystem (A,B,C,c)
% input is a system xdot = Ax + Bu, y = Cx and an order "c"
% output is the (reduced) meta system of meta-level "c"
% (such a system will be used to construct a Lyapunov function of order 2c)

% get system dimension
n = size(A,1);

% Comput Bc from Wc and Ac
Ac = sparse(A);
Bc = B;
Cc = C;
for i = 1:c-1
    % Compute Ac - sparse matrixes are used for computation speed
    Ac = kron(speye(n), Ac) + kron(A, speye(n^i));
    Bc = kron(B,Bc);
    Cc = kron(C,Cc);   
end
% Reduce dimension
if n == 2
    % Calculate Wc from Proposition 2 (See Paper) 
    W = eye(n);
    for i = 2:c
        W = [[W,  zeros(n^(i-1), 1)]; [zeros(n^(i-1), 1), W]];
    end
    Ac = pinv(W)*Ac*W;
    Bc = pinv(W)*Bc;
    Cc = Cc*W;
else
    Ac = full(Ac);
end

end