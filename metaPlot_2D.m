function M = metaPlot_2D(P, x0,x1vec, x2vec)
% plot a polynomial (from a 2x2 system) of 2c-th order
% input P psd matrix, x0 initial point
% x1vec and x2vec are a range of points to meshgrid

% for system of order n=2, the number of (unique) homogeneous monomials of a given meta-level will always
% be c + 1

[X1r, X2r] = meshgrid(x1vec, x2vec); 
c = size(P,1) - 1; % meta level

for i = 1:size(X1r, 1)
    for j = 1:size(X1r, 2)

        x = [];
        for g = 0:c
            x(g+1, 1) = X1r(i, j)^(c - g)*X2r(i, j)^(g);
        end

        Zr(i, j) = x.'*P*x;
    end
end
    x0m = [];
    for g = 0:c
        x0m(g+1, 1) = x0(1)^(c - g)*x0(2)^(g);
    end
    l = x0m.'*P*x0m;
    [M, c] = contour(X1r, X2r, Zr, [l, l]);% 'HandleVisibility', 'off');
    delete(c)
    M = M(:, 2:end); % return level set
end