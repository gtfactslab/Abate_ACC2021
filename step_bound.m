function bound = step_bound(A, b, c)
% given state "x" and meta level "c", compute meta state xm
    a = 'not_solved';
    accuracy = 3;
    while ~isequal(a, 'Solved') && ~(accuracy > 8)
        m = size(A, 1);
        cvx_begin sdp
        variable P(m,m) semidefinite
        minimize matrix_frac(c', P);
        b'*P*b <= 5^accuracy
        A'*P + P*A <= 0
        cvx_end
        
        a = cvx_status;
        accuracy = accuracy + 1;
        pause(1)
    end
        
    bound = (c*inv(P)*c')*(b'*P*b);   
end