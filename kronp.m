function out = kronp(in,c)
% given state "x" and meta level "c", compute meta state xm

xm = in;
for i = 1:c-1
    xm = kron(in,xm); 
end

out = xm;

end