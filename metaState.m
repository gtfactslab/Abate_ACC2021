function xm = metaState (x,c)
% given state "x" and meta level "c", compute meta state xm
  

if c == 1
    xm = x;
	return; 
end

xm = x;
for i = 1:c-1
    xm = kron(x,xm); 
end

end