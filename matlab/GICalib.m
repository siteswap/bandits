
% Assumes (1,1) prior
function gi = GICalib(d,L)
    
l=[(0:L)/L];
% R for final set of states
R= flip((1:d)/(d+1))';
Vd= max(repmat(l,d,1), R*ones(1,L+1));

for r=flip(1:(d-1))
    % Generate transition matrix
    
    P=trans(r+1);
    Vd=P*Vd;
    R= flip((1:r)/(r+1))';
    Vd=max((d-r)*repmat(l,r,1), R*ones(1,L+1) + Vd);
end

[C,I] = min(abs(Vd - l*d));
gi = l(I);

end


