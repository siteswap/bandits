
% Assumes (1,1) prior
function gi = GICalib(d,L)
    
% TODO - we only need to range from
% R(1,1) to R(t+1,1) and may be repeat
% with tighter ranges for more accurate values.
l=[(0:L)/L];
% R for final set of states
R= flip((1:d)/(d+1))';
Vd= max(repmat(l,d,1), R*ones(1,L+1));

for r=flip(1:(d-1))
    % Generate transition matrix
    P=trans(r+1);
    % Need to explicitly state matrix is sparse
    % otherwise there is no performance benefit.
    Vd=P*Vd;
    R= flip((1:r)/(r+1))';
    Vd=max((d-r+1)*repmat(l,r,1), R*ones(1,L+1) + Vd);
end

% Note, assumes l is ordered small to large.
[C,I] = min(abs(Vd - l*d));
gi = l(I);

end


function P = trans(r)
    
% P=zeros(r-1,r);

% for i=1:(r-1)
%     j=i;
%     P(i,j)=(r-i)/r;
%     P(i,j+1)=i/r;
% end

i=1:(r-1);
P=sparse( horzcat(i,i), horzcat(i,i+1), horzcat((r-i)/r, i/r), r-1,r);

end




