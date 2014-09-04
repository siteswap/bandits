% import containers.Map

% TODO
% Interop with R and how to run experiment:
% - hash values
% - reduce R range for greater accuracy
function gi = GIBoth(d,L)

l=(0:L)/L;
% R for final set of states
rstate = genStates(d);
n=nstates(d);
Vd= max(repmat(l,n,1), repmat(rstate(:,5),1,L+1));

for r=flip(1:(d-1))
    % Generate transition matrix
    P=trans(r,rstate);
    Vd=P*Vd;
    rstate = genStates(r);
    Rr = rstate(:,5);
    n=nstates(r);
    fixedArm=(d-r+1)*repmat(l,n,1);
    unknownArm=repmat(Rr,1,L+1) + Vd;
    Vd=max(fixedArm, unknownArm);
end

% Note, assumes l is ordered small to large.
[C,I] = min(abs(Vd - l*d));
gi = l(I);

end

function r = R(r,rstate)
    r = rstate(:,3).*rstate(:,4)
end

% full(P) to see full form of sparse matrix
function P = trans(r,rstate)

s=nstates(r);
s1=nstates(r+1);

% P=zeros(s,s1);
% i=1:s;
% % No clicks keeps same index
% P(i,i)= 1-rstate(i,3); %(r-i)/r;
% % Click but no acquisition advances index by c+1
% P(i,i+rstate(i,1)+1)= rstate(i,3)*(1-rstate(i,4));
% % Click and acquisition advances index by c+1+1
% P(i,i+rstate(i,1)+2)=rstate(i,3)*rstate(i,4);


i=1:s;
P=sparse( horzcat(i,i,i), ...
          horzcat(i,i+rstate(i,1)'+1,i+rstate(i,1)'+2), ...
          horzcat(  1-rstate(i,3)', ...
                    rstate(i,3)'.*(1-rstate(i,4)'), ...
                    rstate(i,3)'.*rstate(i,4)'), ...
          s,s1);

end

function n = nstates(r)
    n=(r+1)*r/2;
end

function rstate = genStates(r)

alphap=1;
betap=1;
alphaq=1;
betaq=1;

% List of states for given round
% Ordered first by increasing clicks, 
% then by increasing acquisitions
% Nicely, index of state is not dependent on current round.
rstate=[];
for c=0:(r-1)
    rstate=vertcat(rstate,[repmat(c,1,c+1);0:c]');
end

% col 3 is E[p] TODO - name cols
n=nstates(r);
rstate(:,3)=(rstate(:,1)+alphap)./(alphap+betap+r-ones(n,1));

% col 4 is E[q]
rstate(:,4)=(rstate(:,2)+alphaq)./(alphaq+betaq+rstate(:,1));

% R value of state
rstate(:,5)=rstate(:,3).*rstate(:,4);

end
