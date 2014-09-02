function P = trans(r)
    
P=zeros(r-1,r);

% Given state in round d, prob 
% you came from state in d-1
% for i=1:(r-1)
%     j=i;
%     P(i,j)=(r-i)/(r-1);
%     P(i,j+1) = j/(r-1);
% end

for i=1:(r-1)
    j=i;
    P(i,j)=(r-i)/r;
    P(i,j+1)=i/r;
end

end