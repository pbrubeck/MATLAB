function G=metric(S)
% Computes the metric tensor of a smooth surface
Su=spPartialD(S,1);
Sv=spPartialD(S,2);

xuu=dot(Su, Su, 3);
xuv=dot(Su, Sv, 3);
xvv=dot(Sv, Sv, 3);

G(:,:,1,1)=xuu;
G(:,:,1,2)=xuv;
G(:,:,2,1)=xuv;
G(:,:,2,2)=xvv;
end