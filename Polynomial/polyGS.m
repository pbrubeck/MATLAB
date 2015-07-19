function P=polyGS(P, a, b)
for i=1:size(P,1)
    q=P(i,:);
    for j=1:i-1
        p=P(j,:);
        q=q-polyInt(polyMult(p,q),a,b)*p;
    end
    P(i,:)=q/sqrt(polyInt(polyMult(q,q),a,b));
end
end