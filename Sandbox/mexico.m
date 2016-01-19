A=imread('mexico.png');
r=A(:,:,1);
g=A(:,:,2);
b=A(:,:,3);
idx=(r==247 & g==0 & b==0);
idx=idx|(r==207 & g==10 & b==0);
idx=idx|(r==153 & g==0 & b==0);
idx=idx|(r==115 & g==0 & b==0);
idx=idx|(r==71 & g==0 & b==0);
idx(300:end,:)=0;
idx=conv2(single(idx), ones(2));
idx=uint8(idx(1:end-1,1:end-1));
B=bsxfun(@times, A, idx);
imshowpair(A, B, 'montage')

count=sum(idx(:));
q=pi/18*(cos(pi/9)-cos(pi/6))*(6367444)^2;
area=q*count/14161;
display(area);
W=(7E3/24)*0.7*0.35*0.1*area;
display(W);

n=W/5E10;
display(n);
