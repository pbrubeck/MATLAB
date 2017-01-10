function [H] = mobius(z,w)
A=zeros(2);
B=zeros(2);
A(:,1)=z(2)-z([3;1]);
A(:,2)=-z([1;3]).*A(:,1);
B(:,1)=w(2)-w([3;1]);
B(:,2)=-w([1;3]).*B(:,1);
H=B\A;
end