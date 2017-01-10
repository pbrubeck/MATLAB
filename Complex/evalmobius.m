function [w] = evalmobius(H,z)
w=(H(1,1)*z+H(1,2))./(H(2,1)*z+H(2,2));
end