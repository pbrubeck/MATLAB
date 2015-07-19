function Q=EulerRotation(alpha, beta, gamma)
%  Returns the rotation matrix described by the Euler angles.
Q=GivensRotation(2,1,gamma)*GivensRotation(2,3,beta)*GivensRotation(2,1,alpha);
end