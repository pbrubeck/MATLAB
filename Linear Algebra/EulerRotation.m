function Q=EulerRotation(alpha, beta, gamma)
Q=GivensRotation(2,1,gamma)*GivensRotation(2,3,beta)*GivensRotation(2,1,alpha);
end