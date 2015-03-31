function Q=yawPitchRoll(yaw, pitch, roll)
% Returns a rotation matrix described by aeronautical angles.
Q=GivensRotation(2,1,roll)*GivensRotation(3,1,-pitch)*GivensRotation(3,2,yaw);
end