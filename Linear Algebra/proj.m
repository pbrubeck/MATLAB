function p = proj(u, v)
% Returns the projection of u over v.
p=v*((v'*v)\(v'*u));
end