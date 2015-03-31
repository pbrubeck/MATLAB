function p=proj(u, v)
% Returns the projection vector of u over v.
p=(v'*v)\(u'*v)*v;
end