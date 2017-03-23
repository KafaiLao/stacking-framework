function p = softmaxfun(z,T)
p = exp(z/T);
p = p/sum(p);
end