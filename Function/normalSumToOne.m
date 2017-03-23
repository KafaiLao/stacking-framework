function p = normalSumToOne(z)
z = z + abs(min(z));
p = z/sum(z);
end