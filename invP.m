% по входному P(s) получаем P(-s)
function P = invP(Ps)

[n1,n]=size(Ps); k=-1;
for i=1:n
    k=-k; P(n-i+1)=Ps(n-i+1)*k;
end

