% сумма двух полиномов
function P=sumP(A,B)

[n1,na]=size(A); [n2,nb]=size(B);
n=max(na,nb); As=zeros(1,n); Bs=zeros(1,n);

for i=1:na
    As(n-na+i)=A(i);
end

for i=1:nb
    Bs(n-nb+i)=B(i);
end
P=As+Bs;