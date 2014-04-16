%факторизация поленома
function P_rez=factor(P)

n=length(P)-1;

r=roots(P); k=0; a1=sqrt(abs(P(1)));%n=0;
% for i=1:n
%     if (real(r(1))<=0)
%         n=n+1; 
%     end
% end
for i=1:n
    if (real(r(i))<=0)
        k=k+1; z(k)=r(i);
    end
end
P_rez=a1*poly(z);