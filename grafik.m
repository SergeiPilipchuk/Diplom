ttz1=conv(He1,sumP(conv(Has,[1 1.20 0.167]),conv(As,[-0.373 -4.12 0])));
ttz2=conv(conv(He2,Ds),[1 1.20 0.167]);

t_netf1=conv(He1,conv(Has,Nas));
t_netf2=conv(He2,conv(Ds,Tas));

figure(3);
cla
hold on  

for w=0.2:0.01:0.6
    dl=abs(polyval(ttz1,w*i)/polyval(ttz2,w*i));
    dl_netf=abs(polyval(t_netf1,w*i)/polyval(t_netf2,w*i));
  
    plot(w,dl,'r*',w,dl_netf,'g*')
end
legend('|He*P1*Sa|','|He*(P1+P2*F)|');

