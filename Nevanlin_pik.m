a11=-0.0454; a12=0.560; a21=0.0267; a22=-0.408;
h1=-0.0648; h2=-0.00456;
b1=-0.0132; b2=-0.00742;
k1_=0.912; k2_=6.11; k3_=-2.22; k0=-0.335; v=3.45;
kR=0.0004;


A=[-0.0454 0.560 0; 0.0267 -0.408 0; 0 1 0];
h=[-0.0648; -0.00456; 0];
c=[0 0 1];
g=[0.0335; 0.00446; 0.0944];
d=0;
E=eye(3);

% построение вспосогательных передаточных функций (в символьном виде)
% He(s)=c((Es-A+gc)^-1)h
tl=ss(A-g*c,h,c,d);
tf(tl)

% A(s)=det(Es-A)
syms s;
A_s=det(E*s-A)

% delta(s)
b=[b1; b2; 0];
k_=[k1_ k2_ k3_];
k=k_+v*c;
D_s=det([E*s-A -b; -k s-k0])

% Ha(s)
Ha_s=det([E*s-A g; -k v])

% теперь эти же передаточный функции в виде полиномов
% He(s)
He1=[-0.00456 -0.001937];
He2=[1 0.5478 0.05083 0.001434];

% A(s)
As=[1 0.4534 0.00357 0];

% Delta(s)
Ds=[1 0.7884 0.2128 0.02323 0.00084785];

% Ha(s)
Has=[3.45 1.7381446 0.0918 0.0017638879];

% теперь задаем весовой частотный множитель Sa(S)=Na(s)/Ta(s)
Nas=[6 0.075 1.5];
Tas=[6.7 7.7 1];

% находим G(s). (факторизация). Gs=G(s) G_s=G(-s)
TT=conv(Tas,invP(Tas));
AA=conv(As,invP(As));
SS=[-1 0 0];
NN=kR*conv(Nas,invP(Nas));
HH=conv(Has,invP(Has));
S1=conv(AA,conv(TT,SS));
S2=conv(NN,HH);
Gs=factor(sumP(S1,S2));
G_s=invP(Gs);

% корни G(-s)
g=roots(G_s);
n=length(G_s)-1;

% находим R_gamma(s). (Факторизация)
GG=1.0455*conv(Gs,G_s);
S3=-kR*conv(HH,TT);
Rg=factor(sumP(GG,S3));

% набор b_i
b_i=-polyval(Tas,g).*(polyval(Tas,-g)).*(polyval(As,-g)).*g; b_i=b_i./(polyval(Nas,g).*polyval(Rg,g));

% находим матрицу L(gamma)

for i=1:n
    for j=1:n
        L(i,j)=(1-b_i(i)*conj(b_i(j)))/(g(i)+conj(g(j)));
    end
end
% проверка существования решения задачи Неванлинны-Пика для наших данных
eig(L);

b=b_i;
% после этого цикла вектор b другой!
for i=2:n
    for j=i:n
        b(j)=((b(j)-b(i-1))/(1-b(j)*conj(b(i-1))))/((g(j)-g(i-1))/(g(j)+conj(g(i-1))));
    end
end

% Последовательное вычисление решения E(s). от E6(s) до E(s). Ei=P1_/P2_;
% Ei*Agi=P1/P2;
% и с учетом, что мы знаем, что в начальных данных есть две пары
% комплексносопряженных чисел
P1_=[1];
P2_=[1];
for i=n:-1:1
    P1=conv(P1_,[1 -g(i)]);
    P2=conv(P2_,[1 conj(g(i))]);
    P1_=sumP(P1,b(i)*P2);
    P2_=sumP(P2,conj(b(i))*P1);
%     if i==5
%         [P1__ P2__]=GR(P1_,P2_);P1_=P1__;P2_=P2__;
%     end
%     if i==3
%         [P1__ P2__]=GR(P1_,P2_);P1_=P1__;P2_=P2__;
%     end
end

% решение задачи Неванлинны-Пика в виде E(s)=m1(s)/m2(s);
m1=P1_;
m2=P2_;


m1R=real(m1);
m1I=imag(m1);
m2R=real(m2);
m2I=imag(m2);

m1_=sumP(conv(m1R,m2R),conv(m1I,m2I));
m2_=sumP(conv(m2R,m2R),conv(m2I,m2I));

% % Строим передаточную функцию искомого фильтра
% F1=conv(Has,sumP(conv(Rg,conv(Nas,m1_)),conv(conv(TT,invP(As)),conv([1 0],m2_))));F1=conv(F1,[1 0]);
% F2=conv(m2_,conv(Gs,G_s));
% % Fs

% F1=poly([-13.8168;-1.0030; -0.9968;-0.4454; -0.4453; -0.1688; -0.1104+0.0512i; -0.1104-0.0512i; -0.0287+0.0192i; -0.0287-0.0192i; -0.0206+0.0168i; -0.0206-0.0168i; 0.0067; -0.0067])
[q r]=deconv(sumP(conv(Rg,conv(Nas,m1_)),conv(conv(TT,invP(As)),conv([1 0],m2_))),G_s);
F1=conv(Has,q);F1=conv(F1,[1 0]);
F2=conv(m2_,Gs);

% редукция до 2-й степени
 tt=tf(F1,F2);
 tt2=balred(tt,2);
