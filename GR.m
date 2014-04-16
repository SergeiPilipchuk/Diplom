
function [m1_ m2_] = GR(m1,m2)

m1R=real(m1);
m1I=imag(m1);
m2R=real(m2);
m2I=imag(m2);

m1_=sumP(conv(m1R,m2R),conv(m1I,m2I));
m2_=sumP(conv(m2R,m2R),conv(m2I,m2I));
