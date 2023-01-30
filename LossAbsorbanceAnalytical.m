function loss = LossAbsorbanceAnalytical(eps,WL_nm,T_experiment)

n1 = 1;
n2 = sqrt(eps);
z = 10; %nm

phi = z*n2*2*pi/WL_nm;

r12 = (n1-n2)./(n1+n2);
%r21 = (n2-n1)/(n1+n2);
t12 = 2*n1./(n1+n2);
t21 = 2*n2./(n1+n2);

r23 = (n2-n1)./(n1+n2);
r21 = (n1-n2)./(n1+n2);
t23 = 2*n2./(n1+n2);
%t32 = 2*n1/(n1+n2);

t13 = t12.*t23.*exp(-1i*phi)./(1-r21.*r23.*exp(-2i*phi));
intencity_trans = abs(t13).^2;

r13 = r12 + t12.*t21.*r23.*exp(-2i*phi)./(1-r21.*r23.*exp(-2i*phi));
I_ref = abs(r13).^2;

T_analytical = 1./intencity_trans

loss = (T_analytical - Absorb_experiment).^2;

end