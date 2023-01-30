clear variables
load CdSe_data.mat


eps_re_data = n.^2-k.^2;
eps_im = 2.*n.*k;

c_nm = 299792458*1e+9; %nm/2
wL = c_nm./(2*pi*omega_data);


eps_re_KK = KK0Eps(eps_im,omega_data);



figure()
hold on
plot(omega_data,eps_re_data,'r-')
plot(omega_data,eps_re_KK,'b-')
