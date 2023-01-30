function [eps_real,eps_imag,n_kk]=KK0(kappa,wl_nm,n_h)
%[eps_real,eps_imag]=KK0_eps(kappa,wl_nm,n_h)
%Computs the real and imaginary part of the permitivity
%eps_real + 1i*eps_imag,
%using measured extiction coeficient
%kappa measured at eavelength
%wl_nm
%n_h real refraction index at high frequency
%

%%  Prepare angular frequency
c_nm_s = 299792458*1e+9; %nm/s
w_data = 2*pi*c_nm_s./wl_nm;
%Behin K-K analys


%% the secondary grid



%% organice vectors to simplify expression

k1 = kappa(1:end-1); % k_l
w1 = w_data(1:end-1); % w'_l
w2 = w_data(2:end); % w'_l+1

w_KK = 0.5*(w1+w2);
L = length(w_KK);

%% 0th Order approximation


n_KK = zeros(L,1);

for i = 1:L
w = w_KK(i);

expression_vector1 = (w2.^2-w^2)./(w1.^2-w^2);
vector_in_sum = k1.*log(abs(expression_vector1));
sum_expression = sum(vector_in_sum);
n_KK(i) = n_h + (1/pi)*sum_expression;
end


%% use interpolation to get a smooth function
[xData, yData] = prepareCurveData(w_KK, n_KK );
% Set up fittype and options.
ft = 'linearinterp';
% Fit model to data.
[fitresult_n, ~] = fit( xData, yData, ft, 'Normalize', 'on' );

n_kk = fitresult_n(w_data);


%% Permittivity
eps_real = n_kk.^2 -kappa.^2;
eps_imag  = 2.*n_kk.*kappa;


end