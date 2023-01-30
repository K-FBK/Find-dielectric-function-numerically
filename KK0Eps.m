function eps_re = KK0Eps(eps_im,omega)

epsIm1 = eps_im(1:end-1); % k_l
w1 = omega(1:end-1); % w'_l
w2 = omega(2:end); % w'_l+1

w_KK = 0.5*(w1+w2);
L = length(w_KK);

%% 0th Order approximation


eps_re = zeros(L,1);

for i = 1:L
w = w_KK(i);

expression_vector1 = (w2.^2-w^2)./(w1.^2-w^2);
vector_in_sum = epsIm1.*log(abs(expression_vector1));
sum_expression = sum(vector_in_sum);
eps_re(i) =  (1/pi)*sum_expression;
end

[xData, yData] = prepareCurveData(w_KK, eps_re);
% Set up fittype and options.
ft = 'linearinterp';
% Fit model to data.
[fitresult_eps_re, ~] = fit( xData, yData, ft, 'Normalize', 'on' );

eps_re = fitresult_eps_re(omega);
