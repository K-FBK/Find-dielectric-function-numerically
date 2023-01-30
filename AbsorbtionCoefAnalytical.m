function alpha = AbsorbtionCoefAnalytical(eps_real,eps_imag,wl)

ns  =1.375; %Hexane

alpha = (2*pi./wl.*ns).*9*ns^4.*eps_imag./(...
    (eps_real+ 2*ns^2).^2 + eps_imag.^2);

%Equation (1) in (Dement, 2018)

end
