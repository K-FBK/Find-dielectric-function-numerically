
clear variables
close all
load StolenData.mat

data = sortrows(StolenData,1);

E_eV = data(:,1);
E_eV = E_eV(E_eV<=4);
alpha_cm = data(:,2);
alpha_cm = alpha_cm(E_eV<=4);

%% Fit: 'untitled fit 1'.
[xData, yData] = prepareCurveData( E_eV, alpha_cm );

% Set up fittype and options.
ft = fittype( 'smoothingspline' );
opts = fitoptions( 'Method', 'SmoothingSpline' );
opts.SmoothingParam = 0.99999979336845;

% Fit model to data.
[fitresult, gof] = fit( xData, yData, ft, opts );

alpha_cm = fitresult(xData) -0.0523594;

hc_um = 1.23984198;
hc_nm = hc_um*1e+3;

wl_nm_stolen = hc_nm./E_eV;
alpha_cm(wl_nm_stolen>617) = 0;
alpha_nm_stolen = alpha_cm*1e-7;
figure()
plot(wl_nm_stolen,alpha_cm,LineWidth=1.5);
set(gca, fontSize= 14)
xlabel('Wave Length [nm]')
ylabel('Absorbtion Coeficient [nm^{-1}]')

save AbsorbanceQDStolen alpha_nm_stolen wl_nm_stolen

[ReEps_bulk, ImEps_bulk, ReEps_inf, n_inf] = PrepareBulkData(wl_nm_stolen);
figure()
hold on
plot(E_eV,ReEps_bulk)
plot(E_eV,ImEps_bulk)

set(gca,"XDir","reverse")