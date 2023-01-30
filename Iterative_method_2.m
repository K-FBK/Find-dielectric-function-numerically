close all
clear variables
load AbsorbanceQDStolen.mat
%% import Bulk data for initial guess


[wl,opticalDensity] = importExperimentData('QD_abs.xlsx');

[ReEps_bulk, ImEps_bulk, ReEps_inf, n_inf] = PrepareBulkData(wl);

% figure()
% plot(wl_bulk,epsRe_bulk, wl_bulk, epsIm_bilk,wl,opticalDensity)

ns  =1.5;

%calculate angular frequency for KK
c_nm = 299792458*1e+9;
omega = 2*pi*c_nm./wl;

%Calculate absorbtion and kappa for bulk
dz = 50; %nm
mu = opticalDensity * log(10)./dz;



%kappa_Bulk = mu.*wl./(2*pi);





%% Bulk
%[ReEps_1,ImEps_1,~]=KK0(kappa_Bulk,wl,n_inf);


ImEps_1 = EpsImFromMu(ReEps_bulk,mu,wl,ns);
ReEps_1 = ReEps_bulk;


%% Initialize plot
f1 = figure();
ax1= axes(f1);
set(ax1,'FontSize',14)
xlabel('Wave Length [nm]')
ylabel('Absorbtion coeficient [nm^{-1}]')
hold on
plot(ax1,wl,mu,'-' ,LineWidth=1.5, DisplayName='Measured absorbtion')
legend show

f2 = figure();
ax2 = axes(f2);
set(ax2,'FontSize',14)
xlabel('Wave Length [nm]')
ylabel('Dielectric function')
hold on
legend show
plot(ax2,wl,ReEps_bulk,'k-',LineWidth=1.5,DisplayName='Re(\epsilon)Bulk')
plot(ax2,wl,ImEps_bulk,'k:',LineWidth=1.5,DisplayName='Im(\epsilon)Bulk')

plot(ax2,wl,ReEps_1,'k+-',LineWidth=1.5,DisplayName='Re(\epsilon):0',MarkerIndices = 2:20:length(wl))
plot(ax2,wl,ImEps_1,'k+:',LineWidth=1.5,DisplayName='Im(\epsilon):0',MarkerIndices = 2:20:length(wl))

plotUpdate = [2 3 5];
showAll = true; %plot every iteration
deletePlots = true; %
wait = true;

%% The itaration


i = 1;
loss = 1;
while loss > 0.001

    ImEps_2 = EpsImFromMu(ReEps_1,mu,wl,ns);
    KK0_eps = KK0Eps(ImEps_2,omega);
    ReEps_inf_2 = ReEps_inf-KK0_eps(1);
    ReEps_2 = ReEps_inf_2 + KK0_eps;
    mu2 = AbsorbtionCoef(ReEps_2,ImEps_2,wl,ns);


    if ismember(i,plotUpdate) | showAll
        [pMu,pRe, pIm] = UpdatePlot(ax1,ax2,wl,mu2,ReEps_2,ImEps_2,i);
        %choose a color
        if wait
            pause(1)
        end

        if deletePlots
            delete(pIm);
            delete(pRe)
            delete(pMu)
        end
    end
    %% Loss function

    loss = sum((mu2-mu).^2)./sum(mu.^2);
    dispStr = strcat('Itaration:  ',num2str(i),'  loss:  ', num2str(loss));
    disp(dispStr);


    %% Send values to next iteration


    ImEps_1 = ImEps_2;
    ReEps_1 = ReEps_2;
    i = i+1;
end




%% Plot final values


%new plot
lineWidth = 1.5;
mIndex = 2:20:length(wl);
muString = strcat('After iteration:', num2str(i));
pMu = plot(ax1,wl,mu2,'r--',LineWidth = lineWidth,DisplayName=muString,MarkerIndices= mIndex);

reStr = strcat('Re(\epsilon):', num2str(i));
imStr = strcat('Im(\epsilon):', num2str(i));
pIm = plot(ax2,wl,ImEps_2,'r*:',LineWidth = lineWidth,DisplayName=imStr,MarkerIndices= mIndex);
pRe = plot(ax2,wl,ReEps_2,'r*-',LineWidth = lineWidth,DisplayName=reStr,MarkerIndices= mIndex);


%% Create fit for the dielectric function of omega
[xData, yData] = prepareCurveData(omega, ReEps_2);
% Set up fittype and options.
ft = 'linearinterp';
% Fit model to data.
[DielectricQD_real, gof_dielQD_real] = fit( xData, yData, ft, 'Normalize', 'on' );


[xData, yData] = prepareCurveData(omega, ImEps_2);
% Set up fittype and options.
ft = 'linearinterp';
% Fit model to data.
[DielectricQD_imag, gof_dielQD_imag] = fit( xData, yData, ft, 'Normalize', 'on' );

save dielectric_fit_QD DielectricQD_real DielectricQD_imag


function ImEps = EpsImFromMu(ReEps,mu,lambda,ns)

A = mu.*(ReEps + 2*ns^2).^2;
B = 2*pi*9*ns^3./lambda;
ImEps = A./(B - mu);
end

function mu = AbsorbtionCoef(ReEps,ImEps,wl,ns)

A = 2*pi*9*ns^3./wl;
B = (ReEps + 2*ns^2).^2 + ImEps.^2;

mu = (A./B).*ImEps;
end

function [pMu,pRe, pIm] = UpdatePlot(ax1,ax2,wl,mu2,ReEps_2,ImEps_2,i)
plotM = {'>','<','+','x','*','square','^','diamond','pentagram','hexagram','v','o'};
r = 0.9;

theta = (pi/2)*rand(1,1);
phi = (pi/2)*rand(1,1);
plotColor = [r*sin(theta)*cos(phi),r*sin(theta)*sin(phi),r*cos(theta)];
lineWidth = 2;

marker = plotM(randi([1,length(plotM)],1));
mIndex = 2:20:length(wl);
muString = strcat('After iteration:', num2str(i));
pMu = plot(ax1,wl,mu2,'--',LineWidth = lineWidth,DisplayName=muString,Marker=marker,Color=plotColor,MarkerIndices= mIndex);


reStr = strcat('Re(\epsilon):', num2str(i));
imStr = strcat('Im(\epsilon):', num2str(i));
pIm = plot(ax2,wl,ImEps_2,':',LineWidth = lineWidth,DisplayName=imStr,Marker=marker,Color=plotColor,MarkerIndices= mIndex);
pRe = plot(ax2,wl,ReEps_2,'-',LineWidth = lineWidth,DisplayName=reStr,Marker=marker,Color=plotColor,MarkerIndices= mIndex);
end