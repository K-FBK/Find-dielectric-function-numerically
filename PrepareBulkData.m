function  [ReEps_bulk, ImEps_bulk, ReEps_inf, n_inf] =PrepareBulkData(wl_new)

[wl_bulk, n_bulk, k_bulk] = ImportRefractiveIndex('Ninomiya-o_CdSe_bulk.txt');

n_inf = n_bulk(1);

ReEps = n_bulk.^2-k_bulk;
ImEps = 2*n_bulk.*k_bulk;

ReEps_inf = ReEps(1);


%% Fit: 'untitled fit 1'.
[xData, yData_Re] = prepareCurveData( wl_bulk, ReEps );
[~, yData_Im] = prepareCurveData( wl_bulk, ImEps );

% Set up fittype and options.
ft = fittype( 'smoothingspline' );

% Fit model to data.
[Re_fit, Re_gof] = fit( xData, yData_Re, ft );
[Im_fit, Im_gof] = fit( xData, yData_Im, ft );


%Get values on new interval
ReEps_bulk = Re_fit(wl_new);
ImEps_bulk = Im_fit(wl_new);

end