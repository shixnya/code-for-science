function ret = minfunc_chi2(y_obs, y_err, y_model)
% primarily used for minuit chi2 optimization.
if isa(y_model, 'ErrorNum')
    ret = chi2(ErrorNum(y_obs, y_err) - y_model);
    return
end
ret = sum((y_obs - y_model).^2 ./ (y_err).^2, 'all');
