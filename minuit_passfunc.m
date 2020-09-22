function retval = minuit_passfunc(p, data)
% this function is a fake function that you can pass to fminuit for
% minimization.
%
% data structure must contain following
%     data.x: the datapoints
%     data.y_obs: the observed function values
%     data.y_err: the observed error of the Y values
%     data.func: A function handle for the fit model.
%                func(p, x) must return something samesize as y
%     data.minfunc: A function handle for minimization.
%                   minfunc(y_obs, y_err, y_model) needs to return a scalar
%                   Usually, this is chi squared.
%     


y_model = data.func(p, data.x);
retval = data.minfunc(data.y_obs, data.y_err, y_model);