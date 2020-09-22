function res = MinuitFitEN(func, initp, x, y, varargin)
% function res = MinuitFitEN(func, initp, x, y, varargin)
%
%    func - fit funciton handle func(p, x)
%    initp - (nPar vec) estimated inital parameters
%    x - (any) 2nd argument of func
%    y - (any) ErrorNum data to fit including error
%
%    -- from here optional parameters (default)
%    'Extradim' - (scalar) extra degree of freedom to subtract (0)
%    'LowerLimit' - (nPar vec) lower bound of the parameters (false)
%    'UpperLimit' - (nPar vec) upper bound of the parameters (false)
%    'StepSize' - (nPar vec) step sizes of the parameters (false)
%    'minfunc' - (String) function to minimize ('PoissonNLL')
%    'command' - (String) command to give to fminuit
%                ('set pri -1; scan; minimize; simplex; improve; minos')
%    
%
% Returns:
%    res.dof - (Scalar) degree of freedom for the fit
%    res.p - (nPar, ErrorNum) Estimated parameters with error
%    res.nll2 - (Scalar) return value from the func the name is from
%                  Negative Log Likelihood * 2, which usually it is.
%    res.errmat - (nPar, nPar) Error matrix of the parameter
%    res.redchi2 - (Scalar) reduced chi2 value
%    res.redchi2sig - (Scalar) p-value of the reduced chi2
%    res.rednll2 - (Scalar) nll2 / DOF
%    res.rednll2sig - (Scalar) significance based on nll2 (not sure if
%    valid)
%    res.passdata - (Structure) entire data passed to minuit
%
%
% Description:
%    This function performs a function fit with given function handle and x
%    y data that include Error. It mainly expects the data to be spike
%    count or firing rate (i.e. non negative number with Poisson
%    statistics), but also works with Gaussian distribution data to compute
%    chi2 as a result.
% 
%    if you give 'LowerLimit', also give 'UpperLimit'.
%    if LowerLimit == UpperLimit == 0 for a parameter, that parameter is
%    considered unbounded.
%
% Example:
%    res = MinuitFitEN(@(p, x) p(1) * x + p(2), [0, 0], 1:30,...
%                      ErrorNum(rand(1, 30) + 3, rand(1, 30) + 0.1))
%
% Requires:
%    ErrorNum.m
%    fminuit.m
%    minuit_passfunc.m
%    minfunc_2nllf_poisson.m (for Poisson fit only (default))
%    minfunc_chi2.m (for chi-squre fit only)
%
% Author: Shinya Ito
%         University of California, Santa Cruz (sito@ucsc.edu)
%
% Created: 9/24/2018
% Modified: 6or7/XX/2020 output messages are concealed by evalc.
% Modified: 7/29/2020 description updated


p = inputParser;
addParameter(p, 'Extradim', 0, @isscalar);
addParameter(p, 'LowerLimit', [], @isvector);
addParameter(p, 'UpperLimit', [], @isvector);
addParameter(p, 'StepSize', false, @isvector);
addParameter(p, 'minfunc', 'PoissonNLL', @ischar);
addParameter(p, 'command', 'set pri -1; scan; minimize; simplex; improve; minos', @ischar);
parse(p, varargin{:})

pr = p.Results;

nPar = numel(initp);

res.dof = numel(y) - nPar - pr.Extradim;

if res.dof <= 0
    error('Degree of freedom <= 0');
end

% construct data to pass to fminuit
passdata.func = func;
passdata.x = x;
passdata.y_obs = y.value;
passdata.y_err = y.err;

if strcmp(pr.minfunc, 'PoissonNLL')
    passdata.minfunc = @minfunc_2nllf_poisson;
elseif strcmp(pr.minfunc, 'Chi2')
    passdata.minfunc = @minfunc_chi2;
elseif strcmp(pr.minfunc, 'builtinchi2')
    passdata.minfunc = @chi2;
end

stepbounds = zeros(nPar, 1);
stepbounds(:, 1) = 1:nPar;
% add step size if specified
if pr.StepSize
    stepbounds(:, 2) = pr.StepSize;
end

% add boundaries if specified
if ~isempty(pr.LowerLimit)
    stepbounds(:, end+1) = pr.LowerLimit;
    stepbounds(:, end+1) = pr.UpperLimit;
end

if size(stepbounds, 2) > 1
    %message = '';
    [message, params, params_err, res.nll2, res.errmat, ep, en] = ...
        evalc("fminuit('minuit_passfunc', initp, 'b', '-c', pr.command, '-s', stepbounds, passdata);");
%    [params, params_err, res.nll2, res.errmat, ep, en] = ...
%        fminuit('minuit_passfunc', initp, 'b', '-c', pr.command, '-s', stepbounds, passdata);
else
    [message, params, params_err, res.nll2, res.errmat, ep, en] = ...
        evalc("fminuit('minuit_passfunc', initp, 'b', '-c', pr.command, passdata);");
end

res.p = ErrorNum(params, params_err);

% if x is numeric, it calculates the chi2 value
if isnumeric(x)
    res.redchi2 = sum((func(params, x) - y.value).^2 ./ y.err.^2, 'all') / res.dof;
    res.redchi2sig = redchi2sig(res.redchi2, res.dof);
else
    res.redchi2 = nan;
    res_redchi2sig = nan;
end
res.passdata = passdata;
res.rednll2 = res.nll2 / res.dof;
res.rednll2sig = redchi2sig(res.rednll2, res.dof);
res.params_err_pos = ep; % saving the assymetric errors
res.params_err_neg = en;
res.message = message;


