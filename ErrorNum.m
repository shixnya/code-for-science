classdef ErrorNum
    % ErrorNum class
    % This class is used for tying error to a number.
    % You can let it do basic arithmetics and also covariance calculation
    % It displays the number with reasonable significant figures.
    % For detailed usage, please see my blog at shixnya.org.
    % Basics:     https://shixnya.org/?p=5
    % Advanced:   https://shixnya.org/?p=82
    %
    % Initialization:
    %    >> en = ErrorNum([3, 2], [0.2, 0.2])
    %
    % returns,
    %
    % en =
    % ErrorNum:
    % 3.0 ± 0.2, 2.0 ± 0.2
    %
    % or,
    %
    %    >> en = ErrorNum.create(1:5, 2); % second argument is dimension
    %
    % returns,
    %
    % en =
    % ErrorNum:
    % 3.0 ± 0.7  % the error is SEM
    %
    % or,
    %
    %    >> en = ErrorNum.init([2, 3]); % argument is size
    %
    % returns,
    %
    % en =
    % ErrorNum:
    % 0 ± 0,  0 ± 0,  0 ± 0
    % 0 ± 0,  0 ± 0,  0 ± 0
    
    properties
        vv
        ee
    end
    methods
        function ret = value(obj)
            % returns the values in matrix format
            s = size(obj);
            ret = zeros(s);
            ret(:) = cat(1, obj(:).vv);
        end
        function ret = err(obj)
            % returns the error in matrix format
            s = size(obj);
            ret = zeros(s);
            ret(:) = cat(1, obj(:).ee);
        end
            
        function obj = ErrorNum(v, e)
            % obj = ErrorNum(v, e)
            % v - (n, m) values
            % e - (n, m) errors for the values
            % returns an n-by-m ErrorNum matrix.
            if nargin == 0
                v = 0;
                e = 0;
            end
            
            if nargin < 2
                e = zeros(size(v));
            end
            
            if numel(v) == 1
                obj.vv = v;
                obj.ee = e;
            else
                s = size(v);
                obj = cell(s);
                for i = 1:numel(v)
                    obj{i} = ErrorNum(v(i), e(i));
                end
                obj = reshape([obj{:}], s);
            end
        end


        function pval = sig(obj)
            %function pval = sig(obj)
            % calculate significance of the value relative to 0
            pval = sigmaToPval(obj.value, obj.err, 0);
        end
        
        function pval = sigpos(obj)
            %function pval = sigpos(obj)
            % determine if the value is significantly positive
            pval = obj.sig / 2; % / 2 because we do one tailed test
            pval(obj.value < 0 ) = 1;
        end
        
        function pval = signeg(obj)
            pval = sigpos(0-obj);
        end
        
        function res = plus(obj, another)
            %function res = plus(obj, another)
            % Addition of two matrices
            % basic arithmetics of ErrorNum are defined by
            % covariance based calculation with 0 covariance.
            res = covadd(obj, another, 0);
        end
        function res = minus(obj, another)
            %function res = minus(obj, another)
            res = covsub(obj, another, 0);
        end
        function res = times(obj, another)
            %function res = times(obj, another)
            res = covmult(obj, another, 0);
        end
        function res = rdivide(obj, another)
            %function res = rdivide(obj, another)
            res = covdiv(obj, another, 0);
        end
        
        % for now, matrix operation is overridden by element-wise
        % calculation
        function res = mtimes(obj, another)
            %function res = mtimes(obj, another)
            res = obj.*another;
        end
        function res = mrdivide(obj, another)
            %function res = mrdivide(obj, another)
            res = obj./another;
        end
        
        % addition with covariance
        function res = covadd(obj, another, cov)
            % function res = covadd(obj, another, cov)
            % calculate addition with a covariance value
            if isnumeric(another)
                another = ErrorNum(another, zeros(size(another)));
            end
            if isnumeric(obj)
                obj = ErrorNum(obj, zeros(size(obj)));
            end
            
            newval = obj.value + another.value;
            newerr = sqrt(obj.err.^2 + another.err.^2 + 2 * cov);
            res = ErrorNum(newval, newerr);
        end
        
        % subtraction with covariance
        function res = covsub(obj, another, cov)
            % function res = covsub(obj, another, cov)
            % calculate subtraction with a covariance value
            if isnumeric(another)
                another = ErrorNum(another, zeros(size(another)));
            end
            if isnumeric(obj)
                obj = ErrorNum(obj, zeros(size(obj)));
            end
            newval = obj.value - another.value;
            newerr = sqrt(obj.err.^2 + another.err.^2 - 2 * cov);
            res = ErrorNum(newval, newerr);
        end
        
        % multiplication with covariance
        function res = covmult(obj, another, cov)
            %function res = covmult(obj, another, cov)
            % calculate multiplication with a covariance value
            if isnumeric(another)
                %another = ErrorNum(another, zeros(size(another)));
                res = ErrorNum(obj.value .* another, abs(obj.err) .* another);
                return
            end
            if isnumeric(obj)
                res = ErrorNum(another.value * obj, abs(another.err) .* obj);
                %obj = ErrorNum(obj, zeros(size(obj)));
                return
            end
            newval = obj.value .* another.value;
            newerr = abs(newval).*...
                sqrt((obj.err./obj.value).^2 + (another.err./another.value).^2 +...
                2 * cov ./ newval);
            res = ErrorNum(newval, newerr);
        end
        
        % division with covariance
        function res = covdiv(obj, another, cov)
            % function res = covdiv(obj, another, cov)
            % calculate division with a covariance value
            if isnumeric(another)
                %another = ErrorNum(another, zeros(size(another)));
                res = ErrorNum(obj.value ./ another, abs(obj.err ./ another));
                return
            end
            if isnumeric(obj)
                obj = ErrorNum(obj, zeros(size(obj)));
            end
            newval = obj.value ./ another.value;
            newerr = abs(newval).*...
                sqrt((obj.err./obj.value).^2 + (another.err./another.value).^2 +...
                - 2 * cov ./ (obj.value .* another.value));
            res = ErrorNum(newval, newerr);
        end
        
        % after some subtraction, one can calculate chi2 value of the
        % vector (beware of the significance. dof needs to be properly
        % calculated
        function [res, sig] = chi2(obj, extradof)
            % function [res, sig] = chi2(obj, extradof)
            % evaluate chi2 value of the distribution.
            % (useful when you have something like
            %  (data - model) as an ErrorNum vector
            % extradof is an extra degree of freedom to subtract.
            % (in case you have constraints in your model)
            if nargin < 2
                extradof = 0;
            end
            res = sum((obj.value).^2 ./ (obj.err).^2);
            if isinf(res)
                sig = nan;
            else
                sig = chi2sig(res, length(obj.value) - extradof);
            end
        end
        
        function res = sum(obj, dim)
            % function res = sum(obj, dim)
            % ErrorNum summation.
            if nargin < 2
                dim = 1;
            end
            v = sum(obj.value, dim);
            e = sqrt(sum((obj.err).^2, dim));
            res = ErrorNum(v, e);
        end
        
        function res = mean(obj, dim)
            % function res = mean(obj, dim)
            % ErrorNum average.
            if nargin < 2
                dim = 1;
            end
            leng = size(obj, dim);
            res = sum(obj, dim) / leng;
        end
       
        
        function [res, sig] = flatchi2(obj, extradof)
            % function [res, sig] = flatchi2(obj, extradof)
            % chi2 value compared to the flat distribution.
            % this function subtracts the average, and adds 1 to extradof.
            % compared to plain chi2.
            if nargin < 2
                extradof = 0;
            end
            meanval = mean(obj.value);
            subobj = obj - meanval;
            [res, sig] = subobj.chi2(1 + extradof); % one extra degree of freedom to take out
        end
        
        function res = abs(obj)
            % function res = abs(obj)
            % ErrorNum absolute value.
            % The error doesn't change, so it might be inaccurate
            % if the error is larger than the value (goes negative).
            res = ErrorNum(abs(obj.value), obj.err);
        end
        
        function res = zscore(obj)
            % function res = zscore(obj)
            % simply, values / error
            res = obj.value ./ obj.err;
        end
        
        function res = pow(obj, p)
            % function res = pow(obj, p)
            % enabling power calculation for the ErrorNum.
            res = ErrorNum(obj.value.^p, obj.err .* abs(p));
        end
        
        function res = power(obj, p)
            % function res = power(obj, p)
            % again, power calculation. I forgot why this is necessary.
            res = obj.pow(p);
        end
        
        function res = sqrt(obj)
            % function res = sqrt(obj)
            % another standard function for the ErrorNum class.
            res = obj.pow(1/2);
        end
        
        function string = to_s(obj, pre)
            % function string = to_s(obj, pre)
            % to make ErrorNum to string.
            % This is where it shines.
            % Through the smartdisplay function, it automatically
            % determines the significant figures based on its error,
            % and display them with appropriate precision.
            % When the leading digit of the error is 1, it displays
            % one more digit.
            if nargin < 2
                pre = 0;
            end
            if length(obj) > 1
                m = mean(obj.value);
                s = std(obj.value) / sqrt(length(obj));
                string = smartdisplay(m, s, pre);
            else
                string = smartdisplay(obj.value, obj.err, pre);
            end
        end
        
        function res = fixnanerror(obj, newerr)
            % function res = fixnanerror(obj, newerr)
            % replace nan error with newerr if the value is not nan.
            if nargin < 2
                newerr = 0;
            end
            v = obj.value;
            e = obj.err;
            e(isnan(e) & ~isnan(v)) = newerr;
            res = ErrorNum(v, e);
        end
        
        function res = setminerr(obj, minerr)
            % function res = setminerr(obj, minerr)
            % replace error smaller than minerr with minerr
            % use when too small error is causing problem.
            % (like in Poisson statistics where error could become 0.)
            res = ErrorNum(obj.value, max(obj.err, minerr));
        end
        
        function smartdisp(obj)
            % function smartdisp(obj)
            % print to_s.
            fprintf('%s', obj.to_s);
        end
        
        function plot(obj)
            % function plot(obj)
            % make a simple error number plot
            % For a better plot, use plotEN function.
            errorbar(obj.value, obj.err);
        end
        
        function res = convert(obj, func, dfunc)
            % function res = convert(obj, func, dfunc)
            % by giving function and derivative, you can convert one number
            % to another number.
            % func can be any function that returns a number.
            % dfunc is the 1st order derivative of func.
            % use anonymous function or function handles for them.
            res = ErrorNum(func(obj.value), abs(obj.err .* dfunc(obj.value)));
        end
        
        function [res, resi] = max(obj, dim)
            % pick up the max of the value, and adhere the error of it.
            % this is not the comparison of two arrays.
            % (not the normal 'max' matlab function.)
            % second argument is dimension.
            if nargin < 2
                dim = 1;
            end
            val = obj.value;
            err = obj.err;
            
            [~, linind] = max(val, [], dim, 'linear');
            [maxval, resi] = max(val, [], dim);
            res = ErrorNum(maxval, err(linind));
            
        end
        
        function dispstring = disp(obj, pre)
            % function dispstring = disp(obj, pre)
            % return a display string.
            
            %disp('ErrorNum, Value:')
            %disp(obj.value)
            %disp('Error:')
            %disp(obj.err)
            if nargin < 2
                pre = 0;
            end
            disp('ErrorNum:')
            osize = size(obj);
            if length(osize) > 2
                fprintf('    size: ');
                disp(osize);
                return
            end
            % otherwise display a matrix.
            scell = arrayfun(@(x) x.to_s(pre), obj, 'UniformOutput', false);
            dispstring = '';
            for j = 1:osize(1)
                l = strjoin(scell(j, :), ',  ');
                dispstring = [dispstring l newline];
                disp(l);
            end
        end
        
    end
    methods (Static)
        function en = create(matrix, dim, poissonimpose)
            % function en = ErrorNum.create(matrix, dir, poissonimpose)
            % create an ErrorNum value based on the matrix.
            % it calculates the mean and SEM in the specified dimension
            % (by dim).
            % if poissonimpose is true, it gives a finite value (that
            % corresponds to one spike) when there is no spike.
            % Note that dim is a necessary argument, unlike matlab standard
            % functions.
            % 
            if nargin < 3
                poissonimpose = 0;
            end
            en = ErrorNum(mean(matrix, dim), std(matrix, 0, dim) / sqrt(size(matrix, dim)));
            if poissonimpose % give finite value when there no spike.
                perr = en.err;
                perr(perr == 0) = 1/size(matrix, dim);
                en = ErrorNum(en.value, perr);
            end
        end
        function en = init(varargin)
            % function en = ErrorNum.init(size)
            % Initialize an ErrorNum variable with a specific size.
            % size is actually varargin that is given to zeros().
            % Useful if you want to preallocate the ErrorNum variable.
            en = ErrorNum(zeros(varargin{:}), zeros(varargin{:}));
        end
    end
end