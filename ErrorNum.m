% definition of a number with an error
% returns the same class with right error
classdef ErrorNum
    % This class is used for tying error to a number.
    % You can let it do basic arithmetics and also covariance calculation
    % It displays the number with reasonable significant figures.
    properties
        vv
        ee
    end
    methods
        
        function ret = value(obj)
            s = size(obj);
            ret = zeros(s);
            ret(:) = cat(1, obj(:).vv);
        end

        function ret = err(obj)
            s = size(obj);
            ret = zeros(s);
            ret(:) = cat(1, obj(:).ee);
        end
            
        
        function obj = ErrorNum(v, e)
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
            pval = sigmaToPval(obj.value, obj.err, 0);
        end
        
        % basic arithmetics of ErrorNum are defined by
        % covariance based calculation with 0 covariance.
        function res = plus(obj, another)
            % assuem that another is also an ErrorNum.
            res = covadd(obj, another, 0);
        end
        function res = minus(obj, another)
            res = covsub(obj, another, 0);
        end
        function res = times(obj, another)
            res = covmult(obj, another, 0);
        end
        function res = rdivide(obj, another)
            res = covdiv(obj, another, 0);
        end
        
        % for now, matrix operation is overridden by element-wise
        % calculation
        function res = mtimes(obj, another)
            res = obj.*another;
        end
        function res = mrdivide(obj, another)
            res = obj./another;
        end
            
        function res = covadd(obj, another, cov)
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
        
        function res = covsub(obj, another, cov)
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
        
        function res = covmult(obj, another, cov)
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
        
        function res = covdiv(obj, another, cov)
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
            if nargin < 2
                dim = 1;
            end
            v = sum(obj.value, dim);
            e = sqrt(sum((obj.err).^2, dim));
            res = ErrorNum(v, e);
        end
        
        function res = mean(obj, dim)
            if nargin < 2
                dim = 1;
            end
            leng = size(obj, dim);
            res = sum(obj, dim) / leng;
        end
       
        
        function [res, sig] = flatchi2(obj, extradof)
            if nargin < 2
                extradof = 0;
            end
            meanval = mean(obj.value);
            subobj = obj - meanval;
            [res, sig] = subobj.chi2(1 + extradof); % one extra degree of freedom to take out
        end
        
        function res = abs(obj)
            res = ErrorNum(abs(obj.value), obj.err);
        end
        
        function res = zscore(obj)
            res = obj.value ./ obj.err;
        end
        
        function res = pow(obj, p)
            res = ErrorNum(obj.value.^p, obj.err .* abs(p));
        end
        
        function res = power(obj, p)
            res = obj.pow(p);
        end
        
        function res = sqrt(obj)
            res = obj.pow(1/2);
        end
        
        function string = to_s(obj, pre)
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
            % replace error smaller than minerr with minerr
            res = ErrorNum(obj.value, max(obj.err, minerr));
        end
        
        function smartdisp(obj)
            fprintf('%s', obj.to_s);
        end
        
        function plot(obj)
            % make a simple error number plot
            errorbar(obj.value, obj.err);
        end
        
        function res = convert(obj, func, dfunc)
            % by giving function and derivative, you can convert one number
            % to another number.
            % func can be any function that returns a number.
            % dfunc is the 1st order derivative of func.
            res = ErrorNum(func(obj.value), abs(obj.err .* dfunc(obj.value)));
        end

        
        function dispstring = disp(obj, pre)
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
        function en = create(matrix, dir, poissonimpose)
            if nargin < 3
                poissonimpose = 0;
            end
            en = ErrorNum(mean(matrix, dir), std(matrix, 0, dir) / sqrt(size(matrix, dir)));
            if poissonimpose % give finite value when there no spike.
                perr = en.err;
                perr(perr == 0) = 1/size(matrix, dir);
                en = ErrorNum(en.value, perr);
            end
        end
        function en = init(varargin)
            en = ErrorNum(zeros(varargin{:}), zeros(varargin{:}));
        end
    end
end