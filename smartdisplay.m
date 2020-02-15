function [resstring, sigstring] = smartdisplay(value, err, precise)
% function [resstring, sigstring] = smartdisplay(value, err, precise)
% a function to make a string of values in a format
% value +/- err, using the leading term of the error.
if nargin < 3
    precise = 0;
end

% this log10(1.95) permit having more digits if the leading order is 1.
digit = floor(log10(err) - log10(1.95));
if precise
    digit = digit - 1;
end



if isinf(digit) % means zero error. special case: follow the num2str rule.
    resstring = sprintf(['%s ' char(177) ' 0'], num2str(value));
elseif digit < 0 % floating number, should be most cases
    format = sprintf('%%.%df', abs(digit));
    resstring = sprintf([format ' ' char(177) ' ' format], value, err);
else
    format = '%.0f';
    % round the value below significant figure
    err2 = round(err .* 10^(-digit)) * 10^(digit);
    resstring = sprintf([format ' ' char(177) ' ' format], value, err2);
end
    
if nargout > 1
    sigval = sigmaToPval(value, err, 0);
    if sigval < 0.0001;
        sigstring = 'p < 0.0001';
    else
        lord = abs(floor(log10(sigval)));
        lord = max(lord, 2); % ensure 2 digits.

        format = sprintf('%%.%df', lord);
        sigstring = sprintf(['p = ' format], sigval);
    end
end