function p = sigmaToPval(diffvals, sig1, sig2)
%function p = sigmaToPval(diffvals, sig1, sig2)
%
%    diffvals - (scalar) difference of the values
%    sig1 - (scalar) error of the first value
%    sig2 - (scalar) error of the second value
%
% Returns:
%    p - (scalar) statistical significance (p-value)
%
% Description:
%    returns p-value (ANOVA) based on difference and two errors.
%    This is the basis for the ErrorNum significance calculation.
%    Probably it is more intuitive to use ErrorNum instead of using this
%    function directly.
%
% Example:
%    value1 = 3;
%    value2 = 10;
%    error1 = 3;
%    error2 = 2;
%    p = sigmaToPval(value2 - value1, error1, error2)
%    (returns 0.0522 or so).
%
% Requires:
%    none
%
% Author: Shinya Ito
%         University of California, Santa Cruz (sito@ucsc.edu)
%
% Created: ??
% Modified: 2020-05-13



p = erfc(abs(diffvals) ./ sqrt(2*(sig1.^2 + sig2.^2)));