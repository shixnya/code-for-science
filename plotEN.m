function h = plotEN(varargin)
% function plotEN(x, EN_y, args for errorbar)
% or
% function plotEN(EN_y, args for errorbar)
%
%    x - (npoint,nseries) x location of the data
%    EN_y - ErrorNum(npoint, nseries) Y data as ErrorNum.
%    for the rest, you can use optional arguments for errorbar.
%    look at the manual of it.
%
% Returns:
%    nothing
%
% Description:
%    Handy function to plot ErrorNum
%
% Example:
%    >> y = ErrorNum([1, 2, 1, 2], [0.2, 0.2, 0.3, 0.3]);
%    >> plotEN(y)
%    >> plotEN([1, 2, 5, 6], y, 'ro')
%
% Requires:
%    ErrorNum.m
%
% Author: Shinya Ito
%         University of California, Santa Cruz (sito@ucsc.edu)
%
% Created: 9/21/2018
% Modified: 9/21/2018, this manual is added

% function plotEN(varargin)
% a function to plot ErrorNum values.

if isa(varargin{1}, 'ErrorNum') % first one is ErorrNum
    % just make up the first argment for this function and call again
    d = varargin{1};
    if isvector(d)
        x = 1:length(d);
    else
        x = repmat((1:length(d))', 1, size(d, 2));
    end
    h = plotEN(x, varargin{:});
elseif isa(varargin{2}, 'ErrorNum')
    h = errorbar(varargin{1}, varargin{2}.value, varargin{2}.err, varargin{3:end});
end
    
    
    
    

