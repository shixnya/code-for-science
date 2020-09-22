function h = plotENEN(xen, yen, varargin)
% let's do it for now just for one line

h = errorbar(xen.value, yen.value, yen.err, yen.err, xen.err, xen.err, varargin{1:end});

