function sig = redchi2sig(redchi2, dof)
% function sig = redchi2sig(redchi2, dof)
% reject if sig < alpha (0.05, 0.01, etc...);

sig = gammainc((redchi2 * dof)/2, dof/2, 'lower');
