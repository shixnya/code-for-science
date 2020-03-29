function p = sigmaToPval(diffvals, sig1, sig2)
%function p = sigmaToPval(diffvals, sig1, sig2)



p = erfc(abs(diffvals) ./ sqrt(2*(sig1.^2 + sig2.^2)));