function y = ideal(r)
%
% IDEAL
%		y = ideal(r);
%
%		calculates the autocorrelation of the cyl(r) function

r = abs(r);
y = zeros(size(r));
i = find(r<1.0);
y(i) = 2.0*(acos(r(i))- r(i).*sqrt(1.0-r(i).*r(i)))/pi;