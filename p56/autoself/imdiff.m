function [d,  scl] = imdiff(a,b)

a = im2double(a);
b = im2double(b);
d = a-b;
mx = max(max(d));
mn = min(min(d));
scl = max(mx,mn);
d = (d/scl+1)/2;

if nargout < 2
   disp(sprintf('maximum difference: %g',scl))
end
