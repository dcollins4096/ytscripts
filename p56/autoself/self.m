function f = self(a)
%
% SELF
%		f = self(a)
%		returns as f the self-convolution of image a.

f = fft2(a);
f = f.*f;
f = ifft2(f);
f = fftshift(f);
f = real(f);
f = f/max(max(f));
