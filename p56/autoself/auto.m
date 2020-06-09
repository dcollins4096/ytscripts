function f = auto(a)
%
% AUTO
%		f = auto(a)
%		returns as f the auto-correlation of image a.

f = fft2(a);
f = f.*conj(f);
f = ifft2(f);
f = fftshift(f);
f = real(f);
f = f/max(max(f));
