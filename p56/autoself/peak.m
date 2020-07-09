function [x, y, v] = peak(inp,threshold,blocksize,max_peaks)
% PEAK
%		count peaks in an image
%
%	[x, y, v] = peak(inp, threshold, blocksize,max_peaks)
%
%		inp			input image
%		threshold	threshold value for peaks
%		blocksize	square region of width 2*blocksize+1 around peak
%		max_peaks	maximum number of peaks to seek
%
%		x	row index of peak
%		y	col index of peak
%		v  value at peak

%	Written by: John Loomis

s = size(inp);

if (nargin<4)
   max_peaks = 20;
end
if (nargin<3)
   blocksize = 5;
end
if (nargin<2)
   threshold = 0.05;
end

np = 0;
disp(sprintf('image size: %d x %d',s(1),s(2)))
disp(sprintf('block size: %d', blocksize))
disp(sprintf('threshold: %g',threshold));
disp(sprintf('max # peaks: %d\n',max_peaks));


while np < max_peaks
   
	[col_val, row_idx] = max(inp);
	[val, col] = max(col_val);
   row = row_idx(col);
   if val<threshold
      break;
   end
   
	% erase region of peak
	ra = max(row-blocksize,1);
	rb = min(row+blocksize,s(1));
	ca = max(col-blocksize,1);
	cb = min(col+blocksize,s(2));
   
	inp(ra:rb,ca:cb) = 0;
   
   % save results
	np = np + 1;
	x(np) = row;
	y(np) = col;
	v(np) = val;
   

	if (np==1)
      disp(sprintf('%4s\t%8s\t%4s\t%4s','peak','value','row','col'))
   end
   disp(sprintf('%4d\t%8.6f\t%4d\t%4d',np,val,row,col))
end

if (np==0)
   disp('no peaks found exceeding threshold')
end

