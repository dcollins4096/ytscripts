a = imread('circ.tif');
out = auto(a);
x = 0:31;
y = out(33,33+x);
xi = linspace(0,31,100);
yi = ideal(xi/32);
plot(x,y,'o',xi,yi);
xlabel('radius (pixels)');
ylabel('values');
