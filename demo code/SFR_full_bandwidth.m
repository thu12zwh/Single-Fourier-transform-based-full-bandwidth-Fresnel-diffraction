clc
clear all
close all

iflag = -1;
eps = 10^(-12);


z = 10;                                            % propagation distance (mm)
n0 = 1000;                                         % sampling number
pitch = 0.002;                                     % sampling interval
lb = n0*pitch;
xb = linspace(-lb/2,lb/2-lb/n0,n0)';

maxf = 150;                                        % maximum spatial frequency of the input chirp grating (1/mm)
t = chirp(xb,0,max(xb),maxf).*cos(0.4*pi*xb);      % input chirp grating
figure,plot(xb,t)

lam = 500e-6;                                      % illumination wavelength (mm)
k = 2*pi/lam;

%% conventional method
tic
l0 = n0*pitch;
x0 = linspace(-l0/2,l0/2-l0/n0,n0)';
chirp0 = exp(1i*k/2/z*x0.^2);
t_c0 = t.*chirp0;
t_pro0 = fftshift(fft(fftshift((t_c0))));
toc
figure,plot(x0*lam*z/pitch/max(abs(x0))/2,abs(t_pro0)/max(abs(t_pro0)))
title('conventional method amplitude')


%% proposed method
tic
t_FT = fftshift(fft(fftshift(t)));
t_FT_pad = padarray(t_FT,n0/2);
t_pad = ifftshift(ifft(ifftshift(t_FT_pad)));
n = size(t_pad,1);
pitchn = lb/n;
l = n*pitchn;
x = linspace(-l/2,l/2-l/n,n)';
chirp1 = exp(1i*k/2/z*x.^2);
t_c = t_pad.*chirp1;
L = lam*z/pitch;
X = linspace(-L/2,L/2-L/n,n)';
fX = X/lam/z;
K = n/2/(1/2/pitchn);
t_pro2 = nufft1d3(n,x/max(abs(x))*pi,t_c,iflag,eps,n,(fX)*K);
toc
t_pro2 = exp(1i*k/2/z*X.^2).*t_pro2/max(abs(t_pro2));
figure,plot(X,abs(t_pro2));title('proposed method amplitude')


%% analytical integral 
uu = zeros(n,1);
tic
for j = 1:n
    fun = @(xn) 1/2/pi*z./sqrt((X(j)-xn).^2+z^2).*(1./sqrt((X(j)-xn).^2+z^2)...
               -1i*k).*exp(1i*k*sqrt((X(j)-xn).^2+z^2))./sqrt((X(j)-xn).^2+z^2).*chirp(xn,0,max(xb),maxf).*cos(0.4*pi*xn);
           uu(j,1) = integral(fun,min(xb),max(xb));
end
toc
uu = uu/max(abs(uu));
amplitude_rsi = abs(uu);
figure,plot(X,amplitude_rsi);title('Analytical integral amplitude')









