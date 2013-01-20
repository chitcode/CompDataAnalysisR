# Computational Methods of data analysis
# Source : page 255

########################SOURCE ##########################
# clear all; close all;
# L=30; % time slot to transform
# n=512; % number of Fourier modes 2^9
# t2=linspace(-L,L,n+1); t=t2(1:n); % time discretization
# k=(2*pi/(2*L))*[0:(n/2-1) -n/2:-1]; % frequency components of FFT
# u=sech(t); % ideal signal in the time domain
# figure(1), subplot(3,1,1), plot(t,u,'k'), hold on
# 
# noise=1;
# ut=fft(u);
# utn=ut+noise*(randn(1,n)+i*randn(1,n));
# 
# un=ifft(utn);
# figure(1), subplot(3,1,2), plot(t,abs(un),’k’), hold on

#####################################################

L <- 30
n <- 512
t2 <- seq(from= -L,to=L,length.out=n+1)
t <- t2[1:n]

k <- (2 * pi/(2 * L)) * c(o:n/2-1, -n/2:-1)

sech <- function(a){ 1/cosh(a)}
u <- sech(t)

layout(1:3)
plot(t,u,type='l')

noise <- 1
ut <- fft(u)
utn <- ut + noise *(rnorm(n)+(0+1i)*rnorm(n)) #white noise is always added to frequency

un <- fft(utn,inverse=T)
plot(t,abs(un),type='l')

noise <- 10
utn <- ut + noise *(rnorm(n)+(0+1i)*rnorm(n))
un <- fft(utn,inverse=T)
plot(t,abs(un),type='l')