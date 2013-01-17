# Computational Methods of data analysis
# Source : page 253



########################### SOURCE ###################
# clear all; close all; % clear all variables and figures

# L=20; % define the computational domain [-L/2,L/2]
# n=128; % define the number of Fourier modes 2^n

# x2=linspace(-L/2,L/2,n+1); % define the domain discretization
# x=x2(1:n); % consider only the first n points: periodicity

# u=exp(-x.*x); % function to take a derivative of
# ut=fft(u); % FFT the function
# utshift=fftshift(ut); % shift FFT
# figure(1), plot(x,u) % plot initial gaussian
# figure(2), plot(abs(ut)) % plot unshifted transform
# figure(3), plot(abs(utshift)) % plot shifted transform
###########################

L <- 20 # define the computational domain [-L/2, L/2]
n <- 128 #define the number of Fourier modes 2^n

x2 <- seq(from=-L/2, to = L/2, length.out = n+1) # define the domain discretization
x <- x2[1:n] # consider only the first n points: periodicity

u <- exp(-x^2) # function to take a derivative of
ut <- fft(u) # FFT the function

plot(x,u,type='l',xlab='x',ylab='u=exp(-x^2)',main='initial Gaussian') #plot initial gaussian


plot(Re(ut),type='l',xlab='Fourier models',ylab='real[U]') # plotting the Fourier transformation
#but the x axis is the index numbers, not yet marked for frequency

#shiftng the values i.e. spliting the values half and moving them other side
k <- (2*pi/L)*c(0:(n/2-1),(-n/2):-1) # rescaled to 2*pi domain
plot(k,Re(ut), type="l", xlab="Frequencies", ylab="Re(fftshift(u)", main="Fourier modes" ) 

#plotting the absolute values
plot(k,abs(ut), type="l", xlab="Frequencies", ylab="Re(fftshift(u)", main="Fourier modes" )