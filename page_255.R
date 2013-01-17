# Computational Methods of data analysis
# Source : page 255

############################ SOURCE ####################
#  clear all; close all; % clear all variables and figures
#  L=20; % define the computational domain [-L/2,L/2]
#  n=128; % define the number of Fourier modes 2^n
#  x2=linspace(-L/2,L/2,n+1); % define the domain discretization
#   x=x2(1:n); % consider only the first n points: periodicity
#   dx=x(2)-x(1); % dx value needed for finite difference
#   u=sech(x); % function to take a derivative of
#   ut=fft(u); % FFT the function
#   k=(2*pi/L)*[0:(n/2-1) (-n/2):-1]; % k rescaled to 2pi domain
#   % FFT calculation of derivatives
#   ut1=i*k.*ut; % first derivative
#   ut2=-k.*k.*ut; % second derivative
#   u1=real(ifft(ut1)); u2=real(ifft(ut2)); % inverse transform
#   u1exact=-sech(x).*tanh(x); % analytic first derivative
#   u2exact=sech(x)-2*sech(x).^3; % analytic second derivative
#   % Finite difference calculation of first derivative
#   % 2nd-order accurate
#   ux(1)=(-3*u(1)+4*u(2)-u(3))/(2*dx);
#   for j=2:n-1
#   ux(j)=(u(j+1)-u(j-1))/(2*dx);
#   end
#   ux(n)=(3*u(n)-4*u(n-1)+u(n-2))/(2*dx);
#   % 4th-order accurate
#   ux2(1)=(-3*u(1)+4*u(2)-u(3))/(2*dx);
#   ux2(2)=(-3*u(2)+4*u(3)-u(4))/(2*dx);
#   for j=3:n-2
#   ux2(j)=(-u(j+2)+8*u(j+1)-8*u(j-1)+u(j-2))/(12*dx);
#   end
#   ux2(n-1)=(3*u(n-1)-4*u(n-2)+u(n-3))/(2*dx);

#########################################

L <- 20 # define the computational domain [-L/2, L/2]
n <- 128 #define the number of Fourier modes 2^n

#x2 <- seq(from=-L/2, to = L/2, length.out = n+1) # define the domain discretization
x2 <- seq(from=1, to = L, length.out = n+1)

x <- x2[1:n] # consider only the first n points: periodicity
dx <- x[2]-x[1] #dx value needed for finite difference

u <- acosh(x) # function to take a derivative of
ut <- fft(u) # FFT the function

k <- (2*pi/L)*c(0:(n/2-1),(-n/2):-1) # k rescaled to 2*pi domain

# FFT calculation of derivatives
i <- 1i # in R, i is not defined directly the way of use is 0+1i
ut1 <- i*k *ut   # first derivative
ut2 <- (i*k)^2 * ut # second derivative

u1 <- Re(fft(ut1,inverse=T)) #inversing the transformation
u2 <- Re(fft(ut2,inverse=T)) #inversing the transformation

u1exact=-acosh(x)*tanh(x); #analytic first derivative
u2exact=acosh(x)-2*acosh(x)^3; # analytic second derivative

#plot(x,u1exact,type='l',col='red')
plot(x,u1,col='green')


