clear all;
close all;
clc

fs = 44100;
k = 1/fs;
durration = 1;          % synthesised sound lenght in s
dur = fs*durration;
H = 0.001;              % plate thickness
rho = 280;              % plate density https://www.wood-database.com/paulownia/
E = 12000;              % ash https://amesweb.info/Materials/Youngs-Modulus-of-Wood.aspx since ash is in the same family as paulownia
nu = 0.3;               % Poisson’s ratio nu < 0.5
L = 1;
r = 1.3;                % grid aspect ratio
mu = 0.25;              % free parameter
D = E*H^3 / (12 * (1 - nu^2)); % plate flexural rigidity
%kappa = sqrt(D / (rho * H* L^4));
kappa = 20;
h = sqrt(kappa*k/mu); % find grid spacing
gamma = 2*200;          %2*f0 tng kurt f0

mu = kappa*k/(h^2);


N = floor(1/h);
h = 1/N;
Nx = floor(sqrt(r)/h);        % number of x-subdivisions of spatial domain
Ny = floor(1/(sqrt(r)*h));    % number of y-subdivisions of spatial domain

uNext = zeros(Nx,Ny);
u = zeros(Nx,Ny);
width = floor(Nx/10);
excitationRange = 1:width;
u(floor(Nx/2),floor(Ny/2)) = 1;
uPrev = u;
figure(1)
plot(u)


l = 3:Nx-2;
m = 3:Ny-2;
for n = 1:dur
    uNext(l,m) = (2-20*mu^2)*u(l,m) + 8*mu^2 * (u(l,m+1) + u(l,m-1) + u(l+1,m) + u(l-1,m)) - ...
        2*mu^2 * (u(l+1,m+1) + u(l-1,m-1) + u(l+1,m-1) + u(l-1,m+1)) - uPrev(l,m) - ...
        mu^2 *(u(l,m+2) + u(l,m-2) + u(l+2,m) + u(l-2,m)); % eq. at pg.347 
    out(n) = uNext(floor(Nx/2),floor(Ny/2));
%     plot(uNext);
%     ylim([-1,1]);
%     drawnow;
    
    uPrev  = u;
    u = uNext;
end
figure(2)
plot(out)