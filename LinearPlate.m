clear all;
close all;
clc

fs = 44100;
k = 1/fs;
durration = 1;          % synthesised sound lenght in s
dur = fs*durration;
f0 = 100;
H = 0.002;              % plate thickness
rho = 1150;             % nylon https://www.engineeringtoolbox.com/engineering-materials-properties-d_1225.html
E = 3e+9;               % nylon https://www.engineeringtoolbox.com/engineering-materials-properties-d_1225.html
nu = 0.4;               % Poisson�s ratio nu < 0.5
L = 0.2;
r = 1.3;                % grid aspect ratio
mu = 0.25;              % free parameter
loss = [100, 10; 1000, 8]; % loss [freq.(Hz), T60(s), freq.(Hz), T60(s)]
D = E*H^3 / (12 * (1 - nu^2)); % plate flexural rigidity pg.341
kappa = sqrt(D / (rho * H* L^4)); % pg.342 eq.12.3
% kappa = 20;

h = sqrt(kappa*k/mu); % find grid spacing
gamma = 2*f0;         

mu = kappa*k/(h^2);

% set scheme for loss parameters 
% NOTE: the equations are the sama as for a string
% zeta0,1  pg.189 Practical setting for decay times
% sigma0,1 pg.189 eq.7.29
zeta1 = (-gamma^2+sqrt(gamma^4+4*kappa^2*(2*pi*loss(1,1))^2))/(2*kappa^2);
zeta2 = (-gamma^2+sqrt(gamma^4+4*kappa^2*(2*pi*loss(2,1))^2))/(2*kappa^2);
sigma0 = 6*log(10)*(-zeta2/loss(1,2)+zeta1/loss(2,2))/(zeta1-zeta2);
sigma1 = 6*log(10)*(1/loss(1,2)-1/loss(2,2))/(zeta1-zeta2);

Nx = floor(sqrt(r)/h);        % number of x-subdivisions of spatial domain
Ny = floor(1/(sqrt(r)*h));    % number of y-subdivisions of spatial domain
h = sqrt(r)/Nx;

% initialize time-step variables
uNext = zeros(Nx,Ny);
u = zeros(Nx,Ny);
uPrev = zeros(Nx,Ny);

% excite at a current time-step with a hamming_3d at the center
u(ceil(Nx/2-Nx/8):floor(Nx/2+Nx/8),ceil(Ny/2-Ny/8):floor(Ny/2+Ny/8)) = ...
    0.1*hamming_3d(length(ceil(Nx/2-Nx/8):floor(Nx/2+Nx/8)),length(ceil(Ny/2-Ny/8):floor(Ny/2+Ny/8)),1);
figure(1)
mesh(u)
out = zeros(length(1:dur),1);

l = 3:Nx-2;
m = 3:Ny-2;
for n = 1:dur
%     uNext(l,m) = (2-20*mu^2)*u(l,m) + 8*mu^2 * (u(l,m+1) + u(l,m-1) + u(l+1,m) + u(l-1,m)) - ...
%         2*mu^2 * (u(l+1,m+1) + u(l-1,m-1) + u(l+1,m-1) + u(l-1,m+1)) - uPrev(l,m) - ...
%         mu^2 *(u(l,m+2) + u(l,m-2) + u(l+2,m) + u(l-2,m)); % eq. at pg.347 
    uNext(l,m) = (1/(k*sigma0 + 1))*(((-(kappa^2)/h^4)*((u(l+2,m) + u(l-2,m) + u(l,m+2) + u(l,m-2)) + ...
        2*(u(l+1,m+1) + u(l+1,m-1) + u(l-1,m+1) + u(l-1,m-1)) - ...
        8*(u(l+1,m) + u(l-1,m) + u(l,m+1) + u(l,m-1)) + 20*u(l,m)) + ...
        (gamma^2/h^2)*(u(l+1,m) + u(l-1,m) + u(l,m+1) + u(l,m-1) - 4*u(l,m)) + ...
        ((2*sigma1)/(k*h^2))*(u(l+1,m) + u(l-1,m) + u(l,m+1) + u(l,m-1) - 4*u(l,m) - ...
        (uPrev(l+1,m) + uPrev(l-1,m) + uPrev(l,m+1) + uPrev(l,m-1) - 4*uPrev(l,m))))*k^2 + ...
        k*sigma0*uPrev(l,m) + 2*u(l,m) - uPrev(l,m));
        
        
    out(n) = uNext(floor(Nx/2),floor(Ny/2));
%     plot(uNext);
%     ylim([-1,1]);
%     drawnow;
    
    uPrev  = u;
    u = uNext;
end
figure(2)
plot(out)