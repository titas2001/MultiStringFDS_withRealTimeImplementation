clear all;
close all;
clc

fs = 44100;
f0 = 100;               % fundamental freq
B = 0.0001;              % inharmonicity parameter (>0)x`
k = 1/fs;
T0 = 60;                % applied string tension
rho = 7700;             % material density
Area = 2.02682992e-7;   % string cross sectional area
c = sqrt(T0/(rho*Area));
h = c*k;                % calculate grid spacing
L = 1;                  % scaling lenght
gamma = 2*f0;           % scaling
durration = 1;          % synthesised sound lenght in s
dur = fs*durration;
theta = 0.575;            % free parameter for the implicit scheme
loss = [100, 10; 1000, 8]; % loss [freq.(Hz), T60(s), freq.(Hz), T60(s)]
%% 
kappa=sqrt(B)*(gamma/pi);

% set scheme for loss parameters 
% zeta0,1  pg.189 Practical setting for decay times
% sigma0,1 pg.189 eq.7.29

zeta1 = (-gamma^2+sqrt(gamma^4+4*kappa^2*(2*pi*loss(1,1))^2))/(2*kappa^2);
zeta2 = (-gamma^2+sqrt(gamma^4+4*kappa^2*(2*pi*loss(2,1))^2))/(2*kappa^2);
sigma0 = 6*log(10)*(-zeta2/loss(1,2)+zeta1/loss(2,2))/(zeta1-zeta2);
sigma1 = 6*log(10)*(1/loss(1,2)-1/loss(2,2))/(zeta1-zeta2);

h = sqrt((gamma^2 * k^2 + sqrt(gamma^4 * k^4 + 16*kappa^2 * k^2 * (2*theta-1)))/(2*(2*theta - 1))); %eq.7.26 pg.188

%h = sqrt(2*kappa*k);
N = floor(1/h);
h = 1/N;

% intialise states of the system

uNext = zeros(N,1);
u = zeros(N,1);
width = floor(N/10);
excitationRange = 1:width;
u(excitationRange + floor(N/5)) = hann(width);
uPrev = u;
plot(u)

% initialise output
out = zeros(dur,1);
outPos = floor(N/2);

l = 3:N-2;
%% loopty loop 
for n = 1:dur
    uNext(l) = (1/(sigma0*k + 1)) * (((gamma^2 * (u(l+1) - 2*u(l) + u(l-1))/h^2) - ...
        (kappa^2 * (u(l+2) - 4*u(l+1) + 6*u(l) - 4*u(l-1) + u(l-2))/h^4) + ...
        (2*sigma1 * (u(l+1) - 2*u(l) + u(l-1) - uPrev(l+1) + 2*uPrev(l) - uPrev(l-1))/k*h^2))*k^2 + ...
        sigma0*k*uPrev(l) + 2*u(l) - uPrev(l)); % eq. 7.30(a) pg.190 
    out(n) = uNext(floor(outPos));
%     plot(uNext);
%     ylim([-1,1]);
%     drawnow;
    
    uPrev  = u;
    u = uNext;
end
plot(out)