clear all;
close all;
clc

fs = 44100;
f0 = 100;               % fundamental freq
k = 1/fs;
durration = 1;          % synthesised sound lenght in s
dur = fs*durration;

Bb = 0.0001;             % inharmonicity parameter (>0)x`
rhoB = 800;              % material density
AreaB = 2.02e-4;         % bridge cross sectional area
LB = 1;                % scaling lenght
EB = 9.5e+9;           % Young's modulus dried Red Alder https://amesweb.info/Materials/Youngs-Modulus-of-Wood.aspx
HB = 0.075;             % thickness
 
kappaB=sqrt((EB*HB^2)/(12*rhoB*LB^4)); % eq. 7.70 pg. 210

% set scheme for loss parameters 
% zeta0,1  pg.189 Practical setting for decay times
% sigma0,1 pg.189 eq.7.29


sigmaB0 = 1.343;
sigmaB1 = 0.00459;

hB = sqrt((4*sigmaB1*k+sqrt((4*sigmaB1*k)^2+16*kappaB^2*k^2))/2);

%h = sqrt(2*kappa*k);
NB = floor(1/hB);
hB = 1/NB;

% intialise states of the system

uBNext = zeros(NB,1);
uB = zeros(NB,1);
widthB = floor(NB/2);
excitationRangeB = 1:widthB;
uB(excitationRangeB + floor(NB/2)) = hann(widthB);
uBPrev = uB;

% initialise output
outB = zeros(dur,1);
outPosB = floor(NB/2);

lB = 3:NB-2;
% figure(3)
% loopty loop 
for n = 1:dur

    
    uBNext(lB) = 1/(hB^4 * (k*sigmaB0 + 1))*(-6*(uB(lB) - 2*uB(lB+1)/3 + uB(lB+2)/6 - 2*uB(lB-1)/3 + uB(lB-2)/6)*kappaB^2 *k^2 + ...
        hB^2 * (hB^2 *uBPrev(lB)*sigmaB0 +4*(uBPrev(lB) - uBPrev(lB+1)/2 - uBPrev(lB-1)/2 - uB(lB) + uB(lB+1)/2 + uB(lB-1)/2)*sigmaB1)*k - ...
        hB^4 * (uBPrev(lB) - 2*uB(lB)));
    % update at boundaries
    uBNext(1) = (-hB*uBNext(3) + 6*uBNext(3) - 4*uBNext(4))/(hB + 2);
    uBNext(2) = -uBNext(4) + 2*uBNext(3);
    uBNext(NB) = (-hB*uBNext(NB-2) + 6*uBNext(NB-2) - 4*uBNext(NB-3))/(hB + 2);
    uBNext(NB-1) = -uBNext(NB-3) + 2*uBNext(NB-2);
    
    outB(n) = uBNext(outPosB);
    plot(uBNext);
    ylim([-1,1]);
    drawnow;
%     
    uBPrev  = uB;
    uB(lB) = uBNext(lB);

end
figure(2)
plot(outB)