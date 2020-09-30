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
    % find virtual grid points
    uB0 = 2*uB(1)-uB(2);                    % uB(0)
    uBm1 = 2*(uB0-uB(2))+uB(3);             % uB(-1)
    uBPrev0 = 2*uBPrev(1)-uBPrev(2);        % uBPrev(0)
    uBNp1 = 2*uB(NB) - uB(NB-1);            % uB(N+1)
    uBNp2 = 2*(uBNp1 - uB(NB-1)) + uB(NB-2);% uB(N+2)
    uBPrevNp1 = 2*uBPrev(NB) - uBPrev(NB-1);% uBPrev(N+1)
    % solve for uBNext at points 1, 2, N-1, N
    uBNext(2) = 1/(hB^4 * (k*sigmaB0 + 1))*(-6*(uB(2) - 2*uB(2+1)/3 + uB(2+2)/6 - 2*uB(2-1)/3 + uB0/6)*kappaB^2 *k^2 + ...
        hB^2 * (hB^2 *uBPrev(2)*sigmaB0 +4*(uBPrev(2) - uBPrev(2+1)/2 - uBPrev(2-1)/2 - uB(2) + uB(2+1)/2 + uB(2-1)/2)*sigmaB1)*k - ...
        hB^4 * (uBPrev(2) - 2*uB(2)));
    uBNext(1) = 1/(hB^4 * (k*sigmaB0 + 1))*(-6*(uB(1) - 2*uB(1+1)/3 + uB(1+2)/6 - 2*uB0/3 + uBm1/6)*kappaB^2 *k^2 + ...
        hB^2 * (hB^2 *uBPrev(1)*sigmaB0 +4*(uBPrev(1) - uBPrev(1+1)/2 - uBPrev0/2 - uB(1) + uB(1+1)/2 + uB0/2)*sigmaB1)*k - ...
        hB^4 * (uBPrev(1) - 2*uB(1)));
    uBNext(NB-1) = 1/(hB^4 * (k*sigmaB0 + 1))*(-6*(uB(NB-1) - 2*uB(NB)/3 + uBNp1/6 - 2*uB(NB-2)/3 + uB(NB-3)/6)*kappaB^2 *k^2 + ...
        hB^2 * (hB^2 *uBPrev(NB-1)*sigmaB0 +4*(uBPrev(NB-1) - uBPrev(NB)/2 - uBPrev(NB-2)/2 - uB(NB-1) + uB(NB)/2 + uB(NB-2)/2)*sigmaB1)*k - ...
        hB^4 * (uBPrev(NB-1) - 2*uB(NB-1)));
    uBNext(NB) = 1/(hB^4 * (k*sigmaB0 + 1))*(-6*(uB(NB) - 2*uBNp1/3 + uBNp2/6 - 2*uB(NB-1)/3 + uB(NB-2)/6)*kappaB^2 *k^2 + ...
        hB^2 * (hB^2 *uBPrev(NB)*sigmaB0 +4*(uBPrev(NB) - uBPrevNp1/2 - uBPrev(NB-1)/2 - uB(NB) + uBNp1/2 + uB(NB-1)/2)*sigmaB1)*k - ...
        hB^4 * (uBPrev(NB) - 2*uB(NB)));
    
    
    
    outB(n) = uBNext(outPosB);
%     plot(uBNext);
%     ylim([-1,1]);
%     drawnow;
     
    uBPrev  = uB;
    uB = uBNext;

end
figure(2)
plot(outB)