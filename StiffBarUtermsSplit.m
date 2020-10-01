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
uBlMult = (2*hB^4 - 4*hB^2*k*sigmaB1 - 6*k^2*kappaB^2)/(hB^4*(k*sigmaB0 + 1));
uBl1Mult = (2*hB^2*k*sigmaB1 + 4*k^2*kappaB^2)/(hB^4*(k*sigmaB0 + 1));
uBl2Mult = (-k^2*kappaB^2)/(hB^4*(k*sigmaB0 + 1));
uBPrevlMult = (hB^4*k*sigmaB0 - hB^4 + 4*hB^2*k*sigmaB1)/(hB^4*(k*sigmaB0 + 1));
uBPrevl1Mult = (-2*sigmaB1*k)*1/(hB^2*(k*sigmaB0 + 1));
tic
for n = 1:dur
    uB0 = 2*uB(1)-uB(2);                    % uB(0)
    uBm1 = 2*(uB0-uB(2))+uB(3);             % uB(-1)
    uBPrev0 = 2*uBPrev(1)-uBPrev(2);        % uBPrev(0)
    uBNp1 = 2*uB(NB) - uB(NB-1);            % uB(N+1)
    uBNp2 = 2*(uBNp1 - uB(NB-1)) + uB(NB-2);% uB(N+2)
    uBPrevNp1 = 2*uBPrev(NB) - uBPrev(NB-1);% uBPrev(N+1)
    
    uBNext(lB) = uB(lB)*uBlMult + (uB(lB-1) + uB(lB+1))*uBl1Mult + (uB(lB-2) + uB(lB+2))*uBl2Mult + ...
        uBPrev(lB)*uBPrevlMult + (uBPrev(lB-1)+ uBPrev(lB+1))*uBPrevl1Mult;
    % find virtual grid points

    % solve for uBNext at points 1, 2, N-1, N
    uBNext(2) = uB(2)*uBlMult + (uB(2-1) + uB(2+1))*uBl1Mult + (uB0 + uB(2+2))*uBl2Mult + ...
        uBPrev(2)*uBPrevlMult + (uBPrev(2-1) + uBPrev(2+1))*uBPrevl1Mult;
    uBNext(1) = uB(1)*uBlMult + (uB0 + uB(1+1))*uBl1Mult + (uBm1 + uB(1+2))*uBl2Mult + ...
        uBPrev(1)*uBPrevlMult + (uBPrev0+ uBPrev(1+1))*uBPrevl1Mult;
    uBNext(NB-1) = uB(NB-1)*uBlMult + (uB(NB-1-1) + uB(NB-1+1))*uBl1Mult + (uB(NB-1-2) + uBNp1)*uBl2Mult + ...
        uBPrev(NB-1)*uBPrevlMult + (uBPrev(NB-1-1) + uBPrev(NB-1+1))*uBPrevl1Mult;
    uBNext(NB) = uB(NB)*uBlMult + (uB(NB-1) + uBNp1)*uBl1Mult + (uB(NB-2) + uBNp2)*uBl2Mult + ...
        uBPrev(NB)*uBPrevlMult + (uBPrev(NB-1) + uBPrevNp1)*uBPrevl1Mult;
        
    
    outB(n) = uBNext(outPosB);
%     plot(uBNext);
%     ylim([-1,1]);
%     drawnow;
     
    uBPrev  = uB;
    uB = uBNext;

end
toc
figure(2)
plot(outB)