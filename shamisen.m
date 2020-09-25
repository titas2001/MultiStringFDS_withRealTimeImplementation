clear all;
close all;
clc

fs = 44100;             % sampling freq
f0 = 100;               % fundamental freq
B = 0.0001;             % inharmonicity parameter (>0)x`
k = 1/fs;               % time step
T0 = 60;                % applied string tension
T = 40;                 % applied plate tension
rhoS = 7700;            % material density of the string
rhoP = 1150;            % nylon https://www.engineeringtoolbox.com/engineering-materials-properties-d_1225.html
Area = 2.02682992e-7;   % string cross sectional area
H = 0.002;              % plate thickness
cS = sqrt(T0/(rhoS*Area));
cP = sqrt(T/(rhoP*H));          
LS = 1;                 % scaling lenght
E = 3e+9;               % nylon https://www.engineeringtoolbox.com/engineering-materials-properties-d_1225.html
nu = 0.4;               % Poissonís ratio nu < 0.5
r = 1.3;                % grid aspect ratio
Lx = r*0.4;             % length of plate in x direction
Ly = (1/r)*0.4;         % length of plate in y direction
durration = 5;          % synthesised sound lenght in seconds
dur = fs*durration;     % synthesised sound lenght in samples
gammaS = 2*f0;          % scaling for a string
gammaP = sqrt(T/(rhoP*H*Lx*Ly));  
theta = 0.575;          % free parameter for the implicit scheme
R = 1;                  % resistance in mass spring
lossS = [100, 10; 1000, 8]; % loss [freq.(Hz), T60(s), freq.(Hz), T60(s)]
lossP = [100, 10; 1000, 8]; % loss [freq.(Hz), T60(s), freq.(Hz), T60(s)]
vP0 = -10;
M = 0.03;
%
kappaS=sqrt(B)*(gammaS/pi);
D = E*H^3 / (12 * (1 - nu^2)); % plate flexural rigidity pg.341
kappaPsq = D / (rhoP * H * Lx^2 * Ly^2); % pg.342 eq.12.3 kappa^2


% zeta0,1  pg.189 Practical setting for decay times
% sigma0,1 pg.189 eq.7.29
% set scheme for loss parameters for String
zetaS1 = (-gammaS^2+sqrt(gammaS^4+4*kappaS^2*(2*pi*lossS(1,1))^2))/(2*kappaS^2);
zetaS2 = (-gammaS^2+sqrt(gammaS^4+4*kappaS^2*(2*pi*lossS(2,1))^2))/(2*kappaS^2);
sigmaS0 = 6*log(10)*(-zetaS2/lossS(1,2)+zetaS1/lossS(2,2))/(zetaS1-zetaS2);
sigmaS1 = 6*log(10)*(1/lossS(1,2)-1/lossS(2,2))/(zetaS1-zetaS2);
% set scheme for loss parameters for Plate
zetaP1 = (-gammaP^2+sqrt(gammaP^4+4*kappaPsq*(2*pi*lossP(1,1))^2))/(2*kappaPsq);
zetaP2 = (-gammaP^2+sqrt(gammaP^4+4*kappaPsq*(2*pi*lossP(2,1))^2))/(2*kappaPsq);
sigmaP0 = 6*log(10)*(-zetaP2/lossP(1,2)+zetaP1/lossP(2,2))/(zetaP1-zetaP2);
sigmaP1 = 6*log(10)*(1/lossP(1,2)-1/lossP(2,2))/(zetaP1-zetaP2);

hS = sqrt((gammaS^2 * k^2 + sqrt(gammaS^4 * k^4 + 16*kappaS^2 * k^2 * (2*theta-1)))/(2*(2*theta - 1))); % set grid spacing for String eq.7.26 pg.188
hP = sqrt((cP^2 * k^2 + 4*sigmaP1*k + sqrt((cP^2 * k^2 + 4*sigmaP1*k)^2 + 16*kappaPsq*k^2))); % set grid spacing for Plate tromba marina paper eq. 20

NS = floor(1/hS);
Nx = floor(sqrt(r)/hP);        % number of x-subdivisions of spatial domain
Ny = floor(1/(sqrt(r)*hP));    % number of y-subdivisions of spatial domain
hP = sqrt(r)/min(Nx, Ny);      % reset grid spacing for Plate
hS = 1/NS;


% intialise states of the system

uSNext = zeros(NS,1);
uS = zeros(NS,1);
width = floor(NS/10);
excitationRange = 1:width;
uS(excitationRange + floor(NS/5)) = hann(width);
uSPrev = uS;

uPNext = zeros(Nx,Ny);
uP = zeros(Nx,Ny);
uPPrev = zeros(Nx,Ny);
% uP(ceil(Nx/2-Nx/8):floor(Nx/2+Nx/8),ceil(Ny/2-Ny/8):floor(Ny/2+Ny/8)) = ... here plate is excited by velocity
%     k*vP0*hamming_3d(length(ceil(Nx/2-Nx/8):floor(Nx/2+Nx/8)),length(ceil(Ny/2-Ny/8):floor(Ny/2+Ny/8)),1);
outPos = floor(NS/2);

uMNext = 0;
uM = zeros(100,1);
uM(1) = 1;
uMPrev = uM;

out = zeros(dur,1);

lP = 3:Nx-2;
mP = 3:Ny-2;
lS = 3:NS-2;
lM = 1;
lMc = 1;
lSc = NS - floor(NS/8);
lPc = Nx - floor(Nx/3);
mPc = Ny - floor(Ny/3);
tic
for n = 1:dur
    % Force from string to mass spring (bridge)
    Fsm = (1/(hS^3*(sigmaS0*k + 1)*(rhoS*Area*hS + M)))*(4*Area*rhoS*M * ((sigmaS0*uM(lMc)*k^3 * pi^2 *f0^2 + ...
        ((R*(uM(lMc)-uMPrev(lMc))*sigmaS0)/4 + pi^2 * f0^2 * uM(lMc))*k^2 + ...
        ((-uM(lMc)/2 + uSPrev(lSc)/4 + uMPrev(lMc)/4)*sigmaS0 + (R*(uM(lMc) - uMPrev(lMc)))/4)*k - ...
        uM(lMc)/2 + uS(lSc)/2 - uSPrev(lSc)/4 + uMPrev(lMc)/4)*hS^4 - ...
        ((gammaS^2 * (uS(lSc) - uS(lSc+1)/2 - uS(lSc-1)/2)*k)/2 + sigmaS1*(uS(lSc) - uS(lSc+1)/2 - uS(lSc-1)/2 - ...
        uSPrev(lSc) + uSPrev(lSc+1)/2 + uSPrev(lSc-1)/2))*k*hS^2 - (3* (uS(lSc) - 2*uS(lSc+1)/3 + uS(lSc+2)/6 - ...
        2*uS(lSc-1)/3 + uS(lSc-2)/6)*kappaS^2 * k^2)/2));
    
    % Force from mass spring (bridge) to plate
    Fmp = (-1/(hP^2 * (k*sigmaP0 + 1)*(rhoP*H*hP^2 + M)))*(4*H*rhoP*(pi^2 * hP^4 * k^3 * sigmaP0 * f0^2 * uM(lMc) + ...
        (((R* (uM(lMc) - uMPrev(lMc))*sigmaP0)/4 + pi^2 * f0^2 * uM(lMc))*hP^4 - ...
        gammaP^2 * (uP(lPc,mPc) - uP(lPc+1,mPc)/4 - uP(lPc,mPc+1)/4 - uP(lPc-1,mPc)/4 - uP(lPc,mPc-1)/4)*hP^2 - ...
        1/4 * (kappaPsq * (uP(lPc+2,mPc) + uP(lPc-2,mPc) + uP(lPc,mPc+2) + uP(lPc,mPc-2) + ...
        2*(uP(lPc+1,mPc+1) + uP(lPc+1,mPc-1) + uP(lPc-1,mPc+1) + uP(lPc-1,mPc-1)) - ...
        8*(uP(lPc+1,mPc) + uP(lPc-1,mPc) + uP(lPc,mPc+1) + uP(lPc,mPc-1)) + 20*uP(lPc, mPc) )))*k^2 - ...
        1/2 *(hP^2 * (((uM(lMc) - uMPrev(lMc)/2 -uPPrev(lPc, mPc)/2)*sigmaP0 - (R*(uM(lMc)-uMPrev(lMc)))/2)*hP^2 + ...
        4*sigmaP1*(uP(lPc,mPc) - uP(lPc+1,mPc)/4 - uP(lPc,mPc+1)/4 - uP(lPc-1,mPc)/4 - uP(lPc,mPc-1)/4 - ...
        uPPrev(lPc,mPc) + uPPrev(lPc+1,mPc)/4 + uPPrev(lPc,mPc+1)/4 + uPPrev(lPc-1,mPc)/4 + uPPrev(lPc,mPc-1)/4))*k)-...
        ((uM(lMc) - uP(lPc, mPc) - uMPrev(lMc)/2 + uPPrev(lPc,mPc)/2)*hP^4)/2)*M);

    % Update function of the String
    uSNext(lS) = (1/(sigmaS0*k + 1)) * (((gammaS^2 * (uS(lS+1) - 2*uS(lS) + uS(lS-1))/hS^2) - ...
        (kappaS^2 * (uS(lS+2) - 4*uS(lS+1) + 6*uS(lS) - 4*uS(lS-1) + uS(lS-2))/hS^4) + ...
        (2*sigmaS1 * (uS(lS+1) - 2*uS(lS) + uS(lS-1) - uSPrev(lS+1) + 2*uSPrev(lS) - uSPrev(lS-1))/k*hS^2))*k^2 + ...
        sigmaS0*k*uSPrev(lS) + 2*uS(lS) - uSPrev(lS)); % eq. 7.30(a) pg.190 
    % Update String function at the localizer point with Fsm (Force loss to the bridge)
    uSNext(lSc) = uSNext(lSc) - Fsm/(rhoS * Area *hS);
    
    % Update function of the Mass spring (bridge)
    uMNext(lM) = -4*pi^2*f0^2*k^2*uM(lM) - R*k*(uM(lM) - uMPrev(lM)) + 2*uM(lM) - uMPrev(lM);
    % Update Mass spring function at the localizer point with Fsm (Force gain from the string)
    uMNext(lMc) = uMNext(lMc) + Fsm/M;
    % Update Mass spring function at the localizer point with Fmp (Force loss to the plate)
    uMNext(lMc) = uMNext(lMc) - Fmp/M;
    
    % Update function of the Plate
    uPNext(lP,mP) = (1/(k*sigmaP0 + 1))*(((-(kappaPsq)/hP^4)*((uP(lP+2,mP) + uP(lP-2,mP) + uP(lP,mP+2) + uP(lP,mP-2)) + ...
        2*(uP(lP+1,mP+1) + uP(lP+1,mP-1) + uP(lP-1,mP+1) + uP(lP-1,mP-1)) - ...
        8*(uP(lP+1,mP) + uP(lP-1,mP) + uP(lP,mP+1) + uP(lP,mP-1)) + 20*uP(lP,mP)) + ...
        (gammaP^2/hP^2)*(uP(lP+1,mP) + uP(lP-1,mP) + uP(lP,mP+1) + uP(lP,mP-1) - 4*uP(lP,mP)) + ...
        ((2*sigmaP1)/(k*hP^2))*(uP(lP+1,mP) + uP(lP-1,mP) + uP(lP,mP+1) + uP(lP,mP-1) - 4*uP(lP,mP) - ...
        (uPPrev(lP+1,mP) + uPPrev(lP-1,mP) + uPPrev(lP,mP+1) + uPPrev(lP,mP-1) - 4*uPPrev(lP,mP))))*k^2 + ...
        k*sigmaP0*uPPrev(lP,mP) + 2*uP(lP,mP) - uPPrev(lP,mP));
    % Update Plate function at the localizer point with Fmp (Force gain from the bridge)
    uPNext(lPc,mPc) = uPNext(lPc,mPc) + Fmp/(rhoP*H*hP^2);
    
    out(n) = uPNext(floor(Nx/2),floor(Ny/2));           % plate
    out(n) = out(n) + uMNext(lM);                       % mass spring
    out(n) = out(n) + uSNext(outPos);                   % spring
    % update the state variables
    uSPrev  = uS;
    uS = uSNext;    
    uPPrev  = uP;
    uP = uPNext;
    uMPrev  = uM;
    uM = uMNext;
end
toc
plot(out);