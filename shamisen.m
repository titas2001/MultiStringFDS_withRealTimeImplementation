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
AreaS = 2.02682992e-7;   % string cross sectional area
HP = 0.002;              % plate thickness
cS = sqrt(T0/(rhoS*AreaS));
cP = sqrt(T/(rhoP*HP));          
LS = 1;                 % scaling lenght
E = 3e+9;               % nylon https://www.engineeringtoolbox.com/engineering-materials-properties-d_1225.html
nu = 0.4;               % Poisson’s ratio nu < 0.5
r = 1.3;                % grid aspect ratio
Lx = r*0.4;             % length of plate in x direction
Ly = (1/r)*0.4;         % length of plate in y direction
durration = 5;          % synthesised sound lenght in seconds
dur = fs*durration;     % synthesised sound lenght in samples
gammaS = 2*f0;          % scaling for a string
gammaP = sqrt(T/(rhoP*HP*Lx*Ly));  
theta = 0.575;          % free parameter for the implicit scheme
R = 1;                  % resistance in mass spring
lossS = [100, 10; 1000, 8]; % loss [freq.(Hz), T60(s), freq.(Hz), T60(s)]
lossP = [100, 10; 1000, 8]; % loss [freq.(Hz), T60(s), freq.(Hz), T60(s)]
vP0 = -10;
M = 0.03;
Bb = 0.0001;             % inharmonicity parameter (>0)x`
rhoB = 800;              % material density
AreaB = 2.02e-4;         % bridge cross sectional area
LB = 1;                % scaling lenght
EB = 9.5e+9;           % Young's modulus dried Red Alder https://amesweb.info/Materials/Youngs-Modulus-of-Wood.aspx
HB = 0.075;             % thickness
 
kappaB=sqrt((EB*HB^2)/(12*rhoB*LB^4)); % eq. 7.70 pg. 210

%
kappaS=sqrt(B)*(gammaS/pi);
D = E*HP^3 / (12 * (1 - nu^2)); % plate flexural rigidity pg.341
kappaPsq = D / (rhoP * HP * Lx^2 * Ly^2); % pg.342 eq.12.3 kappa^2


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
% loss parameters for Bar
sigmaB0 = 1.343;
sigmaB1 = 0.00459;

hB = sqrt((4*sigmaB1*k+sqrt((4*sigmaB1*k)^2+16*kappaB^2*k^2))/2);
hS = sqrt((gammaS^2 * k^2 + sqrt(gammaS^4 * k^4 + 16*kappaS^2 * k^2 * (2*theta-1)))/(2*(2*theta - 1))); % set grid spacing for String eq.7.26 pg.188
hP = sqrt((cP^2 * k^2 + 4*sigmaP1*k + sqrt((cP^2 * k^2 + 4*sigmaP1*k)^2 + 16*kappaPsq*k^2))); % set grid spacing for Plate tromba marina paper eq. 20

NS = floor(1/hS);              % string spatial subdivisions
Nx = floor(sqrt(r)/hP);        % number of x-subdivisions of spatial domain
Ny = floor(1/(sqrt(r)*hP));    % number of y-subdivisions of spatial domain
NB = floor(1/hB);              % bar spatial subdivisions
hP = sqrt(r)/min(Nx, Ny);      % reset grid spacing for Plate
hS = 1/NS;                     % reset grid spacing for String
hB = 1/NB;                     % reset grid spacing for Bar

%% Intialise states of the system

% Strings
uS1Next = zeros(NS,1);
uS1 = zeros(NS,1);
width = floor(NS/10);
excitationRange = 1:width;
uS1(excitationRange + floor(NS/5)) = hann(width);
uS1Prev = uS1;
uS2Next = zeros(NS,1);
uS2 = zeros(NS,1);
uS2Prev = uS2;
uS3Next = zeros(NS,1);
uS3 = zeros(NS,1);
uS3Prev = uS3;


% Plate
uPNext = zeros(Nx,Ny);
uP = zeros(Nx,Ny);
uPPrev = zeros(Nx,Ny);
% uP(ceil(Nx/2-Nx/8):floor(Nx/2+Nx/8),ceil(Ny/2-Ny/8):floor(Ny/2+Ny/8)) = ... here plate is excited by velocity
%     k*vP0*hamming_3d(length(ceil(Nx/2-Nx/8):floor(Nx/2+Nx/8)),length(ceil(Ny/2-Ny/8):floor(Ny/2+Ny/8)),1);
outPos = floor(NS/2);

% Bar
uBNext = zeros(NB,1);
uB = zeros(NB,1);
uBPrev = zeros(NB,1);

% Mass Spring
uMNext = 0;
uM = zeros(100,1);
uM(1) = 1;
uMPrev = uM;


% Output
out = zeros(dur,1);
%% Intialise l for update equations
lP = 3:Nx-2;
mP = 3:Ny-2;
lS = 3:NS-2;
lB = 3:NB-2;

%% Connection points
lBc1 = 5;       % bar connection to the 1st string
lBc2 = 9;       % bar connection to the 2nd string
lBc3 = 13;      % bar connection to the 3rd string
lBcl = 1;       % bar left side connection to the plate
lBcr = 17;      % bar right side connection to the plate

lM = 1;         % will be deprecated
lMc = 1;        % will be deprecated

lS1c = NS - floor(NS/8); % 1st string connection to the bar
lS2c = NS - floor(NS/8); % 2nd string connection to the bar
lS3c = NS - floor(NS/8); % 3rd string connection to the bar

lPc = Nx - floor(Nx/3); % will be deprecated
mPc = Ny - floor(Ny/3); % will be deprecated

lPcl = Nx - floor(2*Nx/5); % Plate connection to the bar on the left side x coordinate
lPcr = Nx - floor(3*Nx/5); % Plate connection to the bar on the right side x coordinate
mPcl = Ny - floor(Ny/4);   % Plate connection to the bar on the left side y coordinate
mPcr = Ny - floor(Ny/4);   % Plate connection to the bar on the right side y coordinate

%% Multipliers

% String
uSlMult = (((-2*gammaS^2)/hS^2 - 6*kappaS^2/hS^4 - 4*sigmaS1/(k*hS^2))*k^2 + 2)/(k*sigmaS0 + 1);
uSl1Mult = (gammaS^2/hS^2 + 4*kappaS^2/hS^4 + 2*sigmaS1/(k*hS^2))*k^2/(k*sigmaS0 + 1);
uSl2Mult = -1*k^2*kappaS^2/((k*sigmaS0 + 1)*hS^4);
uSPrevlMult = ((4*sigmaS1*k^2)/(k*hS^2) + k*sigmaS0 - 1)/(k*sigmaS0 + 1);
uSPrevl1Mult = ((-2*sigmaS1*k^2)*1/(k*hS^2))/(k*sigmaS0 + 1);

% Bar
uBlMult = (2*hB^4 - 4*hB^2*k*sigmaB1 - 6*k^2*kappaB^2)/(hB^4*(k*sigmaB0 + 1));
uBl1Mult = (2*hB^2*k*sigmaB1 + 4*k^2*kappaB^2)/(hB^4*(k*sigmaB0 + 1));
uBl2Mult = (-k^2*kappaB^2)/(hB^4*(k*sigmaB0 + 1));
uBPrevlMult = (hB^4*k*sigmaB0 - hB^4 + 4*hB^2*k*sigmaB1)/(hB^4*(k*sigmaB0 + 1));
uBPrevl1Mult = (-2*sigmaB1*k)*1/(hB^2*(k*sigmaB0 + 1));

% Plate
uPlMult = ((-20*kappaPsq/hP^4 - 4*gammaP^2/hP^2 - 8*sigmaP1/(k*hP^2))*k^2 + 2)/(k*sigmaP0 + 1); 
uPl1Mult = (8*kappaPsq/hP^4 + gammaP^2/hP^2 + 2*sigmaP1/(k*hP^2))*k^2/(k*sigmaP0 + 1);
uPl1dMult = (-2*kappaPsq*k^2)*1/hP^4/(k*sigmaP0 + 1);
uPl2Mult = -kappaPsq*k^2/(hP^4*(k*sigmaP0 + 1));
uPPrevlMult = ((8*sigmaP1*k^2)*1/(k*hP^2) + k*sigmaP0 - 1)/(k*sigmaP0 + 1);
uPPrevl1Mult = (-2*sigmaP1*k^2)*1/(k*hP^2)/(k*sigmaP0 + 1);

% Forces
FsbMult = 1/(1/(rhoB*AreaB*hB) + 1/(rhoS*AreaS*hS));
FbpMult = 1/(-1/(rhoB*AreaB*hB) - 1/(rhoP*HP*hP^2));

%%
tic
for n = 1:dur
        % calculate virtual grid points
    uB0 = 2*uB(1)-uB(2);                    % uB(0)
    uBm1 = 2*(uB0-uB(2))+uB(3);             % uB(-1)
    uBPrev0 = 2*uBPrev(1)-uBPrev(2);        % uBPrev(0)
    uBNp1 = 2*uB(NB) - uB(NB-1);            % uB(N+1)
    uBNp2 = 2*(uBNp1 - uB(NB-1)) + uB(NB-2);% uB(N+2)
    uBPrevNp1 = 2*uBPrev(NB) - uBPrev(NB-1);% uBPrev(N+1)
    % Force from string to mass spring (bridge) 3x
%     Fsm = (1/(hS^3*(sigmaS0*k + 1)*(rhoS*AreaS*hS + M)))*(4*AreaS*rhoS*M * ((sigmaS0*uM(lMc)*k^3 * pi^2 *f0^2 + ...
%         ((R*(uM(lMc)-uMPrev(lMc))*sigmaS0)/4 + pi^2 * f0^2 * uM(lMc))*k^2 + ...
%         ((-uM(lMc)/2 + uSPrev(lS1c)/4 + uMPrev(lMc)/4)*sigmaS0 + (R*(uM(lMc) - uMPrev(lMc)))/4)*k - ...
%         uM(lMc)/2 + uS(lS1c)/2 - uSPrev(lS1c)/4 + uMPrev(lMc)/4)*hS^4 - ...
%         ((gammaS^2 * (uS(lS1c) - uS(lS1c+1)/2 - uS(lS1c-1)/2)*k)/2 + sigmaS1*(uS(lS1c) - uS(lS1c+1)/2 - uS(lS1c-1)/2 - ...
%         uSPrev(lS1c) + uSPrev(lS1c+1)/2 + uSPrev(lS1c-1)/2))*k*hS^2 - (3* (uS(lS1c) - 2*uS(lS1c+1)/3 + uS(lS1c+2)/6 - ...
%         2*uS(lS1c-1)/3 + uS(lS1c-2)/6)*kappaS^2 * k^2)/2));
    % Force from strings to the bridge
    Fsb1 = FsbMult * (uB(lBc1)*uBlMult + (uB(lBc1-1) + uB(lBc1+1))*uBl1Mult + (uB(lBc1-2) + uB(lBc1+2))*uBl2Mult + ...
        uBPrev(lBc1)*uBPrevlMult + (uBPrev(lBc1-1)+ uBPrev(lBc1+1))*uBPrevl1Mult + ...
        uS1(lS1c)*uSlMult + (uS1(lS1c-1) + uS1(lS1c+1))*uSl1Mult + (uS1(lS1c-2) + uS1(lS1c+2))*uSl2Mult + ...
        uS1Prev(lS1c)*uSPrevlMult + (uS1Prev(lS1c-1)+ uS1Prev(lS1c+1))*uSPrevl1Mult);
    
    Fsb2 = FsbMult * (uB(lBc2)*uBlMult + (uB(lBc2-1) + uB(lBc2+1))*uBl1Mult + (uB(lBc2-2) + uB(lBc2+2))*uBl2Mult + ...
        uBPrev(lBc2)*uBPrevlMult + (uBPrev(lBc2-1)+ uBPrev(lBc2+1))*uBPrevl1Mult + ...
        uS2(lS2c)*uSlMult + (uS2(lS2c-1) + uS2(lS2c+1))*uSl1Mult + (uS2(lS2c-2) + uS2(lS2c+2))*uSl2Mult + ...
        uS2Prev(lS2c)*uSPrevlMult + (uS2Prev(lS2c-1)+ uS2Prev(lS2c+1))*uSPrevl1Mult);
    
    Fsb3 = FsbMult * (uB(lBc3)*uBlMult + (uB(lBc3-1) + uB(lBc3+1))*uBl1Mult + (uB(lBc3-2) + uB(lBc3+2))*uBl2Mult + ...
        uBPrev(lBc3)*uBPrevlMult + (uBPrev(lBc3-1)+ uBPrev(lBc3+1))*uBPrevl1Mult + ...
        uS3(lS3c)*uSlMult + (uS3(lS3c-1) + uS3(lS3c+1))*uSl1Mult + (uS3(lS3c-2) + uS3(lS3c+2))*uSl2Mult + ...
        uS3Prev(lS3c)*uSPrevlMult + (uS3Prev(lS3c-1)+ uS3Prev(lS3c+1))*uSPrevl1Mult);
    
    % Force from mass spring (bridge) to plate 2x
%     Fmp = (-1/(hP^2 * (k*sigmaP0 + 1)*(rhoP*H*hP^2 + M)))*(4*H*rhoP*(pi^2 * hP^4 * k^3 * sigmaP0 * f0^2 * uM(lMc) + ...
%         (((R* (uM(lMc) - uMPrev(lMc))*sigmaP0)/4 + pi^2 * f0^2 * uM(lMc))*hP^4 - ...
%         gammaP^2 * (uP(lPc,mPc) - uP(lPc+1,mPc)/4 - uP(lPc,mPc+1)/4 - uP(lPc-1,mPc)/4 - uP(lPc,mPc-1)/4)*hP^2 - ...
%         1/4 * (kappaPsq * (uP(lPc+2,mPc) + uP(lPc-2,mPc) + uP(lPc,mPc+2) + uP(lPc,mPc-2) + ...
%         2*(uP(lPc+1,mPc+1) + uP(lPc+1,mPc-1) + uP(lPc-1,mPc+1) + uP(lPc-1,mPc-1)) - ...
%         8*(uP(lPc+1,mPc) + uP(lPc-1,mPc) + uP(lPc,mPc+1) + uP(lPc,mPc-1)) + 20*uP(lPc, mPc) )))*k^2 - ...
%         1/2 *(hP^2 * (((uM(lMc) - uMPrev(lMc)/2 -uPPrev(lPc, mPc)/2)*sigmaP0 - (R*(uM(lMc)-uMPrev(lMc)))/2)*hP^2 + ...
%         4*sigmaP1*(uP(lPc,mPc) - uP(lPc+1,mPc)/4 - uP(lPc,mPc+1)/4 - uP(lPc-1,mPc)/4 - uP(lPc,mPc-1)/4 - ...
%         uPPrev(lPc,mPc) + uPPrev(lPc+1,mPc)/4 + uPPrev(lPc,mPc+1)/4 + uPPrev(lPc-1,mPc)/4 + uPPrev(lPc,mPc-1)/4))*k)-...
%         ((uM(lMc) - uP(lPc, mPc) - uMPrev(lMc)/2 + uPPrev(lPc,mPc)/2)*hP^4)/2)*M);
    
    % Force from bridge' left and right mounting points to the plate
    Fbpl = FbpMult*(uB(lBcl)*uBlMult + (uB0 + uB(lBcl+1))*uBl1Mult + (uBm1 + uB(lBcl+2))*uBl2Mult + ...
        uBPrev(lBcl)*uBPrevlMult + (uBPrev0+ uBPrev(lBcl+1))*uBPrevl1Mult + ...
        uP(lPcl,mPcl)*uPlMult + (uP(lPcl-1,mPcl) + uP(lPcl+1,mPcl) + uP(lPcl,mPcl+1) + uP(lPcl,mPcl-1))*uPl1Mult + ...
        (uP(lPcl-1,mPcl-1) + uP(lPcl+1,mPcl-1) + uP(lPcl-1,mPcl+1) + uP(lPcl+1,mPcl+1))*uPl1dMult + ...
        (uP(lPcl-2,mPcl) + uP(lPcl+2,mPcl) + uP(lPcl,mPcl-2) + uP(lPcl,mPcl+2))*uPl2Mult + ...
        uPPrev(lPcl,mPcl)*uPPrevlMult + (uPPrev(lPcl+1,mPcl) + uPPrev(lPcl-1,mPcl) + ...
        uPPrev(lPcl,mPcl+1) + uPPrev(lPcl,mPcl-1))*uPPrevl1Mult);
    Fbpr = FbpMult*(uB(lBcr)*uBlMult + (uB(lBcr-1) + uBNp1)*uBl1Mult + (uB(lBcr-2) + uBNp2)*uBl2Mult + ...
        uBPrev(lBcr)*uBPrevlMult + (uBPrev(lBcr-1)+ uBPrevNp1)*uBPrevl1Mult + ...
        uP(lPcr,mPcr)*uPlMult + (uP(lPcr-1,mPcr) + uP(lPcr+1,mPcr) + uP(lPcr,mPcr+1) + uP(lPcr,mPcr-1))*uPl1Mult + ...
        (uP(lPcr-1,mPcr-1) + uP(lPcr+1,mPcr-1) + uP(lPcr-1,mPcr+1) + uP(lPcr+1,mPcr+1))*uPl1dMult + ...
        (uP(lPcr-2,mPcr) + uP(lPcr+2,mPcr) + uP(lPcr,mPcr-2) + uP(lPcr,mPcr+2))*uPl2Mult + ...
        uPPrev(lPcr,mPcr)*uPPrevlMult + (uPPrev(lPcr+1,mPcr) + uPPrev(lPcr-1,mPcr) + ...
        uPPrev(lPcr,mPcr+1) + uPPrev(lPcr,mPcr-1))*uPPrevl1Mult);
    
    
    % UPDATE EQUATIONS
    
    % Update equation for the Strings
%     uSNext(lS) = (1/(sigmaS0*k + 1)) * (((gammaS^2 * (uS(lS+1) - 2*uS(lS) + uS(lS-1))/hS^2) - ...
%         (kappaS^2 * (uS(lS+2) - 4*uS(lS+1) + 6*uS(lS) - 4*uS(lS-1) + uS(lS-2))/hS^4) + ...
%         (2*sigmaS1 * (uS(lS+1) - 2*uS(lS) + uS(lS-1) - uSPrev(lS+1) + 2*uSPrev(lS) - uSPrev(lS-1))/k*hS^2))*k^2 + ...
%         sigmaS0*k*uSPrev(lS) + 2*uS(lS) - uSPrev(lS)); % eq. 7.30(a) pg.190 
    uS1Next(lS) = uS1(lS)*uSlMult + (uS1(lS-1) + uS1(lS+1))*uSl1Mult + (uS1(lS-2) + uS1(lS+2))*uSl2Mult + ...
        uS1Prev(lS)*uSPrevlMult + (uS1Prev(lS-1)+ uS1Prev(lS+1))*uSPrevl1Mult;
    uS2Next(lS) = uS2(lS)*uSlMult + (uS2(lS-1) + uS2(lS+1))*uSl1Mult + (uS2(lS-2) + uS2(lS+2))*uSl2Mult + ...
        uS2Prev(lS)*uSPrevlMult + (uS2Prev(lS-1)+ uS2Prev(lS+1))*uSPrevl1Mult;
    uS3Next(lS) = uS3(lS)*uSlMult + (uS3(lS-1) + uS3(lS+1))*uSl1Mult + (uS3(lS-2) + uS3(lS+2))*uSl2Mult + ...
        uS3Prev(lS)*uSPrevlMult + (uS3Prev(lS-1)+ uS3Prev(lS+1))*uSPrevl1Mult;
    % Update Strings equation at the localizer point with Fsm (Force loss to the bridge)
    uS1Next(lS1c) = uS1Next(lS1c) - Fsb1/(rhoS * AreaS *hS);
    uS2Next(lS2c) = uS2Next(lS2c) - Fsb2/(rhoS * AreaS *hS);
    uS3Next(lS3c) = uS3Next(lS3c) - Fsb3/(rhoS * AreaS *hS);
    
    % Update equation of the Mass spring (bridge)
%     uMNext(lM) = -4*pi^2*f0^2*k^2*uM(lM) - R*k*(uM(lM) - uMPrev(lM)) + 2*uM(lM) - uMPrev(lM);
    % Update Mass spring equation at the localizer point with Fsm (Force gain from the string)
%     uMNext(lMc) = uMNext(lMc) + Fsm/M;
    % Update Mass spring equation at the localizer point with Fmp (Force loss to the plate)
%     uMNext(lMc) = uMNext(lMc) - Fmp/M;
    
    %% Update equation of the Bar (bridge)
   uBNext(lB) = uB(lB)*uBlMult + (uB(lB-1) + uB(lB+1))*uBl1Mult + (uB(lB-2) + uB(lB+2))*uBl2Mult + ...
        uBPrev(lB)*uBPrevlMult + (uBPrev(lB-1)+ uBPrev(lB+1))*uBPrevl1Mult;
    % solve for uBNext at points 1, 2, N-1, N
    uBNext(2) = uB(2)*uBlMult + (uB(2-1) + uB(2+1))*uBl1Mult + (uB0 + uB(2+2))*uBl2Mult + ...
        uBPrev(2)*uBPrevlMult + (uBPrev(2-1) + uBPrev(2+1))*uBPrevl1Mult;
    uBNext(1) = uB(1)*uBlMult + (uB0 + uB(1+1))*uBl1Mult + (uBm1 + uB(1+2))*uBl2Mult + ...
        uBPrev(1)*uBPrevlMult + (uBPrev0+ uBPrev(1+1))*uBPrevl1Mult;
    uBNext(NB-1) = uB(NB-1)*uBlMult + (uB(NB-1-1) + uB(NB-1+1))*uBl1Mult + (uB(NB-1-2) + uBNp1)*uBl2Mult + ...
        uBPrev(NB-1)*uBPrevlMult + (uBPrev(NB-1-1) + uBPrev(NB-1+1))*uBPrevl1Mult;
    uBNext(NB) = uB(NB)*uBlMult + (uB(NB-1) + uBNp1)*uBl1Mult + (uB(NB-2) + uBNp2)*uBl2Mult + ...
        uBPrev(NB)*uBPrevlMult + (uBPrev(NB-1) + uBPrevNp1)*uBPrevl1Mult;
    % Update Bar function at the localizer point with Fs1b (Force gain from the string 1)
    uBNext(lBc1) = uBNext(lBc1) + Fsb1/(rhoB * AreaB *hB);
    % Update Bar function at the localizer point with Fs2b (Force gain from the string 2)
    uBNext(lBc2) = uBNext(lBc2) + Fsb2/(rhoB * AreaB *hB);
    % Update Bar function at the localizer point with Fs3b (Force gain from the string 3)
    uBNext(lBc3) = uBNext(lBc3) + Fsb3/(rhoB * AreaB *hB);
    % Update Bar function at the localizer point lBl with Fbp (Force gain from the string 1)
    uBNext(lBcl) = uBNext(lBcl) - Fbpl/(rhoB * AreaB *hB);
    % Update Bar function at the localizer point lBr with Fbp (Force gain from the string 1)
    uBNext(lBcr) = uBNext(lBcr) - Fbpr/(rhoB * AreaB *hB);

    
    %% Update function of the Plate
%     uPNext(lP,mP) = (1/(k*sigmaP0 + 1))*(((-(kappaPsq)/hP^4)*((uP(lP+2,mP) + uP(lP-2,mP) + uP(lP,mP+2) + uP(lP,mP-2)) + ...
%         2*(uP(lP+1,mP+1) + uP(lP+1,mP-1) + uP(lP-1,mP+1) + uP(lP-1,mP-1)) - ...
%         8*(uP(lP+1,mP) + uP(lP-1,mP) + uP(lP,mP+1) + uP(lP,mP-1)) + 20*uP(lP,mP)) + ...
%         (gammaP^2/hP^2)*(uP(lP+1,mP) + uP(lP-1,mP) + uP(lP,mP+1) + uP(lP,mP-1) - 4*uP(lP,mP)) + ...
%         ((2*sigmaP1)/(k*hP^2))*(uP(lP+1,mP) + uP(lP-1,mP) + uP(lP,mP+1) + uP(lP,mP-1) - 4*uP(lP,mP) - ...
%         (uPPrev(lP+1,mP) + uPPrev(lP-1,mP) + uPPrev(lP,mP+1) + uPPrev(lP,mP-1) - 4*uPPrev(lP,mP))))*k^2 + ...
%         k*sigmaP0*uPPrev(lP,mP) + 2*uP(lP,mP) - uPPrev(lP,mP));
    uPNext(lP,mP) = uP(lP,mP)*uPlMult + (uP(lP-1,mP) + uP(lP+1,mP) + uP(lP,mP+1) + uP(lP,mP-1))*uPl1Mult + ...
        (uP(lP-1,mP-1) + uP(lP+1,mP-1) + uP(lP-1,mP+1) + uP(lP+1,mP+1))*uPl1dMult + ...
        (uP(lP-2,mP) + uP(lP+2,mP) + uP(lP,mP-2) + uP(lP,mP+2))*uPl2Mult + ...
        uPPrev(lP,mP)*uPPrevlMult + (uPPrev(lP+1,mP) + uPPrev(lP-1,mP) + ...
        uPPrev(lP,mP+1) + uPPrev(lP,mP-1))*uPPrevl1Mult;
    % Update Plate function at the localizer points with Fbp (Force gain from the bridge)
    uPNext(lPcl,mPcl) = uPNext(lPcl,mPcl) + Fbpl/(rhoP*HP*hP^2);
    uPNext(lPcr,mPcr) = uPNext(lPcr,mPcr) + Fbpr/(rhoP*HP*hP^2);
    %% Output
    out(n) = uPNext(floor(Nx/2),floor(Ny/2)) + uBNext(floor(NB/2)) + uS1Next(outPos) + uS2Next(outPos) +uS3Next(outPos);  % plate + mass spring + string
    
    %% Update the state variables
    uS1Prev  = uS1;
    uS1 = uS1Next;    
    uS2Prev  = uS2;
    uS2 = uS2Next;
    uS3Prev  = uS3;
    uS3 = uS3Next;  
    uPPrev  = uP;
    uP = uPNext;
    uBPrev  = uB;
    uB = uBNext;
%     uMPrev  = uM;
%     uM = uMNext;

end
toc
plot(out);