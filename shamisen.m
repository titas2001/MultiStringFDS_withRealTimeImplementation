clear all;
close all;
clc

fs = 44100;             % sampling freq
f0 = 150;               % fundamental freq
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
rS = 0.0005;
ES = 3e+9;
IS = (pi*rS^4)/4;
kappaB=sqrt((EB*HB^2)/(12*rhoB*LB^4)); % eq. 7.70 pg. 210

%
kappaS=(ES*IS)/rhoS;
D = E*HP^3 / (12 * (1 - nu^2)); % plate flexural rigidity pg.341
kappaPsq = D / (rhoP * HP * Lx^2 * Ly^2); % pg.342 eq.12.3 kappa^2


% zeta0,1  pg.189 Practical setting for decay times
% sigma0,1 pg.189 eq.7.29
% set scheme for loss parameters for String
zetaS1 = (-gammaS^2+sqrt(gammaS^4+4*kappaS^2*(2*pi*lossS(1,1))^2))/(2*kappaS^2);
zetaS2 = (-gammaS^2+sqrt(gammaS^4+4*kappaS^2*(2*pi*lossS(2,1))^2))/(2*kappaS^2);
sigmaS0 = 0;
% 6*log(10)*(-zetaS2/lossS(1,2)+zetaS1/lossS(2,2))/(zetaS1-zetaS2);
sigmaS1 = 0;
% 6*log(10)*(1/lossS(1,2)-1/lossS(2,2))/(zetaS1-zetaS2);
% set scheme for loss parameters for Plate
zetaP1 = (-gammaP^2+sqrt(gammaP^4+4*kappaPsq*(2*pi*lossP(1,1))^2))/(2*kappaPsq);
zetaP2 = (-gammaP^2+sqrt(gammaP^4+4*kappaPsq*(2*pi*lossP(2,1))^2))/(2*kappaPsq);
sigmaP0 = 6*log(10)*(-zetaP2/lossP(1,2)+zetaP1/lossP(2,2))/(zetaP1-zetaP2);
sigmaP1 = 6*log(10)*(1/lossP(1,2)-1/lossP(2,2))/(zetaP1-zetaP2);
% loss parameters for Bar
sigmaB0 = 1.343;
sigmaB1 = 0.00459;

hB = sqrt((4*sigmaB1*k+sqrt((4*sigmaB1*k)^2+16*kappaB^2*k^2))/2);
hS = sqrt((gammaS^2 * k^2 + 4*sigmaS1*k + sqrt((gammaS^2 * k^2 + 4*sigmaS1*k)^2 + 16*(kappaS^2)*k^2))/2); % set grid spacing for String eq.7.26 pg.188
hP = sqrt((gammaP^2 * k^2 + 4*sigmaP1*k + sqrt((gammaP^2 * k^2 + 4*sigmaP1*k)^2 + 16*kappaPsq*k^2))); % set grid spacing for Plate tromba marina paper eq. 20

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
uSl2Mult = ((-1*k^2*kappaS^2)/(hS^4))/(k*sigmaS0 + 1);
uSPrevlMult = ((4*sigmaS1*k^2)/(k*hS^2) + k*sigmaS0 - 1)/(k*sigmaS0 + 1);
uSPrevl1Mult = ((-2*sigmaS1*k^2)/(k*hS^2))/(k*sigmaS0 + 1);

% Bar
uBlMult = (2*hB^4 - 4*hB^2*k*sigmaB1 - 6*k^2*kappaB^2)/(hB^4*(k*sigmaB0 + 1));
uBl1Mult = (2*hB^2*k*sigmaB1 + 4*k^2*kappaB^2)/(hB^4*(k*sigmaB0 + 1));
uBl2Mult = (-1*k^2*kappaB^2)/(hB^4*(k*sigmaB0 + 1));
uBPrevlMult = (hB^4*k*sigmaB0 - hB^4 + 4*hB^2*k*sigmaB1)/(hB^4*(k*sigmaB0 + 1));
uBPrevl1Mult = (-2*sigmaB1*k)/(hB^2*(k*sigmaB0 + 1));

% Plate
uPlMult = ((-20*kappaPsq/hP^4 - 4*gammaP^2/hP^2 - 8*sigmaP1/(k*hP^2))*k^2 + 2)/(k*sigmaP0 + 1); 
uPl1Mult = (8*kappaPsq/hP^4 + gammaP^2/hP^2 + 2*sigmaP1/(k*hP^2))*k^2/(k*sigmaP0 + 1);
uPl1dMult = (-2*kappaPsq*k^2)/(hP^4)/(k*sigmaP0 + 1);
uPl2Mult = (-1*kappaPsq*k^2)/(hP^4)/(k*sigmaP0 + 1);
uPPrevlMult = ((8*sigmaP1*k^2)/(k*hP^2) + k*sigmaP0 - 1)/(k*sigmaP0 + 1);
uPPrevl1Mult = ((-2*sigmaP1*k^2)/(k*hP^2))/(k*sigmaP0 + 1);

% Forces
FsbMult = 1/(1/(rhoB*AreaB*hB) + 1/(rhoS*AreaS*hS));
FbpMult = 1/(-1/(rhoB*AreaB*hB) - 1/(rhoP*HP*hP^2));

%%
tic
for n = 1:dur
    % Calculate virtual grid points
    uB0 = 2*uB(1)-uB(2);                    % uB(0)
    uBm1 = 2*(uB0-uB(2))+uB(3);             % uB(-1)
    uBPrev0 = 2*uBPrev(1)-uBPrev(2);        % uBPrev(0)
    uBNp1 = 2*uB(NB) - uB(NB-1);            % uB(N+1)
    uBNp2 = 2*(uBNp1 - uB(NB-1)) + uB(NB-2);% uB(N+2)
    uBPrevNp1 = 2*uBPrev(NB) - uBPrev(NB-1);% uBPrev(N+1)
    
    
    
%____UPDATE_EQUATIONS________________________________________________________________________________________
    
    %% Update equation for the Strings
    uS1Next(lS) = uS1(lS)*uSlMult + (uS1(lS-1) + uS1(lS+1))*uSl1Mult + (uS1(lS-2) + uS1(lS+2))*uSl2Mult + ...
        uS1Prev(lS)*uSPrevlMult + (uS1Prev(lS-1)+ uS1Prev(lS+1))*uSPrevl1Mult;
    uS2Next(lS) = uS2(lS)*uSlMult + (uS2(lS-1) + uS2(lS+1))*uSl1Mult + (uS2(lS-2) + uS2(lS+2))*uSl2Mult + ...
        uS2Prev(lS)*uSPrevlMult + (uS2Prev(lS-1)+ uS2Prev(lS+1))*uSPrevl1Mult;
    uS3Next(lS) = uS3(lS)*uSlMult + (uS3(lS-1) + uS3(lS+1))*uSl1Mult + (uS3(lS-2) + uS3(lS+2))*uSl2Mult + ...
        uS3Prev(lS)*uSPrevlMult + (uS3Prev(lS-1)+ uS3Prev(lS+1))*uSPrevl1Mult;

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
 
        
    %% Update equation of the Plate
    uPNext(lP,mP) = (1/(k*sigmaP0 + 1))*(((-(kappaPsq)/hP^4)*((uP(lP+2,mP) + uP(lP-2,mP) + uP(lP,mP+2) + uP(lP,mP-2)) + ...
        2*(uP(lP+1,mP+1) + uP(lP+1,mP-1) + uP(lP-1,mP+1) + uP(lP-1,mP-1)) - ...
        8*(uP(lP+1,mP) + uP(lP-1,mP) + uP(lP,mP+1) + uP(lP,mP-1)) + 20*uP(lP,mP)) + ...
        (gammaP^2/hP^2)*(uP(lP+1,mP) + uP(lP-1,mP) + uP(lP,mP+1) + uP(lP,mP-1) - 4*uP(lP,mP)) + ...
        ((2*sigmaP1)/(k*hP^2))*(uP(lP+1,mP) + uP(lP-1,mP) + uP(lP,mP+1) + uP(lP,mP-1) - 4*uP(lP,mP) - ...
        (uPPrev(lP+1,mP) + uPPrev(lP-1,mP) + uPPrev(lP,mP+1) + uPPrev(lP,mP-1) - 4*uPPrev(lP,mP))))*k^2 + ...
        k*sigmaP0*uPPrev(lP,mP) + 2*uP(lP,mP) - uPPrev(lP,mP));
    uPNext(lP,mP) = uP(lP,mP)*uPlMult + (uP(lP-1,mP) + uP(lP+1,mP) + uP(lP,mP+1) + uP(lP,mP-1))*uPl1Mult + ...
        (uP(lP-1,mP-1) + uP(lP+1,mP-1) + uP(lP-1,mP+1) + uP(lP+1,mP+1))*uPl1dMult + ...
        (uP(lP-2,mP) + uP(lP+2,mP) + uP(lP,mP-2) + uP(lP,mP+2))*uPl2Mult + ...
        uPPrev(lP,mP)*uPPrevlMult + (uPPrev(lP+1,mP) + uPPrev(lP-1,mP) + ...
        uPPrev(lP,mP+1) + uPPrev(lP,mP-1))*uPPrevl1Mult;
    
%% Calculate the forces
    % Force from the Strings to the bridge
    Fs1b = FsbMult * (uBNext(lBc1) + uS1Next(lS1c));
        
    Fs2b = FsbMult * (uBNext(lBc2) + uS2Next(lS2c));
    
    Fs3b = FsbMult * (uBNext(lBc3) + uS3Next(lS3c));
    
    % Force from bridge' left and right mounting points to the plate
    Fbpl = FbpMult*(uBNext(lBcl) + uPNext(lPcl, mPcl));

    Fbpr = FbpMult*(uBNext(lBcr) + uPNext(lPcr, mPcr));
        
    
%% Update equations at localizer points
    % Update Strings equation at the localizer point with Fsm (Force loss to the bridge)
    uS1Next(lS1c) = uS1Next(lS1c) - Fs1b/(rhoS * AreaS *hS);
    uS2Next(lS2c) = uS2Next(lS2c) - Fs2b/(rhoS * AreaS *hS);
    uS3Next(lS3c) = uS3Next(lS3c) - Fs3b/(rhoS * AreaS *hS);
        

    % Update Bar function at the localizer point with Fs1b (Force gain from the string 1)
    uBNext(lBc1) = uBNext(lBc1) + Fs1b/(rhoB * AreaB *hB);
    % Update Bar function at the localizer point with Fs2b (Force gain from the string 2)
    uBNext(lBc2) = uBNext(lBc2) + Fs2b/(rhoB * AreaB *hB);
    % Update Bar function at the localizer point with Fs3b (Force gain from the string 3)
    uBNext(lBc3) = uBNext(lBc3) + Fs3b/(rhoB * AreaB *hB);
    % Update Bar function at the localizer point lBl with Fbp (Force gain from the string 1)
    uBNext(lBcl) = uBNext(lBcl) - Fbpl/(rhoB * AreaB *hB);
    % Update Bar function at the localizer point lBr with Fbp (Force gain from the string 1)
    uBNext(lBcr) = uBNext(lBcr) - Fbpr/(rhoB * AreaB *hB);


    % Update Plate function at the localizer points with Fbp (Force gain from the bridge)
    uPNext(lPcl,mPcl) = uPNext(lPcl,mPcl) + Fbpl/(rhoP*HP*hP^2);
    uPNext(lPcr,mPcr) = uPNext(lPcr,mPcr) + Fbpr/(rhoP*HP*hP^2);
%% plot    
    
    subplot(2,1,1);
    plot(uS1Next);
    subplot(2,1,2);
    plot(uBNext);
    drawnow;
    
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