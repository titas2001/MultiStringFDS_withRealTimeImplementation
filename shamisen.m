clear all;
close all;
clc

fs = 44100;             % sampling freq
k = 1/fs;               % time step
T0 = 60;                % applied string tension
T = 40;                 % applied plate tension
rhoS = 7700;            % material density of the string
rhoP = 1150;            % nylon https://www.engineeringtoolbox.com/engineering-materials-properties-d_1225.html
HP = 0.002;             % plate thickness         
EP = 3e+9;               % nylon https://www.engineeringtoolbox.com/engineering-materials-properties-d_1225.html
nu = 0.4;               % Poisson�s ratio nu < 0.5
r = 1.3;                % grid aspect ratio
Lx = r*0.4;             % length of plate in x direction
Ly = (1/r)*0.4;         % length of plate in y direction
LS = 0.5;               % lenght of the string
LB = 1;                 % lenght of the bridge
durration = 1;          % synthesised sound lenght in seconds
dur = fs*durration;     % synthesised sound lenght in samples
lossS = [100, 10; 1000, 8]; % loss [freq.(Hz), T60(s), freq.(Hz), T60(s)]
lossP = [100, 10; 1000, 8]; % loss [freq.(Hz), T60(s), freq.(Hz), T60(s)]
% vP0 = -10;            % initial velocity of a plate
rhoB = 800;             % material density
AreaS1 = 2.02682992e-7; % string cross sectional area
AreaS2 = 1.52682992e-7; % string cross sectional area
AreaS3 = 0.52682992e-7; % string cross sectional area
AreaB = 2.02e-4;        % bridge cross sectional area
EB = 9.5e+9;            % Young's modulus dried Red Alder https://amesweb.info/Materials/Youngs-Modulus-of-Wood.aspx
HB = 0.075;             % thickness
rS1 = sqrt(AreaS1/pi);  % string1 radius
rS2 = sqrt(AreaS2/pi);  % string2 radius
rS3 = sqrt(AreaS3/pi);  % string2 radius
ES = 3e+9;              % Young's Modulus of nylon

gammaS1 = sqrt(T0/(rhoS*AreaS1*LS^2));        % String tension
gammaS2 = sqrt(T0/(rhoS*AreaS2*LS^2));        % String tension
gammaS3 = sqrt(T0/(rhoS*AreaS3*LS^2));        % String tension
gammaP = sqrt(T/(rhoP*HP*Lx*Ly)); 

% cS = sqrt(T0/(rhoS*AreaS));     % not used
% cP = sqrt(T/(rhoP*HP));         % not used


IS1 = (pi*rS1^4)/4;     % string1 inertia
IS2 = (pi*rS2^4)/4;     % string1 inertia
IS3 = (pi*rS3^4)/4;     % string1 inertia
kappaB=sqrt((EB*HB^2)/(12*rhoB*LB^4)); % eq. 7.70 pg. 210

%
kappaS1=(ES*IS1)/rhoS;
kappaS2=(ES*IS2)/rhoS;
kappaS3=(ES*IS3)/rhoS;
D = EP*HP^3 / (12 * (1 - nu^2)); % plate flexural rigidity pg.341
kappaPsq = D / (rhoP * HP * Lx^2 * Ly^2); % pg.342 eq.12.3 kappa^2


% zeta0,1  pg.189 Practical setting for decay times
% sigma0,1 pg.189 eq.7.29
% set scheme for loss parameters for String
% zetaS1 = (-gammaS^2+sqrt(gammaS^4+4*kappaS^2*(2*pi*lossS(1,1))^2))/(2*kappaS^2);
% zetaS2 = (-gammaS^2+sqrt(gammaS^4+4*kappaS^2*(2*pi*lossS(2,1))^2))/(2*kappaS^2);
sigmaS0 = 1.378027748373650;
% 6*log(10)*(-zetaS2/lossS(1,2)+zetaS1/lossS(2,2))/(zetaS1-zetaS2);
sigmaS1 = 3.570213734102943e-04;
% 6*log(10)*(1/lossS(1,2)-1/lossS(2,2))/(zetaS1-zetaS2);
% set scheme for loss parameters for Plate
zetaP1 = (-gammaP^2+sqrt(gammaP^4+4*kappaPsq*(2*pi*lossP(1,1))^2))/(2*kappaPsq);
zetaP2 = (-gammaP^2+sqrt(gammaP^4+4*kappaPsq*(2*pi*lossP(2,1))^2))/(2*kappaPsq);
sigmaP0 = 1.378062296963499;
% 6*log(10)*(-zetaP2/lossP(1,2)+zetaP1/lossP(2,2))/(zetaP1-zetaP2);
sigmaP1 = 0.096055930949692;
% 6*log(10)*(1/lossP(1,2)-1/lossP(2,2))/(zetaP1-zetaP2);
% loss parameters for Bar
sigmaB0 =  1.343;
sigmaB1 =  0.00459;

hB = sqrt((4*sigmaB1*k+sqrt((4*sigmaB1*k)^2+16*kappaB^2*k^2))/2);
hS1 = sqrt((gammaS1^2 * k^2 + 4*sigmaS1*k + sqrt((gammaS1^2 * k^2 + 4*sigmaS1*k)^2 + 16*(kappaS1^2)*k^2))/2); % set grid spacing for String eq.7.26 pg.188
hS2 = sqrt((gammaS2^2 * k^2 + 4*sigmaS1*k + sqrt((gammaS2^2 * k^2 + 4*sigmaS1*k)^2 + 16*(kappaS2^2)*k^2))/2); % set grid spacing for String eq.7.26 pg.188
hS3 = sqrt((gammaS3^2 * k^2 + 4*sigmaS1*k + sqrt((gammaS3^2 * k^2 + 4*sigmaS1*k)^2 + 16*(kappaS3^2)*k^2))/2); % set grid spacing for String eq.7.26 pg.188
hP = sqrt((gammaP^2 * k^2 + 4*sigmaP1*k + sqrt((gammaP^2 * k^2 + 4*sigmaP1*k)^2 + 16*kappaPsq*k^2))); % set grid spacing for Plate tromba marina paper eq. 20

NS1 = floor(1/hS1);             % string spatial subdivisions
NS2 = floor(1/hS2);             % string spatial subdivisions
NS3 = floor(1/hS3);             % string spatial subdivisions
NB = floor(1/hB);              % bar spatial subdivisions
Nx = floor(sqrt(r)/hP);        % number of x-subdivisions of spatial domain
Ny = floor(1/(sqrt(r)*hP));    % number of y-subdivisions of spatial domain

if Nx > 30
    Nx = 30;
end


hP = sqrt(r)/min(Nx, Ny);      % reset grid spacing for Plate
hS1 = 1/NS1;                     % reset grid spacing for String1
hS2 = 1/NS2;                     % reset grid spacing for String2
hS3 = 1/NS3;                     % reset grid spacing for String3
hB = 1/NB;                     % reset grid spacing for Bar
Ny = floor(1/(sqrt(r)*hP));    % number of y-subdivisions of spatial domain
%% Intialise states of the system

% Strings
uS1Next = zeros(NS1,1);
uS1 = zeros(NS1,1);
width = floor(NS1/10);
excitationRange = 1:width;
uS1(excitationRange + floor(NS1/5)) = hann(width);
uS1Prev = uS1;
uS2Next = zeros(NS2,1);
uS2 = zeros(NS2,1);
% width = floor(NS2/10);
% excitationRange = 1:width;
% uS2(excitationRange + floor(NS2/5)) = hann(width);
uS2Prev = uS2;
uS3Next = zeros(NS3,1);
uS3 = zeros(NS3,1);
% width = floor(NS3/10);
% excitationRange = 1:width;
% uS3(excitationRange + floor(NS3/5)) = hann(width);
uS3Prev = uS3;


% Plate
uPNext = zeros(Nx,Ny);
uP = zeros(Nx,Ny);
uPPrev = zeros(Nx,Ny);
% uP(ceil(Nx/2-Nx/8):floor(Nx/2+Nx/8),ceil(Ny/2-Ny/8):floor(Ny/2+Ny/8)) = ... here plate is excited by velocity
%     k*vP0*hamming_3d(length(ceil(Nx/2-Nx/8):floor(Nx/2+Nx/8)),length(ceil(Ny/2-Ny/8):floor(Ny/2+Ny/8)),1);
outPosS1 = floor(NS1/pi);
outPosS2 = floor(NS2/pi);
outPosS3 = floor(NS3/pi);
outPosP = [floor(2*Nx/(pi)) floor(Ny/(pi))];
outPosB = floor(2*NB/pi);
% Bar
uBNext = zeros(NB,1);
uB = zeros(NB,1);
uBPrev = zeros(NB,1);

% Output
out = zeros(dur,1);
%% Intialise l for update equations
lP = 3:Nx-2;
mP = 3:Ny-2;
lS1 = 3:NS1-2;
lS2 = 3:NS2-2;
lS3 = 3:NS3-2;
lB = 3:NB-2;

%% Connection points
lBc1 = 5;       % bar connection to the 1st string
lBc2 = 9;       % bar connection to the 2nd string
lBc3 = 13;      % bar connection to the 3rd string
lBcl = 1;       % bar left side connection to the plate
lBcr = 17;      % bar right side connection to the plate

lM = 1;         % will be deprecated
lMc = 1;        % will be deprecated

lS1c = NS1 - floor(NS1/8); % 1st string connection to the bar
lS2c = NS2 - floor(NS2/8); % 2nd string connection to the bar
lS3c = NS3 - floor(NS3/8); % 3rd string connection to the bar

lPc = Nx - floor(Nx/3); % will be deprecated
mPc = Ny - floor(Ny/3); % will be deprecated

lPcl = Nx - floor(2*Nx/5); % Plate connection to the bar on the left side x coordinate
lPcr = Nx - floor(3*Nx/5); % Plate connection to the bar on the right side x coordinate
mPcl = Ny - floor(Ny/4);   % Plate connection to the bar on the left side y coordinate
mPcr = Ny - floor(Ny/4);   % Plate connection to the bar on the right side y coordinate

%% Multipliers

% Strings
uS1lMult = (((-2*gammaS1^2)/hS1^2 - 6*kappaS1^2/hS1^4 - 4*sigmaS1/(k*hS1^2))*k^2 + 2)/(k*sigmaS0 + 1);
uS1l1Mult = (gammaS1^2/hS1^2 + 4*kappaS1^2/hS1^4 + 2*sigmaS1/(k*hS1^2))*k^2/(k*sigmaS0 + 1);
uS1l2Mult = ((-1*k^2*kappaS1^2)/(hS1^4))/(k*sigmaS0 + 1);
uS1PrevlMult = ((4*sigmaS1*k^2)/(k*hS1^2) + k*sigmaS0 - 1)/(k*sigmaS0 + 1);
uS1Prevl1Mult = ((-2*sigmaS1*k^2)/(k*hS1^2))/(k*sigmaS0 + 1);

uS2lMult = (((-2*gammaS2^2)/hS2^2 - 6*kappaS2^2/hS2^4 - 4*sigmaS1/(k*hS2^2))*k^2 + 2)/(k*sigmaS0 + 1);
uS2l1Mult = (gammaS2^2/hS2^2 + 4*kappaS2^2/hS2^4 + 2*sigmaS1/(k*hS2^2))*k^2/(k*sigmaS0 + 1);
uS2l2Mult = ((-1*k^2*kappaS2^2)/(hS2^4))/(k*sigmaS0 + 1);
uS2PrevlMult = ((4*sigmaS1*k^2)/(k*hS2^2) + k*sigmaS0 - 1)/(k*sigmaS0 + 1);
uS2Prevl1Mult = ((-2*sigmaS1*k^2)/(k*hS2^2))/(k*sigmaS0 + 1);

uS3lMult = (((-2*gammaS3^2)/hS3^2 - 6*kappaS3^2/hS3^4 - 4*sigmaS1/(k*hS3^2))*k^2 + 2)/(k*sigmaS0 + 1);
uS3l1Mult = (gammaS3^2/hS3^2 + 4*kappaS3^2/hS3^4 + 2*sigmaS1/(k*hS3^2))*k^2/(k*sigmaS0 + 1);
uS3l2Mult = ((-1*k^2*kappaS3^2)/(hS3^4))/(k*sigmaS0 + 1);
uS3PrevlMult = ((4*sigmaS1*k^2)/(k*hS3^2) + k*sigmaS0 - 1)/(k*sigmaS0 + 1);
uS3Prevl1Mult = ((-2*sigmaS1*k^2)/(k*hS3^2))/(k*sigmaS0 + 1);

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
Fs1bMult = 1/(1/(rhoB*AreaB*hB * (sigmaB0 + 1)) + 1/(rhoS*AreaS1*hS1 * (sigmaS0 + 1)));
Fs2bMult = 1/(1/(rhoB*AreaB*hB * (sigmaB0 + 1)) + 1/(rhoS*AreaS2*hS2 * (sigmaS0 + 1)));
Fs3bMult = 1/(1/(rhoB*AreaB*hB * (sigmaB0 + 1)) + 1/(rhoS*AreaS3*hS3 * (sigmaS0 + 1)));
FbpMult = 1/(-1/(rhoB*AreaB*hB * (sigmaB0 + 1)) - 1/(rhoP*HP*hP^2 * (sigmaP0 + 1)));

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
    uS1Next(lS1) = uS1(lS1)*uS1lMult + (uS1(lS1-1) + uS1(lS1+1))*uS1l1Mult + (uS1(lS1-2) + uS1(lS1+2))*uS1l2Mult + ...
        uS1Prev(lS1)*uS1PrevlMult + (uS1Prev(lS1-1)+ uS1Prev(lS1+1))*uS1Prevl1Mult;
    uS2Next(lS2) = uS2(lS2)*uS2lMult + (uS2(lS2-1) + uS2(lS2+1))*uS2l1Mult + (uS2(lS2-2) + uS2(lS2+2))*uS2l2Mult + ...
        uS2Prev(lS2)*uS2PrevlMult + (uS2Prev(lS2-1)+ uS2Prev(lS2+1))*uS2Prevl1Mult;
    uS3Next(lS3) = uS3(lS3)*uS3lMult + (uS3(lS3-1) + uS3(lS3+1))*uS3l1Mult + (uS3(lS3-2) + uS3(lS3+2))*uS3l2Mult + ...
        uS3Prev(lS3)*uS3PrevlMult + (uS3Prev(lS3-1)+ uS3Prev(lS3+1))*uS3Prevl1Mult;

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
    
%% Calculate the forces
    % Force from the Strings to the bridge
    Fs1b = Fs1bMult * (-uBNext(lBc1) + uS1Next(lS1c));
        
    Fs2b = Fs2bMult * (-uBNext(lBc2) + uS2Next(lS2c));
    
    Fs3b = Fs3bMult * (-uBNext(lBc3) + uS3Next(lS3c));
    
    % Force from bridge' left and right mounting points to the plate
    Fbpl = FbpMult*(-uBNext(lBcl) + uPNext(lPcl, mPcl));

    Fbpr = FbpMult*(-uBNext(lBcr) + uPNext(lPcr, mPcr));
        
    
%% Update equations at localizer points
    % Update Strings equation at the localizer point with Fsm (Force loss to the bridge)
    uS1Next(lS1c) = uS1Next(lS1c) - Fs1b/(rhoS * AreaS1 *hS1 * (sigmaS0 + 1));
    uS2Next(lS2c) = uS2Next(lS2c) - Fs2b/(rhoS * AreaS2 *hS2 * (sigmaS0 + 1));
    uS3Next(lS3c) = uS3Next(lS3c) - Fs3b/(rhoS * AreaS3 *hS3 * (sigmaS0 + 1));
        

    % Update Bar function at the localizer point with Fs1b (Force gain from the string 1)
    uBNext(lBc1) = uBNext(lBc1) + Fs1b/(rhoB * AreaB *hB * (sigmaB0 + 1));
    % Update Bar function at the localizer point with Fs2b (Force gain from the string 2)
    uBNext(lBc2) = uBNext(lBc2) + Fs2b/(rhoB * AreaB *hB * (sigmaB0 + 1));
    % Update Bar function at the localizer point with Fs3b (Force gain from the string 3)
    uBNext(lBc3) = uBNext(lBc3) + Fs3b/(rhoB * AreaB *hB * (sigmaB0 + 1));
    % Update Bar function at the localizer point lBl with Fbp (Force gain from the string 1)
    uBNext(lBcl) = uBNext(lBcl) - Fbpl/(rhoB * AreaB *hB * (sigmaB0 + 1));
    % Update Bar function at the localizer point lBr with Fbp (Force gain from the string 1)
    uBNext(lBcr) = uBNext(lBcr) - Fbpr/(rhoB * AreaB *hB * (sigmaB0 + 1));


    % Update Plate function at the localizer points with Fbp (Force gain from the bridge)
    uPNext(lPcl,mPcl) = uPNext(lPcl,mPcl) + Fbpl/(rhoP*HP*hP^2*(sigmaP0 + 1));
    uPNext(lPcr,mPcr) = uPNext(lPcr,mPcr) + Fbpr/(rhoP*HP*hP^2*(sigmaP0 + 1));
%% plot    
%     variable = uPNext;     
%     subplot(9,1,1);
%     plot(uS1Next);
%     subplot(9,1,2);
%     plot(uS2Next);
%     subplot(9,1,3);
%     plot(uS3Next);
%     subplot(6,1,3);
%     plot(uBNext);
%     subplot(3,1,3);
%     mesh(variable);
% %     imagesc(variable);
%     drawnow;
    
    %% Output
    out(n) = uPNext(outPosP(1),outPosP(2)) + uBNext(outPosB) + uS1Next(outPosS1) + uS2Next(outPosS2) +uS3Next(outPosS2);  % plate + mass spring + string
    
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