clear all;
close all;
clc
scale = 100;
shamisenString = 3;
fs = 44100;             % sampling freq
k = 1/fs;               % time step
TS1 = 14.15*9.8;          % applied string tension https://mk0larsenstringsti68.kinstacdn.com/wp-content/uploads/2018/12/Larsen-String-Tension-Charts-18.pdf
TS2 = 14.85*9.8;          % applied string tension
TS3 = 14.36*9.8;          % applied string tension
TP = 4000;            % applied plate tension
rhoS = 1156.48151991993;% material density of the string                        
                        % "Handbook of Fiber Chemistry", Menachem Lewin, Editor, 2nd ed.,1998, Marcel Dekker, pp. 438–441, ISBN 0-8247-9471-0
                        % "ENGINEERING PROPERTIES OF SPIDER SILK"  http://web.mit.edu/course/3/3.064/www/slides/Ko_spider_silk.pdf    
rhoP = 1150;            % nylon https://www.engineeringtoolbox.com/engineering-materials-properties-d_1225.html
HP = 0.0002;             % plate thickness         
EP = 3e+9;              % nylon https://www.engineeringtoolbox.com/engineering-materials-properties-d_1225.html
nu = 0.4;               % Poisson’s ratio nu < 0.5
r = 1.3;                % grid aspect ratio
Lx = r*0.4;             % length of plate in x direction
Ly = (1/r)*0.4;         % length of plate in y direction
LS = 1;               % lenght of the string
LB = 1;                 % lenght of the bridge
durration = 3;          % synthesised sound lenght in seconds
dur = fs*durration;     % synthesised sound lenght in samples
lossS = [100, 10; 1000, 8]; % loss [freq.(Hz), T60(s), freq.(Hz), T60(s)]
lossP = [100, 10; 1000, 8]; % loss [freq.(Hz), T60(s), freq.(Hz), T60(s)]
vP0 = -20;            % initial velocity of a plate
rhoB = 800;             % material density
AreaS1 = 0.00000054129; % string cross sectional area  
AreaS2 = 2.52473376e-7; % string cross sectional area
AreaS3 = 1.398451e-7; % string cross sectional area
AreaB = 2.02e-4;        % bridge cross sectional area
% EB = 9.5e+9;            % Young's modulus dried Red Alder https://amesweb.info/Materials/Youngs-Modulus-of-Wood.aspx
EB = 3.2e+9;            % Young's modulus acrylic https://www.engineeringtoolbox.com/young-modulus-d_417.html
HB = 0.075;             % thickness
rS1 = sqrt(AreaS1/pi);  % string1 radius
rS2 = sqrt(AreaS2/pi);  % string2 radius
rS3 = sqrt(AreaS3/pi);  % string2 radius
ES =  9.9e+9;            % Young modulus "ENGINEERING PROPERTIES OF SPIDER SILK"  http://web.mit.edu/course/3/3.064/www/slides/Ko_spider_silk.pdf    

% gammaS1 = sqrt(TS1/(rhoS*AreaS1));        % String tension
% gammaS2 = sqrt(TS2/(rhoS*AreaS2));        % String tension
% gammaS3 = sqrt(TS3/(rhoS*AreaS3));        % String tension
gammaS1 = sqrt(TS1/(rhoS*AreaS1*LS^2));        % String tension
gammaS2 = sqrt(TS2/(rhoS*AreaS2*LS^2));        % String tension
gammaS3 = sqrt(TS3/(rhoS*AreaS3*LS^2));        % String tension
gammaP = sqrt(TP/(rhoP*HP*Lx*Ly)); 

% cS = sqrt(T0/(rhoS*AreaS));     % not used
% cP = sqrt(T/(rhoP*HP));         % not used


IS1 = (pi*rS1^4)/4;     % string1 inertia
IS2 = (pi*rS2^4)/4;     % string1 inertia
IS3 = (pi*rS3^4)/4;     % string1 inertia
kappaB=sqrt((EB*HB^2)/(12*rhoB*LB^4)); % eq. 7.70 pg. 210

%
kappaS1=sqrt((ES*IS1)/(rhoS*AreaS1*LS^4));
kappaS2=sqrt((ES*IS2)/(rhoS*AreaS2*LS^4));
kappaS3=sqrt((ES*IS3)/(rhoS*AreaS3*LS^4));
D = EP*HP^3 / (12 * (1 - nu^2)); % plate flexural rigidity pg.331
kappaPsq = D / (rhoP * HP * Lx^2 * Ly^2); % pg.332 eq.12.3 kappa^2


% zeta0,1  pg.189 Practical setting for decay times
% sigma0,1 pg.189 eq.7.29
% set scheme for loss parameters for String
% zetaS1 = (-gammaS^2+sqrt(gammaS^4+4*kappaS^2*(2*pi*lossS(1,1))^2))/(2*kappaS^2);
% zetaS2 = (-gammaS^2+sqrt(gammaS^4+4*kappaS^2*(2*pi*lossS(2,1))^2))/(2*kappaS^2);
sigmaS0 = 1.378027748373650;
% 6*log(10)*(-zetaS2/lossS(1,2)+zetaS1/lossS(2,2))/(zetaS1-zetaS2);
sigmaS1 = 3.570213734102943e-03;
% 6*log(10)*(1/lossS(1,2)-1/lossS(2,2))/(zetaS1-zetaS2);
% set scheme for loss parameters for Plate
% zetaP1 = (-gammaP^2+sqrt(gammaP^4+4*kappaPsq*(2*pi*lossP(1,1))^2))/(2*kappaPsq);
% zetaP2 = (-gammaP^2+sqrt(gammaP^4+4*kappaPsq*(2*pi*lossP(2,1))^2))/(2*kappaPsq);
sigmaP0 = 1.378062296963499;
% 6*log(10)*(-zetaP2/lossP(1,2)+zetaP1/lossP(2,2))/(zetaP1-zetaP2);
sigmaP1 = 0.096055930949692;
% 6*log(10)*(1/lossP(1,2)-1/lossP(2,2))/(zetaP1-zetaP2);
% loss parameters for Bar
sigmaB0 =  1.343;
sigmaB1 =  0.00459;

hB = sqrt((4*sigmaB1*k+sqrt((4*sigmaB1*k)^2+16*kappaB^2*k^2))/2); % same method as for the string in pg. 176 but without gamma term 
hS1 = sqrt((gammaS1^2 * k^2 + 4*sigmaS1*k + sqrt((gammaS1^2 * k^2 + 4*sigmaS1*k)^2 + 16*(kappaS1^2)*k^2))/2); % set grid spacing for String eq.7.24-25 pg.176
hS2 = sqrt((gammaS2^2 * k^2 + 4*sigmaS1*k + sqrt((gammaS2^2 * k^2 + 4*sigmaS1*k)^2 + 16*(kappaS2^2)*k^2))/2); % set grid spacing for String eq.7.24-25 pg.176
hS3 = sqrt((gammaS3^2 * k^2 + 4*sigmaS1*k + sqrt((gammaS3^2 * k^2 + 4*sigmaS1*k)^2 + 16*(kappaS3^2)*k^2))/2); % set grid spacing for String eq.7.24-25 pg.176
hP = sqrt((gammaP^2 * k^2 + 4*sigmaP1*k + sqrt((gammaP^2 * k^2 + 4*sigmaP1*k)^2 + 16*kappaPsq*k^2))); % set grid spacing for Plate tromba marina paper eq. 20



NS1 = floor(1/hS1);            % string spatial subdivisions
NS2 = floor(1/hS2);            % string spatial subdivisions
NS3 = floor(1/hS3);            % string spatial subdivisions
NB = floor(1/hB);              % bar spatial subdivisions
Nx = floor(sqrt(r)/hP);        % number of x-subdivisions of spatial domain
Ny = floor(1/(sqrt(r)*hP));    % number of y-subdivisions of spatial domain

if Nx > scale
    Nx = scale;
end


hP = sqrt(r)/min(Nx, Ny);      % reset grid spacing for Plate
hS1 = 1/NS1;                   % reset grid spacing for String1
hS2 = 1/NS2;                   % reset grid spacing for String2
hS3 = 1/NS3;                   % reset grid spacing for String3
hB = 1/NB;                     % reset grid spacing for Bar
Ny = floor(1/(sqrt(r)*hP));    % number of y-subdivisions of spatial domain



%% Connection points
lBc1 = floor((NB/2 - 1)/2);       % bar connection to the 1st string
lBc2 = floor(NB/2);       % bar connection to the 2nd string
lBc3 = ceil((NB/2 - 1)/2)+lBc2;      % bar connection to the 3rd string
lBcl = 1;       % bar left side connection to the plate
lBcr = NB;      % bar right side connection to the plate

lS1c = floor((2*NS1)/(pi*7)); % 1st string connection to the bar
lS2c = floor((2*NS2)/(pi*7)); % 2nd string connection to the bar
lS3c = floor((2*NS3)/(pi*7)); % 3rd string connection to the bar

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
% frettingpos = [0,4, 8,13,17,21,24,27,30,33,36,39,41,43,45,47,49,53,55,56,58,59,60,61,63,64,65];   %1st string
% frettingpos = [0,3,6,9,12,14,16,19,21,23,25,26,28,30,31,32,34,35,36,37,38,39,40,41,42,43];        %2nd string
frettingpos = [0,2,5,7,9,11,12,14,16,17,18,20,21,22,23,24,25,26,27,28,29];                        %3rd string
% noteName = ["C4","Db4","D4","Eb4","E4","F4","Gb4","G4","Ab4","A4","Bb4","B4","C5","Db5","D5"...
%     ,"Eb5","E5","F5","Db6","D6","Eb6","E6","F5","Gb5","G5","Ab5","A5","Bb5","B5","C6","Db6","D6","Eb6"]; %1st string
% noteName = ["G4","Ab4","A4","Bb4","B4","C5","Db5","D5","Eb5","E5","F5","Gb5","G5","Ab5","A5"...          %2nd string
%     ,"Bb5","B5","C6","Db6","D6","Eb6","E6","F6","Gb6","G6","Ab6"];
noteName = ["C5","Db5","D5","Eb5","E5","F5","Gb5","G5","Ab5","A5","Bb5","B5","C6","Db6","D6"...          %3rd string
,"Eb6","E6","F6","Gb6","G6","Ab6"];

for i=1:length(frettingpos)
%% Intialise states of the system

% Strings

uS1Next = zeros(NS1,1);
uS1 = zeros(NS1,1);
if shamisenString == 1
    width = round(NS1/10);
    excitationRange = 1:width;
    uS1(excitationRange + floor((NS1*5)/(pi*6))) = hann(width);
end
uS1Prev = uS1;
uS2Next = zeros(NS2,1);
uS2 = zeros(NS2,1);
if shamisenString == 2
    width = round(NS2/10);
    excitationRange = 1:width;
    uS2(excitationRange + floor((NS2*5)/(pi*6))) = hann(width);
end
uS2Prev = uS2;
uS3Next = zeros(NS3,1);
uS3 = zeros(NS3,1);
if shamisenString == 3
    width = round(NS3/10);
    excitationRange = 1:width;
    uS3(excitationRange + floor((NS3*5)/(pi*6))) = hann(width);
end
uS3Prev = uS3;


% Plate
uPNext = zeros(Nx,Ny);
uP = zeros(Nx,Ny);
uPPrev = zeros(Nx,Ny);
% uP(ceil(Nx/2-Nx/8):floor(Nx/2+Nx/8),ceil(Ny/2-Ny/8):floor(Ny/2+Ny/8)) = ... here plate is excited by velocity
%     k*vP0*hamming_3d(length(ceil(Nx/2-Nx/8):floor(Nx/2+Nx/8)),length(ceil(Ny/2-Ny/8):floor(Ny/2+Ny/8)),1);
% outPosS1 = floor((NS1*5)/(pi*4));
% outPosS2 = floor((NS2*5)/(pi*4));
% outPosS3 = floor((NS3*5)/(pi*4));
outPosP = [floor(2*Nx/(pi)) floor(Ny/(pi))];
outPosB = floor(2*NB/pi);
outPosS1 =floor((2*NS1)/(pi*7))+4;
outPosS2 =floor((2*NS2)/(pi*7))+4;
outPosS3 =floor((2*NS3)/(pi*7))+4;
% Bar
uBNext = zeros(NB,1);
uB = zeros(NB,1);
uBPrev = zeros(NB,1);

% Output
out = zeros(dur,1);
%% Intialise l for update equations
lP = 3:Nx-2;
mP = 3:Ny-2;
lS1 = 3+(0):NS1-2-(0);
lS2 = 3+(0):NS2-2-(frettingpos(i));
lS3 = 3+(0):NS3-2-(0);
lB = 3:NB-2;

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
    out(n) = uPNext(outPosP(1),outPosP(2)) + ...                   % plate
        uBNext(outPosB) + ...                                      % bridge
        uS1Next(outPosS1) + uS2Next(outPosS2) +uS3Next(outPosS2);  % string
    
    %% Update the state variables
    uS1Prev  = uS1;
    uS1 = uS1Next;    
%     uS2Prev  = uS2;
%     uS2 = uS2Next;
%     uS3Prev  = uS3;
%     uS3 = uS3Next;  
    uPPrev  = uP;
    uP = uPNext;
    uBPrev  = uB;
    uB = uBNext;
%      if n == floor(dur/3)
%         uexS2 = zeros(NS2,1);
%         uexS2((1:floor(NS2/10)) + floor(NS2/5)) = hann(floor(NS2/10));
%         uS2((1:floor(NS2/10)) + floor(NS2/5)) = uS2((1:floor(NS2/10)) + floor(NS2/5)) + uexS2((1:floor(NS2/10)) + floor(NS2/5));
% %         uP(ceil(Nx/2-Nx/8):floor(Nx/2+Nx/8),ceil(Ny/2-Ny/8):floor(Ny/2+Ny/8)) = ...
% %             uP(ceil(Nx/2-Nx/8):floor(Nx/2+Nx/8),ceil(Ny/2-Ny/8):floor(Ny/2+Ny/8)) + ... here plate is excited by velocity
% %             k*vP0*hamming_3d(length(ceil(Nx/2-Nx/8):floor(Nx/2+Nx/8)),length(ceil(Ny/2-Ny/8):floor(Ny/2+Ny/8)),1);
%      end
%      if n ==floor(dur*2/3)
%         uexS3 = zeros(NS3,1);
%         uexS3((1:floor(NS3/10)) + floor(NS3/5)) = hann(floor(NS3/10));
%         uS3((1:floor(NS3/10)) + floor(NS3/5)) = uS3((1:floor(NS3/10)) + floor(NS3/5)) + uexS3((1:floor(NS3/10)) + floor(NS3/5));
%      end
    uS2Prev  = uS2;
    uS2 = uS2Next;
    uS3Prev  = uS3;
    uS3 = uS3Next;  

end
toc
write(shamisenString,out,Nx,Ny,noteName(i));
end
plot(out);