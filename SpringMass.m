clear all;
close all;
clc
tic
fs = 44100;
k = 1/fs;
duration = 1;          % synthesised sound duration in s
dur = fs*duration;
f0 = 100;
R = 3;



%initialize u states
uNext = 0;
u = zeros(100,1);
u(1) = 1;
uPrev = u;


l = 1;
for n = 1:dur
    uNext(l) = -4*pi^2*f0^2*k^2*u(1) - R*k*(u(l) - uPrev(l)) + 2*u(l) - uPrev(l);
    
    
    out(n) = uNext(l);
%     plot(uNext);
%     ylim([-1,1]);
%     drawnow;
    uPrev  = u;
    u = uNext;
end
plot(out);
