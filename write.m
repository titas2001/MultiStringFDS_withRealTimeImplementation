function write(shamisenString, out,Nx,Ny,note,i)
M = max(abs(out));
outN = out/M;
filename = strcat("C:\University\SMC\SMC10\multiString\MultiStringFDS_withRealTimeImplementation\audio\samples\",string(i),"#","string",string(shamisenString),note,"_Grid",string(Nx),"x", string(Ny),".wav");
audiowrite(filename, outN, 44100);
end

