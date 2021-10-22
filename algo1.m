car = 2;
cbr = 2;
cer = 2;

cab = 3;
cae = 3;
cbe = 3;

irx = 0;
iry = 0;
irz = -10;
pw = 30;


P1max = dBtoW(pw);
P2max = dBtoW(pw);
num_iter = 5;


rice_factor_r = 0;
rice_factor_d = 8;

thresh1 = 2;
thresh2 = 0.1;

rounds  = 15;
init1    = 1;
init2    = 1;
freq = 750*1000000;  

sig1 =sqrt(dBtoW(-80));
sig2 =sqrt(dBtoW(-80));
sige =sqrt(dBtoW(-80));

sigl1 =sqrt(dBtoW(-75));
sigl2 =sqrt(dBtoW(-75));

fileID = fopen('vals1.txt','w');
%fprintf(fileID,'Case 1\n\n');
%%%%%%%%%%%%%%%%%%%%%%%1
ax = -30;
ay = 0;
az = 0;
bx = 40;
by = 0;
bz = 0;
ex = 35;
ey = 0;
ez = 0;
L0 =-30;
a = 1;
b = 10;
inits = 1;
%RandomInit(fileID,P1max,P2max ,num_iter,rounds,init1,init2,sig1,sig2,sige,L0,thresh1,thresh2)
%LAfromfilesPower(fileID,P1max,P2max ,num_iter,rounds,init1,init2,sig1,sig2,sige,L0,thresh1,thresh2)
LAfromfiles(fileID,P1max,P2max ,num_iter,rounds,sig1,sig2,sige,sigl1,sigl2,thresh1,a,b);
%LAfromfilesRandomPower(inits,fileID,P1max,P2max ,num_iter,rounds,sig1,sig2,sige,sigl1,sigl2,thresh1,a,b);
