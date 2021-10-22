%channel initialization
function [] = LAfromfilesPower(fileID,P1max,P2max ,num_iter,rounds,init1,init2,sig1,sig2,sige,L0,thresh1,thresh2)

%from Alice  to IRS always raleigh fading
%from IRS to Eve and Bob try both Los and NLos
%path loss at 1m = -30dB
%path loss exponents alice -> IRS = 2 IRS -> Bob = 2 IRS -> Eve
%raleigh fading variance 1/2, 1/2
%noise variance -100dBm  10 to 50dBm
%thresh IRS = 0.5
%thresh Power = 10
%convergence 1e-3
%L try from 5 to 50 multiples of 5
%IRS 5*n and n(1 to 10) always lies in a plane perpendicular to z axis
    %center of IRS - (100,5,-5)
%signal 750MHz lamda = 1.33333nm
%spacing between IRS element u=3lamda/8
    %IRS height h=5u
    %IRS length w=nu
    %top corner of IRS = (100-(w/2),5-(h/2),-5) =(a,b)
    %center of RIS element (i,j) = (a+(u/2)+(i-1)u,b+(u/2)+(j-1)u,-5)
%alice fixed position.(0,30,30) move along (0,x,x)
%bob . fixed position (110,0,0)
%eve three places (100,0,0), (-100,0,0),(150,0,0). 
%rician distribution s^2 = 1 sigma^2 = 1
%max iter = 40


direct_channels = zeros(num_iter,3);
RIS_channels    = zeros(num_iter,30,3);

for ui=1:num_iter 
    filex = fopen(strcat('C:\Users\Mevan\Desktop\Uni\Research Project\Scripts\Case2\channel',' ',num2str(ui-1),' ','.txt'),'r');
    num = fscanf(filex,'%e');
    h_ab = num(1) + (1j*num(2));
    h_ae = num(3) + (1j*num(4));
    h_be = num(5) + (1j*num(6));
    for j = 1:30
        RIS_channels(ui,j,1) = num(2*j+5)+(1j*num(2*j+6));
    end
    for j = 1:30
        RIS_channels(ui,j,2) = num(2*j+65)+(1j*num(2*j+66));
    end
    for j = 1:30
        RIS_channels(ui,j,3) = num(2*j+125)+(1j*num(2*j+126));
    end
    direct_channels(ui,1) = h_ab;
    direct_channels(ui,2) = h_ae;
    direct_channels(ui,3) = h_be;
    %[s,s1,s2,optP1,optP2] = LAMNoIRS(rounds,init2,P1max,P2max,sig1,sig2,sige,h_ab,h_ae,h_be,thresh2);
    %sum = sum+s;
    %sum1= sum1+s1;
    %sum2= sum2+s2;
    %fprintf(fileID,'NoIRS, s:%.6f, s1:%.6f, s2:%.6f, P1:%.6f, P2:%.6f\n',s,s1,s2,optP1,optP2);
    
end
%
%fprintf(fileID,'NoIRS, s:%.6f, s1:%.6f, s2:%.6f\n\n',sum,sum1,sum2);

for i=2:7
sum = 0;
sum1= 0;
sum2= 0;
for ui=1:num_iter

L = 20;

h = reshape(RIS_channels(ui,1:L,1),[L,1]);
g = reshape(RIS_channels(ui,1:L,2),[L,1]);
e = reshape(RIS_channels(ui,1:L,3),[L,1]);

h_ab = direct_channels(ui,1);
h_ae = direct_channels(ui,2);
h_be = direct_channels(ui,3);

%user h recieved P2
H1i = diag(h)*g;
H1 = [H1i;h_ab];
%user g recieved P1
G1i = diag(g)*h;
G1 = [G1i;h_ab];
%user e recieved
E11i = diag(e)*g;
E11 = [E11i;h_be];
E22i = diag(e)*h;
E22 = [E22i;h_ae];

H = H1*H1';
G = G1*G1';
E1 = E11*E11';
E2 = E22*E22';


iter = 40;
pw = 5*i;
P1max = dBtoW(pw);
P2max = dBtoW(pw);



%[val1,iter1] = tMet(L,P1,P2,P1max,P2max,sig1,sig2,sige,H,G,E1,E2,iter,t)
%fprintf(fileID,'tMet L:%d, val:%.5f, iter:%d\n',L,val1,iter1);

[s,s1,s2,optP1,optP2] = LAM(rounds,init1,L,P1max,P2max,sig1,sig2,sige,H,G,E1,E2,iter,thresh1,thresh2);
sum = sum+s;
sum1= sum1+s1;
sum2= sum2+s2;
fprintf(fileID,'L:%d, iter:%d, s:%.6f, P1:%.6f, P2:%.6f\n',L,ui,s,optP1,optP2);

end
sum1 = sum1/num_iter;
sum2 = sum2/num_iter;
sum  = sum/num_iter;
fprintf(fileID,'L:%d, s:%.6f\n\n',L,sum);
end


end