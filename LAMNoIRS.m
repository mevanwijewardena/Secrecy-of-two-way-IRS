function [s,optP1,optP2] = LAMNoIRS(P1max,P2max,sig1,sig2,sige,sigl1,sigl2,h_ab,h_ae,h_be)
hab = abs(h_ab)^2;
hbe = abs(h_be)^2;
hae = abs(h_ae)^2;

s= 0;
optP1 = 0;
optP2 = 0;
%stationary point
P2 = (((sig2^2+sigl2^2)*hae)-(sige^2*hab))/(hab*hbe);
P1 = (((sig1^2+sigl1^2)*hbe)-(sige^2*hab))/(hab*hae);
if((P1>=0) && (P2>=0) && (secrecy1(P1,P2,hab,hae,hbe,sig1,sig2,sige,sigl1,sigl2)>0))
    s = secrecy1(P1,P2,hab,hae,hbe,sig1,sig2,sige,sigl1,sigl2);
    optP1 = P1;
    optP2 = P2;
end


%case 1
P2 = P2max;
P1 = P1max;
if(s < secrecy1(P1,P2,hab,hae,hbe,sig1,sig2,sige,sigl1,sigl2))
    s = secrecy1(P1,P2,hab,hae,hbe,sig1,sig2,sige,sigl1,sigl2);
    optP1 = P1;
    optP2 = P2;
end

%case 2
P2 = P2max;
P1 = 0;
if(s < secrecy1(P1,P2,hab,hae,hbe,sig1,sig2,sige,sigl1,sigl2))
    s = secrecy1(P1,P2,hab,hae,hbe,sig1,sig2,sige,sigl1,sigl2);
    optP1 = P1;
    optP2 = P2;
end

%case 2
P1 = P1max;
P2 = 0;
if(s < secrecy1(P1,P2,hab,hae,hbe,sig1,sig2,sige,sigl1,sigl2))
    s = secrecy1(P1,P2,hab,hae,hbe,sig1,sig2,sige,sigl1,sigl2);
    optP1 = P1;
    optP2 = P2;
end
 
end
