function [s,optP1,optP2] = Power(W,H,G,E1,E2,sig1,sig2,sige,sigl1,sigl2,P1max,P2max)

s= -1000000000;
optP1 = 0;
optP2 = 0;
%stationary point
P2 = (((sig2^2+sigl2^2)*trace_mpower(W,1,E1))-(sige^2*trace_mpower(W,1,G)))/(trace_mpower(W,1,G)*trace_mpower(W,1,E2));
P1 = (((sig1^2+sigl1^2)*trace_mpower(W,1,E2))-(sige^2*trace_mpower(W,1,H)))/(trace_mpower(W,1,H)*trace_mpower(W,1,E1));
if((P1>=0) && (P2>=0) && (secrecy(W,P1,P2,G,H,E1,E2,sig1,sig2,sige,sigl1,sigl2)>0))
    s = secrecy(W,P1,P2,G,H,E1,E2,sig1,sig2,sige,sigl1,sigl2);
    optP1 = P1;
    optP2 = P2;
end


%case 1
P2 = P2max;
P1 = P1max;
if(s < secrecy(W,P1,P2,G,H,E1,E2,sig1,sig2,sige,sigl1,sigl2));
    s = secrecy(W,P1,P2,G,H,E1,E2,sig1,sig2,sige,sigl1,sigl2);
    optP1 = P1;
    optP2 = P2;
end

%case 2
P2 = P2max;
P1 = 0;
if(s < secrecy(W,P1,P2,G,H,E1,E2,sig1,sig2,sige,sigl1,sigl2))
    s = secrecy(W,P1,P2,G,H,E1,E2,sig1,sig2,sige,sigl1,sigl2);
    optP1 = P1;
    optP2 = P2;
end

%case 2
P1 = P1max;
P2 = 0;
if(s < secrecy(W,P1,P2,G,H,E1,E2,sig1,sig2,sige,sigl1,sigl2))
   s = secrecy(W,P1,P2,G,H,E1,E2,sig1,sig2,sige,sigl1,sigl2);
   optP1 = P1;
   optP2 = P2;
end
 
optP1 = real(optP1);
optP2 = real(optP2);
end