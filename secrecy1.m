function s = secrecy1(P1,P2,hab,hae,hbe,sig1,sig2,sige,sigl1,sigl2)
    s = log(1+((P2/(sig1^2+sigl1^2))*hab))+log(1+((P1/(sig2^2+sigl2^2))*hab))-log(1+((P1/sige^2)*hae)+((P2/sige^2)*hbe));
    s = max(0,s/(log(2)));
end