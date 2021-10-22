function s = secrecy(W,P1,P2,G,H,E1,E2,sig1,sig2,sige,sigl1,sigl2)
    s = log(1+((P2/(sig1^2+sigl1^2))*trace_mpower(W,1,H)))+log(1+((P1/(sig2^2+sigl2^2))*trace_mpower(W,1,G)))-log(1+((P1/sige^2)*trace_mpower(W,1,E1))+((P2/sige^2)*trace_mpower(W,1,E2)));
    s = s/(log(2));
end