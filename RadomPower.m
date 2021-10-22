function [optimal,optP1,optP2] = RadomPower(inits,rounds,L,P1max,P2max,sig1,sig2,sige,sigl1,sigl2,H,G,E1,E2,thresh1)
optP1 = 0;
optP2 = 0;
optimal = 0;
for i=1:inits
    opt_W   = zeros(L+1,L+1);
    opt_P1  = 0;
    opt_P2  = 0;
    P1 = P1max;
    P2 = P2max;
    wx = reshape(linspace(0,2*pi,L+2),[L+2,1]);
    w0 = exp(1i*wx(2:L+2,:));

    W0 = w0*w0';

    inner_thresh = 0.001;

    steady_num = 2;
    i = 1;
    inner_series = 0;
    value = -10000000;
    while i<rounds && inner_series<steady_num
         div = 1+((P1/(sige^2))*(trace(E1*W0)))+((P2/(sige^2))*(trace(E2*W0)));
         cvx_begin quiet
             cvx_solver mosek
             variable W(L+1,L+1) complex semidefinite 
             maximize(((log(0.0001 + (P2/(10^4*(sig1^2+sigl1^2)))*trace_mpower(W,1,H))+log(0.0001 + (P1/(10^4*(sig2^2+sigl2^2)))*trace_mpower(W,1,G))))-real((((P1/(sige^2))*(trace(E1*W)))+((P2/(sige^2))*(trace(E2*W))))/div))
             subject to
               diag(W)==1;
               %norm(trace((P1*E1+P2*E2).'*(W-W0)),1)<=(sige^2+real(P1*trace(E1*W0))+real(P2*trace(E2*W0)))*thresh1;
               norm((W-W0),1)<=thresh1;
          cvx_end
          W0 = W;
          if(abs(cvx_optval-value)<inner_thresh)
             inner_series = inner_series+1;
          end
          value = cvx_optval;
          i = i+1;
    end
    [U,S,V] = svd(W0);
    w0 = U(:,1)/U(L+1,1);
    w0 = w0./abs(w0);
    Wx = w0*w0';
    s = ((log(1 + (P2/(sig1^2+sigl1^2))*real(trace(H*Wx)))+log(1 + (P1/(sig2^2+sigl2^2))*real(trace(G*Wx)))))-log(1+(P1/sige^2)*real(trace(E1*Wx))+(P2/sige^2)*real(trace(E2*Wx)));
    s = s/(log(2));
    s = max(0,s);
    if(optimal<s)
        optimal = s;
        optP1 = P1;
        optP2 = P2;
    end
end  

end
