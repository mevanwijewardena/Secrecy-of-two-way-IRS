function [s,opt_P1,opt_P2] = LAM(rounds,L,P1max,P2max,sig1,sig2,sige,sigl1,sigl2,H,G,E1,E2,thresh1)

optimal = 0;
opt_W   = zeros(L+1,L+1);
opt_P1  = 0;
opt_P2  = 0;

wx = reshape(linspace(0,2*pi,L+2),[L+2,1]);
w0 = exp(1i*wx(2:L+2,:));
   
W0 = w0*w0';
    
inner_thresh = 0.001;
outer_thresh = 0.001;
steady_num = 1;
series = 0;
ini_value = -1000000;
win_stre = 0;
i = 1;
while i<rounds && series<steady_num
        if(win_stre>10)
            break;
        end
        inner_series = 0;
        value = -10000000;
        %optimizing power
        [s,P1,P2] = Power(W0,H,G,E1,E2,sig1,sig2,sige,sigl1,sigl2,P1max,P2max);
        while inner_series<steady_num
            div = 1+((P1/(sige^2))*(trace(E1*W0)))+((P2/(sige^2))*(trace(E2*W0)));
            cvx_begin quiet
                cvx_solver mosek
                variable W(L+1,L+1) complex semidefinite 
                %minimize(norm(W))
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
        end
   
        
        
        sec = ((log(1 + (P2/(sig1^2+sigl1^2))*real(trace(H*W0)))+log(1 + (P1/(sig2^2+sigl2^2))*real(trace(G*W0)))))-log(1+(P1/sige^2)*real(trace(E1*W0))+(P2/sige^2)*real(trace(E2*W0)));

        if(sec>optimal)
            optimal = sec;
            opt_W   = W0;
            opt_P1  = P1;
            opt_P2  = P2;
            win_stre =0;
        else
            win_stre = win_stre+1;
        end
        if(abs(sec-ini_value)<outer_thresh)
            series = series+1;
        end
        ini_value = sec;
       %x(i)= i;
        %y(i)= sec;
        %y1(i)=s1;
        %y2(i)=s2;
        %plot(x(1:i),y(1:i));
        
        %hold on
        %plot(x(1:i),y1(1:i));
        %hold on
        %plot(x(1:i),y2(1:i));
        %hold off
        %xlim([0 100])
        %ylim([0 40])
        i = i+1;
        %drawnow
end  

[U,S,V] = svd(opt_W);
w0 = U(:,1)/U(L+1,1);
w0 = w0./abs(w0);
Wx = w0*w0';

s = ((log(1 + (opt_P2/(sig1^2+sigl1^2))*real(trace(H*Wx)))+log(1 + (opt_P1/(sig2^2+sigl2^2))*real(trace(G*Wx)))))-log(1+(opt_P1/sige^2)*real(trace(E1*Wx))+(opt_P2/sige^2)*real(trace(E2*Wx)));
s = s/(log(2));
s = max(0,s);
end
