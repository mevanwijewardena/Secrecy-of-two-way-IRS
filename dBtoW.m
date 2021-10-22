function [watt] = dBtoW(x)
    watt = (10^(0.1*x))/1000;
end