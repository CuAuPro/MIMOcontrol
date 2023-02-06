function [c,ceq]=pogoji_uglasevanje(param)

    global out u_max u_min
 
c=    [u_min - min(out.u1_lin.Data);
    max(out.u1_lin.Data) - u_max;
    u_min - min(out.u2_lin.Data);
    max(out.u2_lin.Data) - u_max;]

ceq=[];