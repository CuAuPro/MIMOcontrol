function J = cenilka_maciejowski(param)

    global model out G_PI

    gamma = param(1);
    delta = param(2);
    
    G1 = inv(model.C*inv((-model.A))*model.B+model.D);

    w = 10^-3; % Izberemo!!
    resp = freqresp(model,w);
    G2 = align(resp);
    

    Kp = gamma*G1;
    Ki = delta*G2;
    s  = tf('s');

    G_PI = Kp + Ki/(s);

      
    out = sim("uglasevanje.slx");
    
    
    J = sum(abs(out.r1.Data - out.y1_lin.Data)) + sum(abs(out.r2.Data - out.y2_lin.Data));
  