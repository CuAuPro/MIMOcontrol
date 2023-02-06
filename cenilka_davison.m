function J = cenilka_davison(param)

    global model out G_PI

    gamma = param(1);
    delta = param(2);
    
    G0_inv = inv(model.C*inv((-model.A))*model.B+model.D);
    Kp = gamma*G0_inv;
    Ki = delta*G0_inv;
    s  = tf('s');

    G_PI = Kp + Ki/(s);

      
    out = sim("uglasevanje.slx");
    
    
    J = sum(abs(out.r1.Data - out.y1_lin.Data)) + sum(abs(out.r2.Data - out.y2_lin.Data));
  