function J = cenilka_penttinen(param)

    global model out G_PI

    gamma = param(1);
    delta = param(2);
    
    if det(model.C*model.B) == 0
    G1 = inv(model.C*model.B+model.C*model.A*model.B);
    else
        G1 = inv(model.C*model.B);
    end
    
    G2 = inv(model.C*inv((-model.A))*model.B+model.D);

    Kp = gamma*G1;
    Ki = delta*G2;
    s  = tf('s');

    G_PI = Kp + Ki/(s);

      
    out = sim("uglasevanje.slx");
    
    
    J = sum(abs(out.r1.Data - out.y1_lin.Data)) + sum(abs(out.r2.Data - out.y2_lin.Data));
  