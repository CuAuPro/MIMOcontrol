%%
clear all;close all;

%% System Parameters

% Physical constants
r = 0.07; %m;
S = 0.0152; %m

% Steady states
u1_ss = -4; %V (napetost vzbujanja prve crpalke v delovni tocki)
u2_ss = -3.5; %V (napetost vzbujanja druge crpalke v delovni tocki)
x1_ss = 0.7712; %V (napetost prvega senzorja v delovni tocki)
x2_ss = 3.2841; %V (napetost drugega senzorja v delovni tocki)
x3_ss = 4.8098; %V (napetost tretjega senzorja v delovni tocki)
% oziroma
h1_ss = 0.3030; % m (normalna visina vode v prvem shranjevalniku)
h2_ss = 0.2290; % m  (normalna visina vode v drugem shranjevalniku)
h3_ss = 0.1520; % m  (normalna visina vode v tretjem shranjevalniku)

integrator1_init = h1_ss; % zacetno stanje integratorja
integrator2_init = h2_ss; % zacetno stanje integratorja
integrator3_init = h3_ss; % zacetno stanje integratorja

% Constants
Ks1 = -27.6937; %V/m
Ks2 = -27.5827; %V/m
Ks3 = -27.5087; %V/m
Kos1 = 9.1500; %V
Kos2 = 9.6000; %V
Kos3 = 9.0000; %V
K32 = -0.1107*10e-3; % m/s
K31 = 0.1747*10e-3; % m/s
K30 = 0.0360*10e-3; % m/s

%% Curve fitting
pump1_u = [-10 -9 -8 -7 -6 -5 -4.5 -4.4 -4 -3.85 -3.6 -3.5 -3.15 -3 -2.5 -2.3 -2.1 -1.9 -1.7 -1.5 -1.3 -1.1 -0.9 -0.7 -0.5 -0.3 -0.1 0 1 2 3 4 5 6 7 8 9 10];
pump1_flow = [0.00000000000 0.00000000000 0.00000689390 0.00001238400 0.00001775100 0.00002311700 0.00002571800 0.00002630800 0.00002836600 0.00002892000 0.00003038100 0.00003094700 0.00003231900 0.00003196800 0.00003715200 0.00003888600 0.00004000100 0.00004024900 0.00004123900 0.00004161100 0.00004173500 0.00004179700 0.00004185800 0.00004303500 0.00004532600 0.00004644100 0.00004829800 0.00004935100 0.00005461400 0.00005981600 0.00006390200 0.00006922700 0.00007554300 0.00008024900 0.00008384100 0.00008507900 0.00009312900 0.00009783500];

pump2_u = [-10 -9 -8 -7 -6 -5 -4.5 -4.4 -4 -3.85 -3.6 -3.5 -3.15 -3 -2.5 -2.3 -2.1 -1.9 -1.7 -1.5 -1.3 -1.1 -0.9 -0.7 -0.5 -0.3 -0.1 0 1 2 3 4 5 6 7 8 9 10];
pump2_flow = [0.0000000000 0.0000000000 0.0000077952 0.0000136100 0.0000192980 0.0000248600 0.0000278520 0.0000277270 0.0000298140 0.0000305550 0.0000321400 0.0000325130 0.0000338560 0.0000350370 0.0000428520 0.0000440000 0.0000439770 0.0000448750 0.0000449380 0.0000450010 0.0000450010 0.0000450010 0.0000445510 0.0000450010 0.0000472770 0.0000496780 0.0000517010 0.0000523960 0.0000580210 0.0000630780 0.0000678180 0.0000738220 0.0000786260 0.0000863370 0.0000892440 0.0000915200 0.0001006200 0.0001041600];

valve3_h = [0.0000000 0.0079540 0.0146260 0.0218480 0.0294610 0.0374640 0.0457510 0.0543390 0.0632290 0.0724740 0.0820920 0.0919580 0.1021100 0.1125600 0.1233500 0.1344200 0.1458500 0.1575600 0.1695700 0.1819400 0.1945900 0.2075700 0.2208900 0.2345000 0.2484300 0.2627400 0.2772300 0.2920500 0.3072600 0.3227100 0.3385600 0.3547400 0.3711600 0.3878900 0.4050500 0.4224900 0.4400600 0.4580200 0.4762400 0.4948000 0.5137500 0.5330300 0.5523100 0.5716200 0.5907100 0.6000000];
valve3_flow = [0.0000000000 0.0000324280 0.0000360720 0.0000380760 0.0000400800 0.0000420840 0.0000429950 0.0000451810 0.0000460920 0.0000488250 0.0000499180 0.0000513750 0.0000528320 0.0000544720 0.0000562940 0.0000579610 0.0000599380 0.0000615080 0.0000630350 0.0000639460 0.0000659500 0.0000672250 0.0000695930 0.0000701400 0.0000728720 0.0000739650 0.0000748760 0.0000772450 0.0000788840 0.0000797950 0.0000828920 0.0000832570 0.0000852610 0.0000865360 0.0000896330 0.0000894510 0.0000909080 0.0000934590 0.0000936410 0.0000969200 0.0000976490 0.0000978000 0.0000988000 0.0000988000 0.0001004000 0.0001006000];
%% Calculate Valve parameters K1 and K2
%Valve parameters: Flows in steady state are exracted from measurement
%tables.

%flow3 in steady state is same as flow1
flow3_ss = 2.8366e-05;
K1 = flow3_ss/(sqrt(h1_ss-h2_ss));

%Equation on page 3 (fi3 - fi4 = Sdh_2/dt... Steady state --> fi3-fi4 = 0)
K2 = K1*(sqrt(h1_ss-h2_ss))/(sqrt(h2_ss-h3_ss));

%% Constrains
global u_min u_max
u_min = -10;
u_max = 10;
%% Linearized model
A = [-0.0125,  0.0126,  0.0000;
      0.0125, -0.0246,  0.0121;
      0.0000,  0.0120, -0.0212];
B = [-0.0091,  0.0000;
      0.0000,  0.0000;
      0.0000, -0.0092];
C = [1.0, 0.0, 0.0;
     0.0, 0.0, 1.0];
D = [0.0, 0.0;
     0.0, 0.0];
%x_initial = [x1_ss;
%             x2_ss;
%             x3_ss];
global model;
model = ss(A,B,C,D);
x_initial = [0;
             0;
             0];


%% Simulate parameters
tsim = 14000; %s

%PAZI: To je dejansko napetost na IZHODNEM senzorju. Torej referenca je
%napetost na izhodnem senzorju, ki meri visino.
r1_ss = x1_ss;
r2_ss = x3_ss;
d_r1 = 0.3; %V
d_r2 = 0.3; %V

r1_in = {};
r1_in.Vals = [r1_ss, r1_ss+d_r1, r1_ss-d_r1, r1_ss];
r1_in.Times = [0, 4000, 7000, 10000];
% u2_in = {};
% u2_in.Vals = [r2_ss, r2_ss+d_r2, r2_ss-d_r2, r2_ss];
% u2_in.Times = [0, 4000, 7000, 10000];
r2_in = {};
r2_in.Vals = [r2_ss];
r2_in.Times = [0];

regulator_type = 'davison' % izbire so 'davison', 'pettinen', 'maciejowski'
optimize_regulator = 1



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%        REGULATORS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Davison Regulator

if strcmp(regulator_type, 'davison')
    fprintf("Using regulator: %s.\n", regulator_type)
%Check conditions:
%1. m==l? YES - nr. columns in C is same as nr. rows in B
%2. Is system stable? Yes, confirmed in 'analysys.m'
%3. Is G(0) singular? No
if det(dcgain(A,B,C,D)) == 0
    fprintf("Matrika G(0) je singularna.\n")
else
    fprintf("Matrika G(0) ni singularna.\n")
end

if optimize_regulator == 1
    fprintf("Optimizing regulator.\n")
    %Optimisation
    global out G_PI
    gamma = 10;
    delta = 0.1;
    param0 = [gamma, delta];
    options = optimset('MaxIter', 200, ...
                      'Display', 'iter', ...
                      'TolFun', 1e-4, ...
                      'TolX', 1e-4);
    [paramOpt, J] = fmincon('cenilka_davison',...
       param0,[],[],[],[],[0, 0],[20, 1],...
       [],options)
    %Pred options gre lahko se 'pogoji_uglasevanje' (itak imamo saturation...)
    
    %Optimalni parametri za izbrano cenilko
    gamma = paramOpt(1) %15.0484
    delta = paramOpt(2) %0.0910
else
    fprintf("Not optimizing regulator.\n")
    gamma = 15.0484;
    delta = 0.0910;
end
G0_inv = inv(C*inv((-A))*B+D);

Kp = gamma*G0_inv;
Ki = delta*G0_inv;

s = tf('s');

G_PI = Kp + Ki/(s);
end

%% Penttinen Koivo Regulator
if strcmp(regulator_type, 'pettinen')
    fprintf("Using regulator: %s.\n", regulator_type)
%Check conditions:
%1. m==l? YES - nr. columns in C is same as nr. rows in B
%2. Is system stable? Yes, confirmed in 'analysys.m'
%3. Is G(0) singular? No
if det(dcgain(A,B,C,D)) == 0
    fprintf("Matrika G(0) je singularna.\n")
else
    fprintf("Matrika G(0) ni singularna.\n")
end
%4. Is D = 0? Yes
%5. Is C*B singular? No
if det(C*B) == 0
    G1 = inv(C*B+C*A*B);
else
    G1 = inv(C*B);
end

if optimize_regulator == 1
    fprintf("Optimizing regulator.\n")
    %%Optimisation
    global out G_PI
    gamma = 0.28;
    delta = 0.02;
    param0 = [gamma, delta];
    options = optimset('MaxIter', 200, ...
                       'Display', 'iter', ...
                       'TolFun', 1e-4, ...
                       'TolX', 1e-4);
    [paramOpt, J] = fmincon('cenilka_davison',...
        param0,[],[],[],[],[0, 0],[2, 1],...
        [],options)
    
    % Optimalni parametri za izbrano cenilko
    gamma = paramOpt(1) %0.2893
    delta = paramOpt(2) %0.0207
else
    fprintf("Not optimizing regulator.\n")
    gamma = 0.2893;
    delta = 0.0207;
end
G2 = inv(C*inv((-A))*B+D);

Kp = gamma*G1;
Ki = delta*G2;

s = tf('s');

G_PI = Kp + Ki/(s);
end
%% Maciejowski Regulator

if strcmp(regulator_type, 'maciejowski')
    fprintf("Using regulator: %s.\n", regulator_type)
%Check conditions:
%1. m==l? YES - nr. columns in C is same as nr. rows in B
%2. Is system stable? Yes, confirmed in 'analysys.m'
%3. Is G(0) singular? No
G1 = inv(-C*inv((A))*B+D);

if G1 == 0
    fprintf("Matrika G(0) je singularna.\n")
else
    fprintf("Matrika G(0) ni singularna.\n")
end
%3. Is G(jw) singular? No

w = 10^-3; % Izberemo!!
resp = freqresp(ss(A,B,C,D),w)
G2 = align(resp);
if det(G2) == 0
    fprintf("Matrika G(jw) je singularna.\n")
else
    fprintf("Matrika G(jw) ni singularna.\n")
end
if optimize_regulator == 1
    fprintf("Optimizing regulator.\n")
    %Optimisation
    global out G_PI
    gamma = 0.1;
    delta = 0.01;
    param0 = [gamma, delta];
    options = optimset('MaxIter', 200, ...
                       'Display', 'iter', ...
                       'TolFun', 1e-4, ...
                       'TolX', 1e-4);
    [paramOpt, J] = fmincon('cenilka_maciejowski',...
        param0,[],[],[],[],[0, 0],[1, 1],...
        [],options)
    
    % Optimalni parametri za izbrano cenilko
    gamma = paramOpt(1) %0.7549
    delta = paramOpt(2) %0.0540
else
    fprintf("Not optimizing regulator.\n")
    gamma = 0.1;
    delta = 0.0540;
end

Kp = gamma*G1;
Ki = delta*G2;

s = tf('s');

G_PI = Kp + Ki/(s);

end
%% Simulate


out = sim("uglasevanje.slx");
%% Results for validation

t = out.r1.Time;

%Referenci
r1 = out.r1.Data;
r2 = out.r2.Data;

%Regulirni signali
u1lin = out.u1_lin.Data;
u2lin = out.u2_lin.Data;
u1nelin = out.u1_nonlin.Data;
u2nelin = out.u2_nonlin.Data;

y1lin = out.y1_lin.Data;
y2lin = out.y2_lin.Data;
y1nelin = out.y1_nonlin.Data;
y2nelin = out.y2_nonlin.Data;

Jlin = sum(abs(r1 - y1lin)) + sum(abs(r2 - y2lin))
Jnelin = sum(abs(r1 - y1nelin)) + sum(abs(r2 - y2nelin))

%%
figure;

subplot(2,2,1);
hold on;
plot(t, r1);
plot(t, y1lin);
plot(t, y1nelin);
ylim([0,2]);
xlabel('Cas [s]', 'Interpreter', 'latex');
ylabel('Napetost [V]', 'Interpreter', 'latex');
legend('$r_1$', '$y_{1lin}$',' $y_{1nelin}$', 'Interpreter', 'latex');

subplot(2,2,2);
hold on;
plot(t, r2);
plot(t, y2lin);
plot(t, y2nelin);
ylim([4,6]);
xlabel('Cas [s]', 'Interpreter', 'latex');
ylabel('Napetost [V]', 'Interpreter', 'latex');
legend('$r_2$', '$y_{2lin}$',' $y_{2nelin}$', 'Interpreter', 'latex')

subplot(2,2,3);
hold on;
plot(t, u1lin);
plot(t, u1nelin);
yline([u_min u_max],'--');
ylim([-11, 11]);
xlabel('Cas [s]', 'Interpreter', 'latex');
ylabel('Napetost [V]', 'Interpreter', 'latex');
legend('$u_{1lin}$', '$u_{1nelin}$', 'Interpreter', 'latex');

subplot(2,2,4);
hold on;
plot(t, u2lin);
plot(t, u2nelin);
yline([u_min u_max],'--');
ylim([-11, 11]);
xlabel('Cas [s]', 'Interpreter', 'latex');
ylabel('Napetost [V]', 'Interpreter', 'latex');
legend('$u_{1lin}$', '$u_{2nelin}$', 'Interpreter', 'latex');