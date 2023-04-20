%% Lecture 2

% Plotting a function and its Laplace transform
t = linspace(0,10,1000);
y = sin(pi*t) + sin(2*pi*t) + sin(4*pi*t) + sin(8*pi*t) + sin(10*pi*t);
w = fft(y);
% figure; stem(abs(w));

figure; 
subplot(1,2,1);
plot(t,y);
xlabel('Time');
title('Time domain');
subplot(1,2,2);
stem(fftshift(abs(w)));
xlabel('Frequency');
title('Frequency domain');

% Use Laplace function to compute Laplace transform of Exercise 2.1
syms  t
f = t*exp(-5*t); 
F = laplace(f);

% Use residue function to compute partial fraction decomposition of (-4s + 8)/(s^2 + 6s + 8) (Ex. 2.2)
% Multiply polynomials via conv function. 
a = [1 6 8];
b = [-4 8];
[r,p,k] = residue(b,a);



%% Lecture 3

syms s
A = [0 1; -2 0];
B = [0; 1];
C = [1 0];
D = 0;
transfer_function = C * inv(s * eye(2) - A)* B + D;
figure;
plot(transfer_function);



%% Lecture 4

syms s A y x u1 u2 B
A = [5 s 0 0 ; 0 s -1*(10/s + 2) 0 ; 1 -1 0 -1 ; 0 0 1 -1];
x = [u1 ; 0 ; 0 ; u2];
B = inv(A);
y = B*x;
transfer_function = s*y(2);
figure;
plot(transfer_function);



%% Lecture 5

sys = tf([1 3], [1 5]);

[num,den] = ss2tf([1 1; 0 1], [0 ; 1], [1 0], 0);

[A,B,C,D] = tf2ss(1, [1 1]);

% Pole-Zero Map
figure;
pzmap(sys);

% Impulse Response
figure;
impulse(sys);

% Step Response 
figure;
step(sys);

% Ramp Response
figure;
t = linspace(0,1,1000);
ramp = t;
y_ramp = lsim(g0,ramp,t);  % response of ramp input in time domain
plot(t,y_ramp);
xlabel('Time (s)');
ylabel('Output y(t)');
title('Ramp response');

% Series/cascade interconnection
sys1 = tf(1, [1 1]);
sys2 = tf(1, [1 4]);
sys3 = tf([1 3], [1 1 5]);
cascade_interconnect = series(sys1,series(sys2,sys3));

% Parallel interconnection
sys1 = tf(1, [1 1]);
sys2 = tf(1, [1 4]);
sys3 = tf([1 3], [1 1 5]);
parallel_interconnect = parallel(sys1, parallel(sys2,sys3));

% Feedback interconnection
sys1 = tf(1, [1 1]);
sys1_unity_FB = feedback(sys1,1);



%% Lecture 6

a = conv([1 0], conv([1 2], conv([1 4], [1 5])));
b = [1 3];
[r,p,k] = residue(b,a); 



%% Lecture 9

N = 100*conv([1 2], [1 6]);
D = conv([1 0], conv([1 3], [1 4]));
g = tf(N,D);
g0 = feedback(g,1);




%% Lecture 10

sys = tf(1, [1 10 0]);
rlocus(sys);

K = 0;
z = roots([1 10 K]);
figure;
plot(real(z), imag(z), 'kx');
title('Poles when K=0');



%% Lecture 11

plant = tf(1, conv([1 1], conv([1 2], [1 10])));
figure;
rlocus(plant);
sgrid;

figure;
syms K
closed_plant = feedback(K * plant, 1);
plot(closed_plant);

% search all K candidates
num_poles = 3;
zeta_target = 0.174;
K_candidates = linspace(1, 200, 1000);
[pole_final, zeta_final, K_final] = find_pole_zeta(plant, num_poles, zeta_target, K_candidates); %% check this function

% Find pole location and gain with damping ratio of 0.174
K = 165;
Kp = 91.6;
Kp_current = K/(2*10);
zc_pc = Kp/Kp_current;
pc = 0.01;
zc = zc_pc*pc;
controller = tf([1 zc], [1 pc]);
compensated_plant = series(controller,plant);
figure;
rlocus(compensated_plant);

plant = tf(1, conv([1 0], conv([1 4], [1 6])));
os_percent = 16;
zeta_target = zeta_from_OS(os_percent);
num_poles = 3;
K_candidates = linspace(1, 100, 2000);
[pole_final, zeta, K] = find_pole_zeta(plant, num_poles, zeta_target, K_candidates); 
figure;
rlocus(plant);
% Then find pole location with 16% overshoot
dominant_pole = -1.205+2.0645*j;
Ts = 4/1.205;
Ts_prime = Ts/3;
new_dominant_pole = -4/(Ts_prime) + j*(2.0645/1.205)*(4/Ts_prime);
angle_offset = angle(new_dominant_pole) + angle(new_dominant_pole+4) + angle(new_dominant_pole+6)-pi;
omega_d = imag(new_dominant_pole); sigma_d = real(new_dominant_pole);
z_c = omega_d/tan(angle_offset)-sigma_d;
compensated_plant = tf([1 z_c], conv([1 0], conv([1 4], [1 6])));
figure;
rlocus(compensated_plant);

% Lead compensator design example
figure;
plant = tf(1, conv([1 0], conv([1 4], [1 6])));
rlocus(plant);
% Find pole location with 30% overshoot
os_percent = 30;
zeta_target = zeta_from_OS(os_percent);
num_poles = 3;
K_candidates = linspace(1, 200, 1000);
[pole, zeta, K] = find_pole_zeta(plant, num_poles, zeta_target, K_candidates); 
dominant_pole = -1.01+j*2.63;
Ts = 4/1.01;
Ts_prime = Ts/2;
new_dominant_pole = -4/(Ts_prime) + j*(2.63/1.01)*(4/Ts_prime);
% Put new compensator zero at -5
z_c = 5;
angle_offset = -1*angle(new_dominant_pole+5)+angle(new_dominant_pole)+angle(new_dominant_pole+4)+angle(new_dominant_pole+6)-pi;
omega_d = imag(new_dominant_pole); sigma_d = real(new_dominant_pole);
p_c = -1*omega_d/tan(angle_offset)-sigma_d;
compensated_plant = tf([1 z_c], conv([1 p_c], conv([1 0], conv([1 4], [1 6]))));
figure;
rlocus(compensated_plant);



%% Lecture 12

plant = tf(1, conv([1 1], conv([1 2], [1 10])));
PI_controller = tf([1 0.1], [1 0.05]);
PI_plant = series(PI_controller, plant);
rlocus(PI_plant);

plant = tf(1, conv([1 1], conv([1 2], [1 10])));
PD_controller = tf([1 3], 1);
PD_plant = series(PD_controller, plant);
rlocus(PD_plant);

z_c = 57.3;
K = 5.08;
plant = tf([1 8], conv([1 3], conv([1 6], [1 10])));
compensator = tf(conv([1 0.5], [1 z_c]), [1 0]);
compensated_plant = series(compensator,plant);
rlocus(compensated_plant);

lead_compensated_tf = tf(1,conv([1 0], conv([1 10], [1 29])));
rlocus(lead_compensated_tf);



%% Lecture 13

K = 150;
sys = tf(K, conv([1 2], conv([1 3], [1 5])));
nyquist(sys);



%% Lecture 14

a=1;
system = tf(1, [1 a]);
figure;
bode(system);

zeta = 2;
omega_n = 1;
system = tf([1 zeta*omega_n omega_n^2], 1);
figure;
bode(system);

system = tf([1 3], conv([1 0], conv([1 1], [1 2])));
figure;
bode(system);

system = tf(6, conv([1 2], [1 2 2]));
figure;
bode(system);
[gm,pm,wcg,wcp] = margin(system);










