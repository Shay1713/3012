%% QUESTION 1

% table parameters & measurements
g = 9.8; L = 11.2; mp = 381;
Ip = 0.00616; mw = 36; icmw = 0.00000746;
rw = 2.1; R = 4.4; kb = 0.444;
kt = 0.470;
% state space equations
q1 = icmw + rw * (mp * (L + rw) + mw*rw); % 1.0800e+04
q2 = icmw * (Ip + (L^2)*mp) + (rw^2) * (Ip * (mp + mw)) + (L^2) * mp * mw; % 1.7205e+06
q3 = Ip + L * mp * (L + rw); % 5.6754e+04

A = [0 1 0 0; (g*L*mp*((mp+mw)*(rw^2)+icmw))/q2 -1*(q1*kb*kt)/(q2*R) 0 -1*(q1*kb*kt)/(q2*R*rw); ...
    0 0 0 1; (g*(L^2)*(mp^2)*(rw^2))/q2 -1*(q3*kb*kt*rw)/(q2*R) 0 -1*(q3*kb*kt)/(q2*R)];
B = [0; -1*(q1*kt)/(q2*R); 0; -1*(q3*kt*rw)/(q2*R)];
C = [1 0 0 0];
D = 0;
% computing the values of the state model that are formulas
A1 = (g*L*mp*((mp+mw)*(rw^2)+icmw))/q2; % 44.6969
A2 = -1*(q1*kb*kt)/(q2*R); % -2.9771e-04
A4 = -1*(q1*kb*kt)/(q2*R*rw); % -1.4177e-04
A11 = (g*(L^2)*(mp^2)*(rw^2))/q2; % 457.3874
A22 = -1*(q3*kb*kt*rw)/(q2*R); % -0.0033
A44 = -1*(q3*kb*kt)/(q2*R); % -0.0016
B2 = -1*(q1*kt)/(q2*R); % -6.7051e-04
B4 = -1*(q3*kt*rw)/(q2*R); % -0.0074



%% QUESTION 2

format long
[num,den] = ss2tf(A, B, C, D);
transfer_function = tf(num,den);
transfer_function2 = C *inv(s * eye(4) - A)* B + D;
% disp(transfer_function2); 
% (12368749889199202*(2305843009213693952*s + 3607331900216889))/
% (- 42535295865117307932921825928971026432*s^3 - 79206583938566159591905553771659264*s^2 
% + 1901195190866511078379084520764361250635*s + 216230285280176978917243307261689856) 
% - 44618186041112393819948449944256/(- 42535295865117307932921825928971026432*s^3 
% - 79206583938566159591905553771659264*s^2 + 1901195190866511078379084520764361250635*s 
% + 216230285280176978917243307261689856)

% getting TF num values
num1 = num(1); num2 = num(2); % ... through num5 = num(5);
% s^4 0 , s^3 0 , s^2 -6.705112750399810e-04 , 
% s^1 2.554888831734579e-19 , s^0 -4.532455538564416e-19
% getting TF den values
den1 = den(1); den2 = den(2); den3 = den(3); den4 = num(4);
% s^3 0.001862137839356 , s^2 -44.696884133481682 , 
% s^1 -0.005083549576471 , s^0 0



%% QUESTION 3

figure;
impulse(transfer_function); 
title('Open-Loop Impulse Response');
% if need actual impulse function 
% disp(ilaplace(transfer_function2)); 



%% QUESTION 4

figure;
step(transfer_function); % plotting unit step
legend('step');title('Step Responce');xlabel('Time'); ylabel('Amplitude');
% if need actual impulse function
% disp(transfer_function2/s);

figure;
t = linspace(0,1,1000);
ramp = t;
y_ramp = lsim(transfer_function,ramp,t);
plot(t,y_ramp,'-'); hold on; plot(t,ramp,'--'); hold off;
title('Ramp Response'); xlabel('Time (s)'); ylabel('Output y(t)'); 
% if need actual ramp function
% disp(transfer_function2/(s^2));



%% QUESTION 5

figure;
pzmap(transfer_function);

[z,p,k] = tf2zp(num,den);
disp(p); % poles are ...
% 0
% -6.6864
% 6.6847
% -0.0001
%disp(z); % zeros are ...
%   1.0e-07 *
% (0.0000 + 0.2600i)
% (0.0000 - 0.2600i)



%% QUESTION 6

% Refrence Table 1



%% QUESTION 7

K = 1; % also tried K = 5, 50, and 1/4
z = roots([0 -6.6864 6.6847 -0.0001 K]);
figure;
plot(real(z), imag(z), 'kx');
title('Poles when K=1'); 

figure;
sys = tf(num,den);
rlocus(sys);



%%% QUESTION 8

% Plugged values from first column of RH Table
% into variables based on column 1 and row 
% number and tried different gain values of K
% to see if could make all column 1 positive 

% Example
% K = 100; % no matter what K was could not stabilize system 
% col1row3 = ((0.001862137839356)*(-44.6969-K*(6.705112750399810e-04))-(1)
%          *(-0.0051+K*(2.554888831734579e-19)))/(0.001862137839);
% col2row4 = ( (( (0.001862137839356)*(-44.6969 - K *(6.705112750399810e-04)) 
% - (1)*(-0.0051 + K * (2.554888831734579e-19)) ) / (0.001862137839))*(-0.0051 
% + K * (2.554888831734579e-19)) - (0.001862137839)*(( (-0.0051 + K * (2.554888831734579e-19))
% *(-K * (4.532455538564416e-19)) - 0 ) / (-0.0051 + K * (2.554888831734579e-19))) ) / 
% (((0.001862137839356)*(-44.6969-K*(6.705112750399810e-04))-(1)*(-0.0051+K*(2.554888831734579e-19)))/(0.001862137839));
% col1row5 = ( (( (( (0.001862137839356)*(-44.6969 - K *(6.705112750399810e-04)) - (1)*(-0.0051 + K 
% * (2.554888831734579e-19)) ) / (0.001862137839))*(-0.0051 + K * (2.554888831734579e-19)) - (0.001862137839)*
% (((-0.0051 + K * (2.554888831734579e-19))*(-K * (4.532455538564416e-19)) - 0 ) / (-0.0051 + K * (2.554888831734579e-19))) ) 
% / (( (0.001862137839356)*(-44.6969 - K *(6.705112750399810e-04)) - (1)*(-0.0051 + K * (2.554888831734579e-19)))/(0.001862137839)))
% *(( (-0.0051 + K * (2.554888831734579e-19))*(-K * (4.532455538564416e-19)) - 0 ) / (-0.0051 + K * (2.554888831734579e-19))) - 0 ) 
* / (( (( (0.001862137839356)*(-44.6969 - K *(6.705112750399810e-04)) - (1)*(-0.0051 + K * (2.554888831734579e-19)) ) / 
* (0.001862137839))*(-0.0051 + K * (2.554888831734579e-19)) - (0.001862137839)*(( (-0.0051 + K * (2.554888831734579e-19))*
% (-K*(4.532455538564416e-19)) - 0 ) / (-0.0051 + K * (2.554888831734579e-19))) ) / (( (0.001862137839356)*(-44.6969 - K *(6.705112750399810e-04)) 
% - (1)*(-0.0051 + K * (2.554888831734579e-19)) ) / (0.001862137839)) );

% Reference Table 2


