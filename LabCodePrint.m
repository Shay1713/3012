clear all
close all

%% LAB 1

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

% Reference Table 2

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
% * / (( (( (0.001862137839356)*(-44.6969 - K *(6.705112750399810e-04)) - (1)*(-0.0051 + K * (2.554888831734579e-19)) ) / 
% * (0.001862137839))*(-0.0051 + K * (2.554888831734579e-19)) - (0.001862137839)*(( (-0.0051 + K * (2.554888831734579e-19))*
% (-K*(4.532455538564416e-19)) - 0 ) / (-0.0051 + K * (2.554888831734579e-19))) ) / (( (0.001862137839356)*(-44.6969 - K *(6.705112750399810e-04)) 
% - (1)*(-0.0051 + K * (2.554888831734579e-19)) ) / (0.001862137839)) );



%% LAB 2

%% QUESTION 1

% Found TF in lab 1 



%% QUESTION 2

% ROUTH TABLE CRITERION
% K > 1/G(s)
% K*z > 0
% Reference Table 3

% chose (s-z) as (s-6.68) , to cancel out 
% the pole in the right half of G(s)'s P-Z plot map
% kept K the same at -6.705
% ... in doing so, both conditons in R-H table are satisfied

figure;
rlocus(transfer_function);
title('Transfer Function Lab1 Root Locus');



%% QUESTION 3

[b,a] = eqtflength(num,den);
[z,p,k] = tf2zp(b,a);
%{
z = zeros = 
   1.0e-07 *
  0.000000001905180 + 0.259994244262245i
  0.000000001905180 - 0.259994244262245i
p = poles = 
                   0
  -6.686447109116259
   6.684698705145661
  -0.000113733868756
k = gain = 
    -6.705112750399810e-04   
%}

%tfH = 1*(s-6.686447109116259) * (10^(-7) *(0.000000001905180 + 0.259994244262245*i) 
% *(0.000000001905180 - 0.259994244262245*i)) / (s *s* (s + 6.686447109116259) * 
% (s - 6.684698705145661)*(s +0.000113733868756));

numTF = -6.705*10^(-7) *((s -0.000000001905180 + 0.259994244262245*i) *(s - 0.000000001905180-0.259994244262245*i));
numH = sym2poly(numTF);
denTF = (s *s* (s + 6.68) *(s +0.0001137));
denH = sym2poly(denTF);
transfer_functionH = tf(numH,denH);
%   -6.705e-07 s^2 + 2.555e-15 s - 4.532e-08
%  ----------------------------------------
%       s^4 + 6.68 s^3 + 0.0007595 s^2

figure;
pzmap(transfer_functionH);% plot poles\zeros
title('Transfer FunctionH Pole-Zero MapH');

[z,p,k] = tf2zp(numH,denH);
%disp('poles');disp(p); % poles are ...
%                   0
%                   0
%  -6.680000000000000
%  -0.000113700000000
%disp('zeros'); disp(z); % zeros are ...
%  0.000000001905180 + 0.259994244262245i
%  0.000000001905180 - 0.259994244262245i



%% QUESTION 4 

% Refernece Table 4

figure;
step(transfer_functionH); 
title('Step Responce');xlabel('Time'); ylabel('Amplitude');



%% QUESTION 5

figure;
rlocus(transfer_functionH);
title('Transfer FunctionH Root Locus');



%% QUESTION 6

% Reference Table 5

[b,a] = eqtflength(numH,denH);
[z,p,k] = tf2zp(b,a);
% z = 0.000000001905180 + 0.259994244262245i
%     0.000000001905180 - 0.259994244262245i
% p = -6.680000000000000
%     -0.000113700000000
% k = -6.705000000000000e-07

numPD = (-6.705*10^(-7) *s^2 + 2.555*10^(-15) *s - 4.532* 10^(-8) )*(s + 5.764);
numPDH = sym2poly(numPD);
transfer_functionPD = tf(numPDH,denH);
% -6.705e-07*s^3 + 3.865e-06*s^2 - 4.532e        -08 s + 2.612e-07
% -----------------------------------------
% s^4 + 6.68*s^3 + 0.0007595*s^2

numPDK = numPD * (-6.705000000000000e-07); % complex so using same K as before -6.705000000000000e-07
numPDHK = sym2poly(numPDK);
transfer_functionPD = tf(numPDHK,denH); 
%4.496e-13*s^3 + 2.591e-12*s^2 + 3.039e-14*s + 1.752e-13
% --------------------------------------
% s^4 + 6.68*s^3 + 0.0007595*s^2



%% QUESTION 7

closedLoopPDHK = transfer_functionPD / (1+ transfer_functionPD);
figure;
pzmap(closedLoopPDHK);
figure;
step(closedLoopPDHK);
figure;
rlocus(closedLoopPDHK);



%% QUESTION 8

[z,p,k] = tf2zp(numPDHK, denH);
% z =
%  -5.764000000000000 + 0.000000000000000i
%   0.000000001905295 + 0.259983364350492i
%   0.000000001905295 - 0.259983364350492i
% p =
%                    0
%                    0
%   -6.680000000000000
%   -0.000113700000000
% k =
%      4.495702500000000e-13

H = zpk([5.764000000000000 -0.000000001905295+0.259983364350492i ...
    -0.000000001905295-0.259983364350492i],[ 6.680000000000000 0.000113700000000] , 4.495702500000000e-13);

C = pidtune(H,'PID');
% C =
%             1          
%  Kp + Ki * --- + Kd * s
%             s          
%  with Kp = -3.45e+08, Ki = -699, Kd = -4.25e+13








