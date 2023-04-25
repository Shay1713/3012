%% QUESTION 1

% uncompensated
num = 6.128;
den = conv([1 0], conv([1 1], [1 3]));
tf0 = tf(num, den);
rlocus(tf0);
title('uncompensated TF root locus'); 

% compensated
figure;
num = conv(2.116, conv([1 6], [1 0.53215]));
den = conv([1 0.01], conv([1 1.28], conv([1 0], conv([1 1], [1 3]))));
tf1 = tf(num , den); 
rlocus(tf1);
title('compensated TF root locus'); 

%% QUESTION 3
num = 50;
den = conv([1 0], conv([1 3], [1 6]));
sys1 = tf(num, den);
figure;
nyquist(sys1);
title('nyqusit sys1');
[z1,p1,k1] = tf2zp(num,den);
% z1 = n/a % p1 = (0,0) (-6,0) (-3,0) % k1 = 50

num = conv(50, [1, 4]);
den = conv([1 0 0], [1 1]);
sys2 = tf(num, den);
figure;
nyquist(sys2);
title('nyqusit sys2');
[z2,p2,k2] = tf2zp(num,den);
% z2 = (0,-4) % p2 = (0,0) (0,0) (-1,0) % k2 = 50;

num = conv(20, [1 3]);
den = conv([1 0], conv([1 1], [1 4]));
sys3 = tf(num, den);
figure;
nyquist(sys3);
title('nyqusit sys3');
[z3,p3,k3] = tf2zp(num,den);
% z3 = (-3,0) % p3 = (0,0) (0,-4) (0,-1) % k3 = 20

num = conv(100, [1 5]); 
den = conv([1 0], conv([1 3], [1 0 4]));
sys4 = tf(num, den);
figure;
nyquist(sys4);
title('nyqusit sys4');
[z4,p4,k4] = tf2zp(num,den); 
% z4 = (-5,0) % p4 = (0,0) (-3,0) (0, j2) (0, -j2) % k4 = 100
