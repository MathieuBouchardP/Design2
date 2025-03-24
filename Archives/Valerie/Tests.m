% Définir une fonction de transfert
f = tf([0 2],[1 3]);

% Définir le temps d'échantillonage (inverse de fréquence f = 1/T)
Ts = 0.1;

% Conversion Convert the transfer function to discrete-time using zero-order hold (ZOH)
H_z = c2d(f, Ts, 'tustin');


% Voir ce que ça fait 
% (en bleu la tf et en rouge la tz)
step(f,'-',H_z,'--r')

num = cell2mat(H_z.Numerator);
den = cell2mat(H_z.Denominator);
info = [num;den];

a = info(1,:);
b = info(2,:);
% disp(a)
% disp(b)
% test de la transformee en z
tz = TransfoZ(f,10);
disp(tz)

num_z = tz(1,:);
den_z = tz(2,:);
% disp(num_z)

% test pour inverser des arrays et un matrices
B = [0 0.1728];
A = flip(B);

D = [1 3 5; 2 4 6; 7 8 10];
C = flip(D,2);

disp(class(tz))