% Définir une fonction de transfert
f = tf(2,[1 3]);

% Définir le temps d'échantillonage (inverse de fréquence f = 1/T)
Ts = 0.1;

% Conversion Convert the transfer function to discrete-time using zero-order hold (ZOH)
H_z = c2d(f, Ts, 'zoh');


% Voir ce que ça fait (en bleu la tf et en rouge la tz)
step(f,'-',H_z,'--r')

num = cell2mat(H_z.Numerator);
den = cell2mat(H_z.Denominator);
info = [num;den];

a = info(1,:);
b = info(2,:);
disp(a)

% test de la transformee en z
tz = TransfoZ(f,10);
num_z = tz(1,:);
den_z = tz(2,:);
disp(num_z)

