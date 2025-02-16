% Définir une fonction de transfert

f = tf(2,[1 3]);

Ts = 0.1;

H_z = c2d(f, Ts, 'zoh')

disp(H_z.Numerator)
disp(H_z.Denominator)

% voir reponse a un echelon de puissance
step(f,'-', H_z, '--r')

% verification que la fonction marche

tz = TransfoZ(f, 10)
tz_num = tz(1,:);
tz_den = tz(2,:);

% display pour voir si meme valeur
% disp()