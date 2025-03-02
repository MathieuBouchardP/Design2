function info = TransfoZ(fonction, frequence)
% argument la fonction de transfert et la fréquence d'échantillonnage
% la fonction redonne le numérateur et le dénominateur en transformée en z

% Définir la fréquence d'échantillonnage
Ts = (1/frequence);

% Conversion en transformée en z
TransZ = c2d(fonction, Ts, 'tustin');
    
num = cell2mat(TransZ.Numerator);
den = cell2mat(TransZ.Denominator);

% matrice avec coefficient de puissance decroissante
flipp = flip([num;den],2);
% diviser par le coefficient de z à la 0
a = flipp(1,1);
info = flipp/a;
% conversion de matrice en tableau
% utiliser la methode : array2table(info); mais bof cest mieux une matrice
end

