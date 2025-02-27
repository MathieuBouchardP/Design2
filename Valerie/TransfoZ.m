function info = TransfoZ(fonction, frequence)
% argument la fonction de transfert et la fréquence d'échantillonnage
% la fonction redonne le numérateur et le dénominateur en transformée en z

% Définir la fréquence d'échantillonnage
Ts = (1/frequence);

% Conversion en transformée en z
TransZ = c2d(fonction, Ts, 'zoh');
    
num = cell2mat(TransZ.Numerator);
den = cell2mat(TransZ.Denominator);
info = [num;den];
end

