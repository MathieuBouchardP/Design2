function [info, TransZ] = s2z(fonction, frequence, methode)
% Transforme une fonction de transfert continue en discrète (Z-transformée)
% Entrées :
%   - fonction : fonction de transfert continue (tf)
%   - frequence : fréquence d'échantillonnage en Hz
% Sortie :
%   - info : matrice contenant les coefficients normalisés du numérateur et du dénominateur

    % Définir la période d'échantillonnage
    Ts = 1 / frequence;

    % Conversion en transformée en Z (Tustin)
    TransZ = c2d(fonction, Ts, methode);
    %disp(TransZ);
    % Extraction des coefficients (SISO uniquement)
    num = cell2mat(TransZ.Numerator); 
    den = cell2mat(TransZ.Denominator);
    retard = TransZ.InputDelay + TransZ.IODelay + TransZ.OutputDelay;
    if (retard ~= 0)
    num = [zeros(1, retard), num];
    den = [den, zeros(1, retard)];
    end
    info = [num;den];
    %TransZ = tf(num, den);
end