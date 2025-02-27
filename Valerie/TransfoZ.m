function [res] = TransfoZ(fonction, frequence)
% argument la fonction de transfert et la frequence dechantillonnage
% la fonction retourne le denominateur et le numerateur dans une table

    % Definir la frequence d'echantillonnage
    Ts = (1/frequence);

    % Conversion en transformee en z
    TransZ = c2d(fonction, Ts, 'zoh');

    % Mettre dans la table
    res = [];
    res(1,:) = cell2mat(TransZ.Numerator);
    res(2,:) = cell2mat(TransZ.Denominator);


end
