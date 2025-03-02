function info = TransfoZ(fonction, frequence)
% Transforme une fonction de transfert continue en discrète (Z-transformée)
% Entrées :
%   - fonction : fonction de transfert continue (tf)
%   - frequence : fréquence d'échantillonnage en Hz
% Sortie :
%   - info : matrice contenant les coefficients normalisés du numérateur et du dénominateur

    % Définir la période d'échantillonnage
    Ts = 1 / frequence;

    % Conversion en transformée en Z (Tustin)
    TransZ = c2d(fonction, Ts, 'zoh');
    disp(TransZ);
    % Extraction des coefficients (SISO uniquement)
    num = cell2mat(TransZ.Numerator); 
    den = cell2mat(TransZ.Denominator);
    retard = TransZ.InputDelay;
    if (retard ~= 0)
    num = [zeros(1, retard), num];
    den = [den, zeros(1, retard)];
    end
    info = [num;den];
end


%% Fonctions de transfert de T1, T2 et T3
G1 = tf(11,[105 1]);
G2 = tf(0.78,[37 1],'InputDelay',5);
G3 = tf(0.87,[6.7 1],'InputDelay',15);
%% Fonctions de transfert pour le régulateur
DF = tf([10.203, 1],[5.1017, 1]);
Ti = tf(1, [138.62, 1]);
K = tf(0.05818, 1);
% Transformée des transformées de laplace en Z avec une fréquence de 10Hz
f = 1/5;

% première rangée = numérateur, deuxième rangée = dénominateur
Z2 = TransfoZ(G2*G3, f);
Z3 = TransfoZ(G3, f);
Zdf = TransfoZ(DF, f);
Zti = TransfoZ(Ti ,f);
Zk = TransfoZ(K, f);