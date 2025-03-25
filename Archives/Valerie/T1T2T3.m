%% Fonctions de transfert de T1, T2 et T3
G1 = tf(11,[105 1]);
G2 = tf(0.78,[37 1],'InputDelay',5);
G3 = tf(0.87,[6.7 1],'InputDelay',13);

% Transformée des transformées de laplace en Z avec une fréquence de 10Hz
f = 10;
% redonne l'info en matrice (double)
% première rangée = numérateur, deuxième rangée = dénominateur
Z1 = TransfoZ(G1, f);
Z2 = TransfoZ(G2, f);
Z3 = TransfoZ(G3, f);
