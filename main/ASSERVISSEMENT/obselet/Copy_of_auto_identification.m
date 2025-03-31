clear
%% Identification automatique de procédé
%load data_log.mat
%test = table2array(datalog);
addpath("support")
addpath("Data")
%% Sortir les données 
cut = 1;
fin = 259;
res_exp = read_csv("Essai1_0.6V.csv");
y3 = res_exp.T3(cut:fin, 1);
y2 = res_exp.T2(cut:fin, 1);
y1 = res_exp.T1(cut:fin, 1);
t = res_exp.Temps(cut:fin, 1);
u = ones(size(y1))*0.6;
u(1,1) = 0;

Gain_tec = 2.2;

%remettre à 0
y1 = y1 - y1(1, 1); % Retirer les point d'opération
y2 = y2 - y2(1, 1); % Retirer les point d'opération
y3 = y3 - y3(1, 1); % Retirer les point d'opération

%% Identification des modèles

% P -> T1
[gp_1, modele_1] = identify(y1, u, t, 1, 0, true, NaN);
    assignin('base', 'gp1', gp_1);
gp_1 = gp_1/Gain_tec; % On divise par gain tec pour isoler gp_1

% T1 -> T2
[gp_2, modele_2] = identify(y2, y1, t, 1, 0, true, 5);
    assignin('base', 'gp2', gp_2);

% T2 -> T3
[gp_3, modele_3] = identify(y3, y2, t, 1, 0, false, 15);
    assignin('base', 'gp3', gp_3);

%bode(gp_3);
[gp_4, modele_4] = identify(y3, y1, t, 2, 0, false, 20);
    assignin('base', 'gp4', gp_4);

%% Identifier le controleur
chose = identify(y3, u, t, 1, 0, false, 20); % La ft tension -> T3

%% Avec placement de pôles
[Controleur, Ti, Td, Tf, kc] = pidfPoleCancellation(chose, 1);

%assignin('base', 'procede', procede);
assignin('base', 'Controleur', Controleur);

assignin('base', 'Ti', Ti);
assignin('base', 'Td', Td);
assignin('base', 'Tf', Tf);
assignin('base', 'kc', kc);

%% Coeff de S-H
R_nom = 10000;
% Coeffe S-H r2t
A_r2t = 0.00335401643468053;
B_r2t = 0.000256523550896126;
C_r2t = 0.00000260597012072052;
D_r2t = 0.000000063292612648746;
assignin('base', 'A_r2t', A_r2t);
assignin('base', 'B_r2t', B_r2t);
assignin('base', 'C_r2t', C_r2t);
assignin('base', 'D_r2t', D_r2t);
assignin('base', 'R_nom', R_nom);
% Coeffe S-H t2r
A_t2r =  -14.65719769;
B_t2r = 4798.84200000;
C_t2r = -115334.00000000;
D_t2r = -3730535.00000000;
assignin('base', 'A_t2r', A_t2r);
assignin('base', 'B_t2r', B_t2r);
assignin('base', 'C_t2r', C_t2r);
assignin('base', 'D_t2r', D_t2r);
assignin('base', 'R_nom', R_nom);
%% Procédé en Z
f = 1/5;
assignin('base', 'f', f);

[~, Z_T1T3] = s2z(gp_4, f, "zoh");
assignin('base', 'Z_T1T3', Z_T1T3);

[~, Z_T1] = s2z(gp_1, f, "zoh");
assignin('base', 'Z_T1', Z_T1);

[~, Z_T2] = s2z(gp_2, f, "zoh");
assignin('base', 'Z_T2', Z_T2);

[~, Z_T3] = s2z(gp_3, f, "zoh");
assignin('base', 'Z_T3', Z_T3);
%bode(Z_T3);
%% Controleur en Z
% TF I: 1 / (Ti + 1)
bloc_Ti = tf(1, [Ti, 1]);
[~, Tiz] = s2z(bloc_Ti, f, "zoh");
assignin('base', 'Tiz', Tiz);
% TF DF (Td + 1) / (Tf + 1)
bloc_DF = tf([Td, 1],[Tf, 1]);
[~, DFz] = s2z(bloc_DF, f, "tustin");
assignin('base', 'DFz', DFz);
% K = k_c
x = 2;
disp(x)
indent_pertub();
%assignin('base', 'filtre_consigne', filtre_consigne);

