clear
%% Identification automatique de procédé
load data_log.mat
test = table2array(datalog);
addpath("support")
%% Dossier d'enregistrement
save_in = "Identified_models";
base_file_name = "Identifié";

%% Sortir les données 
cut = 10;
fin = 1440+cut;
res_exp = read_csv("data_new.csv");
y3 = res_exp.Temp_0__C_(cut:fin, 1);
y2 = res_exp.Temp_1__C_(cut:fin, 1);
y1 = res_exp.Temp_2__C_(cut:fin, 1);

courant = res_exp.Courant_V_(cut:fin, 1);
k_c = 1.34375;
courant = courant/k_c;
tension = res_exp.DeltaV_V_(cut:fin, 1);
%u = times(courant, tension);
Gain_tec = 2.1;
u = tension;
assignin('base', 'TEC', Gain_tec);

%normaliser les valeurs
offset_t2 = y2(1, 1) - y1(1, 1);
y2 = y2 - offset_t2;
offset_t3 = y3(1, 1) - y1(1, 1);
y3 = y3 - offset_t3;

%remettre à 0
y1 = y1 - y1(1, 1); % Retirer les point d'opération
y2 = y2 - y2(1, 1); % Retirer les point d'opération
y3 = y3 - y3(1, 1); % Retirer les point d'opération

t =  res_exp.Temps_s_(cut:end, 1) - res_exp.Temps_s_(cut, 1);
%u = res_exp.Echelon_V_(cut:end, 1);
%u = u - u(end-10, 1);
%% Aller chercher les valeurs
bin = false;
if bin == true
    cut = 1;                                % L'échantillonnage a commencé avant l'échelon
    t = test(cut:end, 1) - test(cut, 1); % Le temps
    y1 = test(cut:end, 2) - test(cut, 2);    % La température T1
    y1 = y1 - y1(1, 1); % Retirer les point d'opération
    y2 = test(cut:end, 3) - test(cut, 3);    % la température T2
    y2 = y2 - y2(1, 1); % Retirer les point d'opération
    y3 = test(cut:end, 4) - test(cut, 4);    % La température T3
    y3 = y3 - y3(1, 1); % Retirer les point d'opération

    
    %% initialisation de la consigne
    n_zero = 21;
    echelon = 2.2;
    N = size(y1, 1)-1;
    u = [zeros(21, 1) ; ones(N-n_zero+1, 1)] * echelon; % création d'un vecteur de la consigne
end
%% Identification des modèles

display = true;
chose = true;
%% début
% y1, y2, y3, u, t
if chose == true
% P -> T1
[gp_1, modele_1] = identify(y1, u, t, 1, 0, false, NaN);
    assignin('base', 'gp1', gp_1);
gp_1 = gp_1 / Gain_tec;
% T1 -> T2
[gp_2, modele_2] = identify(y2, y1, t, 1, 0, false, 5);
    assignin('base', 'gp2', gp_2);
 
% T2 -> T3
[gp_3, modele_3] = identify(y3, y2, t, 1, 0, false, 15);
    assignin('base', 'gp3', gp_3);

[gp_4, modele_4] = identify(y3, y1, t, 2, 0, false, 20);
    assignin('base', 'gp3', gp_3);
end

%% Identifier le controleur
%[procede, ~] = identify(y3, u, t, 2, 0, false, 20);

procede = gp_1 * gp_2 * gp_3;
procede.IODelay = gp_1.IODelay + gp_2.IODelay + gp_3.IODelay;
order = 2;
%chose = balred(procede, order);
chose = identify(y3, u, t, 2, 0, false, 20);
%chose.IODelay = 20;

%% Avec placement de pôles
[Controleur, Ti, Td, Tf, kc] = pidfPoleCancellation(chose, 1.22);

assignin('base', 'procede', procede);
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

[~, Z_T1T3] = s2z(gp_4, f, "tustin");
assignin('base', 'Z_T1T3', Z_T1T3);

[~, Z_T1] = s2z(gp_1, f, "tustin");
assignin('base', 'Z_T1', Z_T1);

[~, Z_T2] = s2z(gp_2, f, "tustin");
assignin('base', 'Z_T2', Z_T2);

[~, Z_T3] = s2z(gp_3, f, "tustin");
assignin('base', 'Z_T3', Z_T3);
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


%assignin('base', 'filtre_consigne', filtre_consigne);

