clear
filename = "data_log.mat";
load(filename);
vars = who('-file', filename);
disp(vars)
%     test = table2array(datalog);
%     addpath("support")
%     %% Dossier d'enregistrement
% save_in = "Identified_models";
% base_file_name = "Identifié";

% %% Sortir les données 
% cut = 10;
% fin = 1440+cut;
% res_exp = read_csv("data_new.csv");
% y3 = res_exp.Temp_0__C_(cut:fin, 1);
% y2 = res_exp.Temp_1__C_(cut:fin, 1);
% y1 = res_exp.Temp_2__C_(cut:fin, 1);
% 
% courant = res_exp.Courant_V_(cut:fin, 1);
% k_c = 1.34375;
% courant = courant/k_c;
% tension = res_exp.DeltaV_V_(cut:fin, 1);
% %u = times(courant, tension);
% Gain_tec = 2.1;
% u = tension;
% assignin('base', 'TEC', Gain_tec);
% 
% %normaliser les valeurs
% offset_t2 = y2(1, 1) - y1(1, 1);
% y2 = y2 - offset_t2;
% offset_t3 = y3(1, 1) - y1(1, 1);
% y3 = y3 - offset_t3;
% 
% %remettre à 0
% y1 = y1 - y1(1, 1); % Retirer les point d'opération
% y2 = y2 - y2(1, 1); % Retirer les point d'opération
% y3 = y3 - y3(1, 1); % Retirer les point d'opération
% 
% t =  res_exp.Temps_s_(cut:end, 1) - res_exp.Temps_s_(cut, 1);
% %u = res_exp.Echelon_V_(cut:end, 1);
% %u = u - u(end-10, 1);
% %% Aller chercher les valeurs
%     cut = 1;                                % L'échantillonnage a commencé avant l'échelon
%     t = test(cut:end, 1) - test(cut, 1); % Le temps
%     y1 = test(cut:end, 2) - test(cut, 2);    % La température T1
%     y1 = y1 - y1(1, 1); % Retirer les point d'opération
%     y2 = test(cut:end, 3) - test(cut, 3);    % la température T2
%     y2 = y2 - y2(1, 1); % Retirer les point d'opération
%     y3 = test(cut:end, 4) - test(cut, 4);    % La température T3
%     y3 = y3 - y3(1, 1); % Retirer les point d'opération
% 
% 
%     %% initialisation de la consigne
%     n_zero = 21;
%     echelon = 2.2;
%     N = size(y1, 1)-1;
%     u = [zeros(21, 1) ; ones(N-n_zero+1, 1)] * echelon; % création d'un vecteur de la consigne
