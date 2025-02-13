clear
%% Identification automatique de procédé
load data_log.mat
test = table2array(datalog);
%% Dossier d'enregistrement
save_in = "Identified_models";
base_file_name = "Identifié";



%% Aller chercher les valeurs
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

%% Identifier le modèle
gp = identify(y2, u, t, 1, 0, true);

gp_num = gp(1, :);
gp_den = gp(2, :);
gp_retard = gp(3, 1);
procede = tf(gp_num, gp_den, 'InputDelay', gp_retard);

%% Sauvgarder
save_model(gp)


[gc, Ti]= actually_calculates_PID(gp);
gc_num = gc(1, :);
gc_den = gc(2, :);
Controleur = tf(gc_num, gc_den); % Création du TF dans MATLAB
filtre_consigne = tf(1, [Ti 1]);

save_controler(gc)

assignin('base', 'procede', procede);
assignin('base', 'Controleur', Controleur);
assignin('base', 'filtre_consigne', filtre_consigne);

