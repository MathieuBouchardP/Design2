function auto_ident_for_app_v2(csvfile, start_time, end_time, power_in, list_of_options)
%% Identification automatique de procédé
addpath("support");
addpath("Data");
% List of options:
% 1: (bool) ident. de T1 
% 2: (bool) ident. de T2
% 3: (bool) ident. de T3 
% 4: (bool) ident. de perturb_T1
% 5: (bool) ident. de perturb_T2
% 6: (bool) ident. de perturb_T3
% 7: (bool) calcul du régulateur

%% Sortir les données 
res_exp = read_csv(csvfile);
y3  = res_exp.T3(1:end, 1);
y2  = res_exp.T2(1:end, 1);
y1  = res_exp.T1(1:end, 1);
u   = power_in;
t   = res_exp.Temps(1:end, 1);

if ~isnan(start_time)
    id_start = find(t >= start_time, 1);
    y3  = y3(id_start:end, 1);
    y2  = y2(id_start:end, 1);
    y1  = y1(id_start:end, 1);
    u   = u(id_start:end, 1);
    t   = t(id_start:end, 1);
end
if ~isnan(end_time)
    id_end = find(t >= end_time, 1);
    y3  = y3(1:id_end, 1);
    y2  = y2(1:id_end, 1);
    y1  = y1(1:id_end, 1);
    u   = u(1:id_end, 1);
    t   = t(1:id_end, 1);
end


%% Retirer les points d'opération
y1      = y1 - y1(1, 1); % Retirer les point d'opération
y2      = y2 - y2(1, 1); % Retirer les point d'opération
y3      = y3 - y3(1, 1); % Retirer les point d'opération
t       = t - t(1,1);     % idem
u(1,1)  = 0;         % pour l'identification il faut un zéro au moins au départ
%% Identification des modèles
f =  evalin('base', 'f');
%%%%%%%%%%% P -> T1 %%%%%%%%%%%
if list_of_options(1)
    gp_1 = identify(y1, u, t, 1, 0, false, NaN);
        assignin('base', 'gp1', gp_1);
end
%%%%%%%%%%% T1 -> T2 %%%%%%%%%%%
if list_of_options(2)
    gp_2 = identify(y2, y1, t, 1, 0, false, NaN);
        assignin('base', 'gp2', gp_2);
% else
%     gp_2 = evalin('base', 'gp_2');
end
%%%%%%%%%%% T2 -> T3 %%%%%%%%%%%
if list_of_options(3)
    gp_3 = identify(y3, y2, t, 1, 0, false, NaN);
    assignin('base', 'gp3', gp_3);
    % T2 -> T3 (en Z)
    [~, Z_T3] = s2z(gp_3, f, "zoh");
    assignin('base', 'Z_T3', Z_T3);
% else
%     gp_3 = evalin('base', 'gp_3');
end

%% Identifier le controleur
if list_of_options(7)
    
    % tension -> T3
    gp_for_controleur = identify(y3, u, t, 2, 0, false, NaN);
    %% Avec placement de pôles
    [~, Ti, Td, Tf, kc] = pidfPoleCancellation(gp_for_controleur, 1);
    assignin('base', 'kc', kc);

    %% Controleur en Z
    % Bloc Ti: 1 / (Ti + 1)
    bloc_Ti = tf(1, [Ti, 1]);
    [~, Tiz] = s2z(bloc_Ti, f, "zoh");
    assignin('base', 'Tiz', Tiz);

    % Bloc DF: (Td + 1) / (Tf + 1)
    bloc_DF = tf([Td, 1],[Tf, 1]);
    [~, DFz] = s2z(bloc_DF, f, "tustin");
    assignin('base', 'DFz', DFz);
    % k_cz = k_c
end

end