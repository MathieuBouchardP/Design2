function auto_ident_for_app(csvfile, u_echelon, deb, fin)
%% Identification automatique de procédé
addpath("support");
addpath("Data");

%% Sortir les données 
res_exp = read_csv(csvfile);
y3  = res_exp.T3(1:end, 1);
y2  = res_exp.T2(1:end, 1);
y1  = res_exp.T1(1:end, 1);
t   = res_exp.Temps(1:end, 1);

if nargin < 2 || isnan(u_echelon) 
    % Si on donne pas de u_échelon, ça signifie que l'échelon est dans le fichier excel
    u = res_exp.Commande(1:end, 1);
else
    % Si on le spécifie, on l'initialise.
    u = ones(size(t))*u_echelon;
    u(1,1) = 0;
end

if nargin >= 3 % Début donné, mais pas fin, donc on prend jusqu'à la fin
    if ~isnan(deb)
    y3  = y3(deb:end,1);
    y2  = y2(deb:end,1);
    y1  = y1(deb:end,1);
    t   = t(deb:end,1);
    u   = u(deb:end,1);
    end
elseif nargin == 4 % Début et fin précisé, on prend ce qui est demandé
    if ~isnan(fin)
    y3  = y3(1:fin, 1);
    y2  = y2(1:fin, 1);
    y1  = y1(1:fin, 1);
    t   = t(1:fin, 1);
    u   = u(1:fin, 1);
    end
end

Gain_tec = 2.2;

%% Retirer les points d'opération
y1 = y1 - y1(1, 1); % Retirer les point d'opération
y2 = y2 - y2(1, 1); % Retirer les point d'opération
y3 = y3 - y3(1, 1); % Retirer les point d'opération

%% Identification des modèles

% P -> T1
gp_1 = identify(y1, u, t, 1, 0, true, NaN);
    assignin('base', 'gp1', gp_1);
gp_1 = gp_1/Gain_tec; % On divise par gain tec pour isoler gp_1

% T1 -> T2
gp_2 = identify(y2, y1, t, 1, 0, true, 5);
    assignin('base', 'gp2', gp_2);

% T2 -> T3
gp_3 = identify(y3, y2, t, 1, 0, false, 15);
    assignin('base', 'gp3', gp_3);

% T1 -> T3
gp_4 = identify(y3, y1, t, 2, 0, false, 20);
    assignin('base', 'gp4', gp_4);

%% Identifier le controleur
% tension -> T3
gp_controleur = identify(y3, u, t, 1, 0, false, 20);

%% Avec placement de pôles
[Controleur, Ti, Td, Tf, kc] = pidfPoleCancellation(gp_controleur, 1);

%% Assigner le contrôleur dans base
assignin('base', 'Controleur', Controleur);
assignin('base', 'Ti', Ti);
assignin('base', 'Td', Td);
assignin('base', 'Tf', Tf);
assignin('base', 'kc', kc);

%% Procédé en Z
f = 1/5;
assignin('base', 'f', f);

% Puissance -> T1 (en Z)
[~, Z_T1] = s2z(gp_1, f, "zoh");
assignin('base', 'Z_T1', Z_T1);

% T1 -> T2 (en Z)
[~, Z_T2] = s2z(gp_2, f, "zoh");
assignin('base', 'Z_T2', Z_T2);

% T2 -> T3 (en Z)
[~, Z_T3] = s2z(gp_3, f, "zoh");
assignin('base', 'Z_T3', Z_T3);

% T1 -> T3 (en Z)
[~, Z_T1T3] = s2z(gp_4, f, "zoh");
assignin('base', 'Z_T1T3', Z_T1T3);

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
indent_pertub()
end