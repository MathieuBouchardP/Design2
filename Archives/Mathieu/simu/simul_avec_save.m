function simul_avec_save(json_path)
    tic;
    params = load_json_params(json_path);
        % Assigner les valeurs aux variables dans l'espace de travail
    TempsTotal       = params.simulation.TempsTotal;
    Lx               = params.geometrie.Lx;
    Ly               = params.geometrie.Ly;
    epaisseur        = params.geometrie.epaisseur;
    
    Nx               = params.grille.Nx;
    Ny               = params.grille.Ny;
    
    k                = params.materiau.k;
    rho              = params.materiau.rho;
    cp               = params.materiau.cp;
    materiau         = params.materiau.nom;
    couplage_tec     = params.puissance.couplage_TEC;
    
    alpha            = k / (rho * cp);
    
    dx               = Lx / Nx;
    dy               = Ly / Ny;
    dz               = epaisseur;
    
    dt               = (1/(4*alpha)) * (dx^2 * dy^2) / (dx^2 + dy^2);
    Nt               = round(TempsTotal / dt);
    
    h_conv           = params.conditions_limites.h_conv;
    convection_activee = params.conditions_limites.convection_activee;
    Pin              = params.puissance.pin;
    Pin_start_time   = params.puissance.pin_start_time;
    
    Pin_end_time     = params.puissance.pin_end_time;
    if isnan(Pin_end_time)
       Pin_end_time = TempsTotal + 1;
    end
    
    Pin_loc_x_min    = fix(params.puissance.pin_loc_x_min / dx) + 1;
    Pin_loc_x_max    = fix(params.puissance.pin_loc_x_max / dx) + 1;
    Pin_loc_y_min    = fix(params.puissance.pin_loc_y_min / dy) + 1;
    Pin_loc_y_max    = fix(params.puissance.pin_loc_y_max / dy) + 1;
    
    T_piece          = params.conditions_initiales.T_piece;
    T_loc_x          = fix(params.conditions_initiales.T_loc_x / dx) + 1;
    T_loc_y          = fix(params.conditions_initiales.T_loc_y / dy) + 1;
    
    Therm1_loc_x     = params.conditions_initiales.Therm1_loc_x;
    Therm1_loc_y     = params.conditions_initiales.Therm1_loc_y;
    Therm2_loc_x     = params.conditions_initiales.Therm2_loc_x;
    Therm2_loc_y     = params.conditions_initiales.Therm2_loc_y;
    Therm3_loc_x     = params.conditions_initiales.Therm3_loc_x;
    Therm3_loc_y     = params.conditions_initiales.Therm3_loc_y;
    
    pert_loc_x_min   = fix(params.pertub.pert_loc_x_min / dx) + 1;
    pert_loc_x_max   = fix(params.pertub.pert_loc_x_max / dx) + 1;
    pert_loc_y_min   = fix(params.pertub.pert_loc_y_min / dy) + 1;
    pert_loc_y_max   = fix(params.pertub.pert_loc_y_max / dy) + 1;
    pert_pow         = params.pertub.pert_pow;
    t_pert_deb       = params.pertub.t_pert_deb;
    t_pert_fin       = params.pertub.t_pert_fin;
    
%% Pour la convection
    aire_sides_y = dy*dz;      % aire des côtés sur la largeur y
    aire_sides_x = dx*dz;      % aire des côtés sur x
    aire_top = dx*dy;          % aire du dessus et du dessous

    volume = dx*dy*dz;         % Volume d'un élément

    Temps = (0:Nt-1)*dt;       % Vecteur de temps pour la simulation
    x = ((1:Nx) - 0.5).*dx;    % Vecteur de coordonnées en x
    y = ((1:Ny) - 0.5).*dy;    % Vecteur de coordonnées en y
    [Y, X] = meshgrid(y, x);   % Génération des matrices coordonnées 2D
    
    %% Puissance déposée par l'actuateur
    P = zeros(Nx,Ny);          % Matrice de puissance à ajouter à chaque élément

    nb_elts_pin = (Pin_loc_y_max - Pin_loc_y_min + 1) * (Pin_loc_x_max - Pin_loc_x_min + 1);
    P(Pin_loc_x_min:Pin_loc_x_max, Pin_loc_y_min:Pin_loc_y_max) = Pin * couplage_tec / nb_elts_pin;
    
    %% Conditions initiales
    T_piece = 273.15 + T_piece;    % conversion Température pièce en [K]
    T = T_piece .* ones(Nx, Ny);   % Température de tous les éléments
     
    Therm1_loc = [(fix(Therm1_loc_x/dx) + 1), (fix(Therm1_loc_y/dy) + 1)]; % Thermistance 1
    Therm2_loc = [(fix(Therm2_loc_x/dx) + 1), (fix(Therm2_loc_y/dy) + 1)]; % Thermistance 2
    Therm3_loc = [(fix(Therm3_loc_x/dx) + 1), (fix(Therm3_loc_y/dy) + 1)]; % Thermistance 3

   %% Préallocation des vecteurs utilisés dans la boucle
    energy_added = zeros(1, Nt);
    energy_loss  = zeros(1, Nt);

    thermistance_1 = T_piece .* ones(1, Nt);
    thermistance1  = zeros(1, Nt);
    thermistance_2 = T_piece .* ones(1, Nt);
    thermistance2  = zeros(1, Nt);
    thermistance_3 = T_piece .* ones(1, Nt);
    thermistance3  = zeros(1, Nt);
    Tnew = T;
    
    % Préallocation pour la sauvegarde des données toutes les 50 itérations
    num_save = floor(Nt / 50);
    T_record = zeros(Nx, Ny, num_save);
    therm_record = zeros(3, num_save);
    energy_record_added = zeros(1, num_save);
    energy_record_loss = zeros(1, num_save);
    time_record = zeros(1, num_save);
    record_idx = 0;
    
%% Précalcul des constantes
    dt_dx2 = (alpha * dt) / dx^2;
    dt_dy2 = (alpha * dt) / dy^2;
    conv_term_top = (2 * aire_top * h_conv * dt) / (volume * rho * cp);
    conv_term_sides_y = (aire_sides_y * h_conv * dt) / (volume * rho * cp);
    conv_term_sides_x = (aire_sides_x * h_conv * dt) / (volume * rho * cp);
    power_term = dt / (rho * cp * volume);
    deposited = 0;
    
% Critère de fin
    t_end = Pin_end_time * Nt / TempsTotal;
    is_end_contition_met = false;
    
%% Boucle principale de simulation
for t = 1:Nt
    Tnew = T;
    if t > t_end
        P = 0;
        is_end_contition_met = true;
    end
    
    % Conduction interne (milieu)
    Tnew(2:Nx-1, 2:Ny-1) = T(2:Nx-1, 2:Ny-1) ...
        + dt_dx2 * (T(1:Nx-2, 2:Ny-1) - 2*T(2:Nx-1, 2:Ny-1) + T(3:Nx, 2:Ny-1)) ...
        + dt_dy2 * (T(2:Nx-1, 1:Ny-2) - 2*T(2:Nx-1, 2:Ny-1) + T(2:Nx-1, 3:Ny));
    
    % Conduction aux bords
    Tnew(1, 2:Ny-1) = T(1, 2:Ny-1) ...
        + dt_dx2 * (T(2, 2:Ny-1) - T(1, 2:Ny-1)) ...
        + dt_dy2 * (T(1, 1:Ny-2) - 2*T(1, 2:Ny-1) + T(1, 3:Ny));
    
    Tnew(Nx, 2:Ny-1) = T(Nx, 2:Ny-1) ...
        + dt_dx2 * (T(Nx-1, 2:Ny-1) - T(Nx, 2:Ny-1)) ...
        + dt_dy2 * (T(Nx, 1:Ny-2) - 2*T(Nx, 2:Ny-1) + T(Nx, 3:Ny));
    
    Tnew(2:Nx-1, 1) = T(2:Nx-1, 1) ...
        + dt_dy2 * (T(1:Nx-2, 1) - 2*T(2:Nx-1, 1) + T(3:Nx, 1)) ...
        + dt_dx2 * (T(2:Nx-1, 2) - T(2:Nx-1, 1));
    
    Tnew(2:Nx-1, Ny) = T(2:Nx-1, Ny) ...
        + dt_dy2 * (T(1:Nx-2, Ny) - 2*T(2:Nx-1, Ny) + T(3:Nx, Ny)) ...
        + dt_dx2 * (T(2:Nx-1, Ny-1) - T(2:Nx-1, Ny));
    
    % Conduction aux coins
    Tnew(1, 1) = T(1, 1) ...
        + dt_dx2 * (T(2, 1) - T(1, 1)) ...
        + dt_dy2 * (T(1, 2) - T(1, 1));
    
    Tnew(1, Ny) = T(1, Ny) ...
        + dt_dx2 * (T(2, Ny) - T(1, Ny)) ...
        + dt_dy2 * (T(1, Ny-1) - T(1, Ny));
    
    Tnew(Nx, 1) = T(Nx, 1) ...
        + dt_dx2 * (T(Nx-1, 1) - T(Nx, 1)) ...
        + dt_dy2 * (T(Nx, 2) - T(Nx, 1));
    
    Tnew(Nx, Ny) = T(Nx, Ny) ...
        + dt_dx2 * (T(Nx-1, Ny) - T(Nx, Ny)) ...
        + dt_dy2 * (T(Nx, Ny-1) - T(Nx, Ny));
    
    % Convection
    Tnew(1:Nx,1:Ny) = Tnew(1:Nx,1:Ny) - conv_term_top .* (T(1:Nx,1:Ny) - T_piece);
    Tnew(1, :) = Tnew(1, :) - conv_term_sides_y .* (T(1, :) - T_piece);
    Tnew(Nx, :) = Tnew(Nx, :) - conv_term_sides_y .* (T(Nx, :) - T_piece);
    Tnew(:, 1) = Tnew(:, 1) - conv_term_sides_x .* (T(:, 1) - T_piece);
    Tnew(:, Ny) = Tnew(:, Ny) - conv_term_sides_x .* (T(:, Ny) - T_piece);
    
    % Ajout des perturbations
    if pert_pow ~= 0 && deposited == 0 
        if (round(t_pert_deb/dt) == t)
            nb_elts_pert = (pert_loc_y_max - pert_loc_y_min + 1) * (pert_loc_x_max - pert_loc_x_min + 1);
            ajout = (pert_pow / nb_elts_pert);
            P(pert_loc_x_min:pert_loc_x_max, pert_loc_y_min:pert_loc_y_max) = ...
                P(pert_loc_x_min:pert_loc_x_max, pert_loc_y_min:pert_loc_y_max) + ajout;
            deposited = deposited + 1;
        end
    end
  
    if pert_pow ~= 0 && deposited == 1
        if (t == round(t_pert_fin/dt))
            P(pert_loc_x_min:pert_loc_x_max, pert_loc_y_min:pert_loc_y_max) = ...
                P(pert_loc_x_min:pert_loc_x_max, pert_loc_y_min:pert_loc_y_max) - ajout;
            deposited = deposited + 1;
        end 
    end
    % Ajout de la puissance
    Tnew = Tnew + power_term .* P;
    
    % Mise à jour de la température et enregistrement des valeurs
    T = Tnew;
    thermistance1(t) = T(Therm1_loc(1), Therm1_loc(2));
    thermistance2(t) = T(Therm2_loc(1), Therm2_loc(2));
    thermistance3(t) = T(Therm3_loc(1), Therm3_loc(2));
    
    % Calcul des bilans d'énergie
    % energy_added(t) = sum(P(:)) * dt;
    
    % energy_loss_sides_ligne_1 = (h_conv * aire_sides_y * dt) * sum(T(1, 1:Ny) - T_piece);
    % energy_loss_sides_ligne_f = (h_conv * aire_sides_y * dt) * sum(T(Nx, 1:Ny) - T_piece);
    % energy_loss_top_down = (h_conv * 2 * aire_top * dt) * sum(sum(T(1:Nx, 1:Ny) - T_piece));
    % energy_loss_sides_colonne_1 = (h_conv * aire_sides_x * dt) * sum(T(1:Nx, 1) - T_piece);
    % energy_loss_sides_colonne_f = (h_conv * aire_sides_x * dt) * sum(T(1:Nx, Ny) - T_piece);
    
    % energy_loss(t) = energy_loss_sides_ligne_1 + energy_loss_sides_ligne_f + ...
    %                  energy_loss_top_down + energy_loss_sides_colonne_1 + energy_loss_sides_colonne_f;
    
    % Sauvegarde des données toutes les 50 itérations
    if mod(t,750) == 0
        record_idx = record_idx + 1;
        T_record(:,:,record_idx) = T;
        therm_record(:,record_idx) = [thermistance1(t); thermistance2(t); thermistance3(t)];
        % energy_record_added(record_idx) = energy_added(t);
        % energy_record_loss(record_idx) = energy_loss(t);
        time_record(record_idx) = t * dt;
    end
    
    % Critère de fin
    % if abs((energy_added(t) - energy_loss(t)) / (energy_added(t)+1e-9)) < 1e-3 && is_end_contition_met
    %     runtime = toc;
    %     fprintf('Temps d''exécution : %.6f secondes\n', runtime);
    %     break;
    % end
end

% Ajustement de la taille des enregistrements en cas d'arrêt prématuré
T_record = T_record(:,:,1:record_idx);
therm_record = therm_record(:,1:record_idx);
% energy_record_added = energy_record_added(1:record_idx);
% energy_record_loss = energy_record_loss(1:record_idx);
time_record = time_record(1:record_idx);
runtime = toc;
fprintf('Temps d''exécution : %.6f secondes\n', runtime);


% Sauvegarde des données enregistrées dans un fichier .mat avec un nom unique
filename = sprintf('simulation_data_%s.mat', datestr(now, 'yyyymmdd_HHMMSS'));
save(filename, 'T_record', 'therm_record', 'energy_record_added', 'energy_record_loss', 'time_record', '-v7.3');
fprintf('Données sauvegardées dans le fichier : %s\n', filename);

end

function params = load_json_params(filename)
    % Lire le contenu du fichier JSON
    fid = fopen(filename, 'r');
    raw = fread(fid, inf, 'uint8=>char')';
    fclose(fid);
    
    % Décoder le JSON en structure MATLAB
    params = jsondecode(raw);
end
