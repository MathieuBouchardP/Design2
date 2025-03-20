function sim_with_snapshots_old(json_path)
    %% Chargement des paramètres
    tic;
    params = load_json_params(json_path);
    %% Affectation des paramètres (extrait du JSON)
    TempsTotal         = params.simulation.TempsTotal;
    Lx                 = params.geometrie.Lx;
    Ly                 = params.geometrie.Ly;
    epaisseur          = params.geometrie.epaisseur;
    
    Nx                 = params.grille.Nx;
    Ny                 = params.grille.Ny;
    
    k                  = params.materiau.k;
    rho                = params.materiau.rho;
    cp                 = params.materiau.cp;
    materiau           = params.materiau.nom;
    couplage_tec       = params.puissance.couplage_TEC;
    
    alpha              = k / (rho * cp);
    
    dx                 = Lx / Nx;
    dy                 = Ly / Ny;
    dz                 = epaisseur;
    
    dt                 = (1/(4*alpha)) * (dx^2 * dy^2) / (dx^2 + dy^2);
    Nt                 = round(TempsTotal / dt);
    
    h_conv             = params.conditions_limites.h_conv;
    convection_activee = params.conditions_limites.convection_activee;
    Pin                = params.puissance.pin;
    Pin_start_time     = params.puissance.pin_start_time;
    
    Pin_end_time       = params.puissance.pin_end_time;
    if isnan(Pin_end_time)
       Pin_end_time = TempsTotal + 1;
    end
    
    Pin_loc_x_min      = fix(params.puissance.pin_loc_x_min / dx) + 1;
    Pin_loc_x_max      = fix(params.puissance.pin_loc_x_max / dx) + 1;
    Pin_loc_y_min      = fix(params.puissance.pin_loc_y_min / dy) + 1;
    Pin_loc_y_max      = fix(params.puissance.pin_loc_y_max / dy) + 1;
    
    T_piece            = params.conditions_initiales.T_piece;
    T_loc_x            = fix(params.conditions_initiales.T_loc_x / dx) + 1;
    T_loc_y            = fix(params.conditions_initiales.T_loc_y / dy) + 1;
    
    Therm1_loc_x       = params.conditions_initiales.Therm1_loc_x;
    Therm1_loc_y       = params.conditions_initiales.Therm1_loc_y;
    Therm2_loc_x       = params.conditions_initiales.Therm2_loc_x;
    Therm2_loc_y       = params.conditions_initiales.Therm2_loc_y;
    Therm3_loc_x       = params.conditions_initiales.Therm3_loc_x;
    Therm3_loc_y       = params.conditions_initiales.Therm3_loc_y;
    
    pert_loc_x_min     = fix(params.pertub.pert_loc_x_min / dx) + 1;
    pert_loc_x_max     = fix(params.pertub.pert_loc_x_max / dx) + 1;
    pert_loc_y_min     = fix(params.pertub.pert_loc_y_min / dy) + 1;
    pert_loc_y_max     = fix(params.pertub.pert_loc_y_max / dy) + 1;
    pert_pow           = params.pertub.pert_pow;
    t_pert_deb         = params.pertub.t_pert_deb;
    t_pert_fin         = params.pertub.t_pert_fin;
    
    %% Pour la convection
    aire_sides_y = dy * dz;      
    aire_sides_x = dx * dz;
    aire_top     = dx * dy;
    volume       = dx * dy * dz;
    
    Temps        = (0:Nt-1) * dt;
    x            = ((1:Nx) - 0.5) .* dx;
    y            = ((1:Ny) - 0.5) .* dy;
    [Y, X]       = meshgrid(y, x);
    
    %% Puissance déposée par l'actuateur
    P = zeros(Nx,Ny);
    nb_elts_pin = (Pin_loc_y_max - Pin_loc_y_min + 1) * (Pin_loc_x_max - Pin_loc_x_min + 1);
    P(Pin_loc_x_min:Pin_loc_x_max, Pin_loc_y_min:Pin_loc_y_max) = Pin * couplage_tec / nb_elts_pin;
    
    %% Conditions initiales
    T_piece = 273.15 + T_piece;  % conversion en Kelvin
    T       = T_piece .* ones(Nx, Ny);
    
    Therm1_loc = [(fix(Therm1_loc_x/dx) + 1), (fix(Therm1_loc_y/dy) + 1)];
    Therm2_loc = [(fix(Therm2_loc_x/dx) + 1), (fix(Therm2_loc_y/dy) + 1)];
    Therm3_loc = [(fix(Therm3_loc_x/dx) + 1), (fix(Therm3_loc_y/dy) + 1)];
    
    %% Préallocation pour diagnostics
    energy_added   = zeros(1, Nt);
    energy_loss    = zeros(1, Nt);
    thermistance1  = zeros(1, Nt);
    thermistance2  = zeros(1, Nt);
    thermistance3  = zeros(1, Nt);
    Tnew           = T;
    
    %% Précalcul des constantes
    dt_dx2 = (alpha * dt) / dx^2;
    dt_dy2 = (alpha * dt) / dy^2;
    conv_term_top     = (2 * aire_top * h_conv * dt) / (volume * rho * cp);
    conv_term_sides_y = (aire_sides_y * h_conv * dt) / (volume * rho * cp);
    conv_term_sides_x = (aire_sides_x * h_conv * dt) / (volume * rho * cp);
    power_term        = dt / (rho * cp * volume);
    deposited         = 0;
    
    %% Configuration de la figure
    f1 = figure(1);
    sgtitle(strcat("Distribution de température sur une plaque d'", materiau));
    set(gcf, 'Units', 'normalized', 'OuterPosition', [0 0 1 1]);
    set(f1, 'Color', 'w');
    rotate3d on;
    
    % Graphique 1 : Surface
    subplot(131)
    f1_surf = meshc(1000*X, 1000*Y, T-273.15);
    shading interp;
    xlabel('Longueur de la plaque (mm)')
    ylabel('Largeur de la plaque (mm)')
    zlabel('Température en °C')
    title(['Température à t = ', num2str(0), ' s'],'FontSize',16);
    timeText = title(['Temps : ' num2str(0, '%.2f') ' s'], 'FontSize', 16);
    colorbar;
    colormap jet; 
    clim([10 50]);
    grid on;
    view(3); 
    zticks(floor(0):1:ceil(100));
    axis manual;
    pbaspect([Lx Ly min([Lx Ly])]);
    axis 'auto z'
    
    % Graphique 2 : Thermistances
    subplot(132);
    hold on;
    f2_t1 = plot(Temps(1:2), thermistance1(1:2)-273.15, 'r', 'DisplayName', 'Thermistance 1');
    f2_t2 = plot(Temps(1:2), thermistance2(1:2)-273.15, 'g', 'DisplayName', 'Thermistance 2');
    f2_t3 = plot(Temps(1:2), thermistance3(1:2)-273.15, 'b', 'DisplayName', 'Thermistance 3');
    grid on; 
    ax = gca;
    ax.FontSize = 16;
    xlabel('Temps [s]', 'FontSize', 16);
    ylabel('Température [°C]', 'FontSize', 16);
    title('Température aux thermistances', 'FontSize', 16);
    legend('show', 'FontSize', 14, 'Location', 'best');
    
    % Graphique 3 : Énergie
    subplot(133);
    hold on;
    f3_add = plot(Temps(1:2), energy_added(1:2));
    f3_loss = plot(Temps(1:2), energy_loss(1:2));
    xlabel('Temps [s]', 'FontSize', 16);
    ylabel('Énergie dans l''itération', 'FontSize', 16);
    legend('Energie déposée', 'Energie dissipée par convection', 'FontSize', 16, 'Location', 'southeast');
    grid on;
    
    %% Préallocation des snapshots pour sauvegarde (facultatif)
    snapshotPeriod = 750;  % Sauvegarder toutes les 10 itérations (modifiez selon vos besoins)
    nbSnapshots = floor(Nt / snapshotPeriod) + 1;
    snapshots = zeros(Nx, Ny, nbSnapshots);
    snapshotIndex = 1;
    snapshots(:, :, snapshotIndex) = T;
    
    %% Création du DataQueue pour la communication
    dq = parallel.pool.DataQueue;
    % Le callback updateSnapshot reçoit une structure contenant le snapshot et l'indice (ou temps)
    afterEach(dq, @updateSnapshot);
    
    %% Lancement de la simulation en arrière-plan (parfeval)
    f = parfeval(@simulationWithSnapshots, 0, T, Nt, snapshotPeriod, dq);
    
    %% Boucle de traitement léger pour permettre aux callbacks de s'exécuter
    while ~strcmp(f.State, 'finished')
        drawnow;   % Permet de traiter les callbacks et rafraîchir la figure
        pause(0.05);
    end
    wait(f);  % On s'assure que la simulation est terminée
    
    runtime = toc;
    fprintf('Simulation complète en %.6f secondes\n', runtime);
    
    %% --- Fonction de simulation (exécutée en parallèle) ---
    function simulationWithSnapshots(T0, Nt_local, snapshotPeriod_local, dq_local)
        % La simulation s'exécute ici en arrière-plan.
        T_sim = T0;
        localSnapshotIndex = 1;
        % Préallocation locale (optionnelle) pour sauvegarder les snapshots
        localSnapshots = zeros(Nx, Ny, floor(Nt_local/snapshotPeriod_local)+1);
        localSnapshots(:, :, localSnapshotIndex) = T_sim;
        
        for t = 1:Nt_local
            % Calcul de Tnew en suivant le schéma de conduction/convection (exemple simplifié)
            Tnew = T_sim;
            % Conduction interne
            Tnew(2:Nx-1, 2:Ny-1) = T_sim(2:Nx-1, 2:Ny-1) ...
                + dt_dx2 * (T_sim(1:Nx-2, 2:Ny-1) - 2*T_sim(2:Nx-1, 2:Ny-1) + T_sim(3:Nx, 2:Ny-1)) ...
                + dt_dy2 * (T_sim(2:Nx-1, 1:Ny-2) - 2*T_sim(2:Nx-1, 2:Ny-1) + T_sim(2:Nx-1, 3:Ny));
            % Conduction aux bords (exemples simplifiés)
            Tnew(1, 2:Ny-1) = T_sim(1, 2:Ny-1) ...
                + dt_dx2 * (T_sim(2, 2:Ny-1) - T_sim(1, 2:Ny-1)) ...
                + dt_dy2 * (T_sim(1, 1:Ny-2) - 2*T_sim(1, 2:Ny-1) + T_sim(1, 3:Ny));
            Tnew(Nx, 2:Ny-1) = T_sim(Nx, 2:Ny-1) ...
                + dt_dx2 * (T_sim(Nx-1, 2:Ny-1) - T_sim(Nx, 2:Ny-1)) ...
                + dt_dy2 * (T_sim(Nx, 1:Ny-2) - 2*T_sim(Nx, 2:Ny-1) + T_sim(Nx, 3:Ny));
            Tnew(2:Nx-1, 1) = T_sim(2:Nx-1, 1) ...
                + dt_dy2 * (T_sim(1:Nx-2, 1) - 2*T_sim(2:Nx-1, 1) + T_sim(3:Nx, 1)) ...
                + dt_dx2 * (T_sim(2:Nx-1, 2) - T_sim(2:Nx-1, 1));
            Tnew(2:Nx-1, Ny) = T_sim(2:Nx-1, Ny) ...
                + dt_dy2 * (T_sim(1:Nx-2, Ny) - 2*T_sim(2:Nx-1, Ny) + T_sim(3:Nx, Ny)) ...
                + dt_dx2 * (T_sim(2:Nx-1, Ny-1) - T_sim(2:Nx-1, Ny));
            % Conduction aux coins
            Tnew(1, 1) = T_sim(1, 1) + dt_dx2 * (T_sim(2, 1) - T_sim(1, 1)) + dt_dy2 * (T_sim(1, 2) - T_sim(1, 1));
            Tnew(1, Ny) = T_sim(1, Ny) + dt_dx2 * (T_sim(2, Ny) - T_sim(1, Ny)) + dt_dy2 * (T_sim(1, Ny-1) - T_sim(1, Ny));
            Tnew(Nx, 1) = T_sim(Nx, 1) + dt_dx2 * (T_sim(Nx-1, 1) - T_sim(Nx, 1)) + dt_dy2 * (T_sim(Nx, 2) - T_sim(Nx, 1));
            Tnew(Nx, Ny) = T_sim(Nx, Ny) + dt_dx2 * (T_sim(Nx-1, Ny) - T_sim(Nx, Ny)) + dt_dy2 * (T_sim(Nx, Ny-1) - T_sim(Nx, Ny));
            % Convection
            Tnew(1:Nx,1:Ny) = Tnew(1:Nx,1:Ny) - conv_term_top .* (T_sim(1:Nx,1:Ny) - T_piece);
            Tnew(1, :)      = Tnew(1, :) - conv_term_sides_y .* (T_sim(1, :) - T_piece);
            Tnew(Nx, :)     = Tnew(Nx, :) - conv_term_sides_y .* (T_sim(Nx, :) - T_piece);
            Tnew(:, 1)      = Tnew(:, 1) - conv_term_sides_x .* (T_sim(:, 1) - T_piece);
            Tnew(:, Ny)     = Tnew(:, Ny) - conv_term_sides_x .* (T_sim(:, Ny) - T_piece);
            
            % Ajout des perturbations (si applicable)
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
            
            % Mise à jour de la température et des mesures
            T_sim = Tnew;
            thermistance1(t) = T_sim(Therm1_loc(1), Therm1_loc(2));
            thermistance2(t) = T_sim(Therm2_loc(1), Therm2_loc(2));
            thermistance3(t) = T_sim(Therm3_loc(1), Therm3_loc(2));
            
            % Calcul des bilans d'énergie (pour suivi)
            energy_added(t) = sum(P(:)) * dt;
            energy_loss_sides_ligne_1 = (h_conv * aire_sides_y * dt) * sum(T_sim(1, 1:Ny) - T_piece);
            energy_loss_sides_ligne_f = (h_conv * aire_sides_y * dt) * sum(T_sim(Nx, 1:Ny) - T_piece);
            energy_loss_top_down = (h_conv * 2 * aire_top * dt) * sum(sum(T_sim(1:Nx, 1:Ny) - T_piece));
            energy_loss_sides_colonne_1 = (h_conv * aire_sides_x * dt) * sum(T_sim(1:Nx, 1) - T_piece);
            energy_loss_sides_colonne_f = (h_conv * aire_sides_x * dt) * sum(T_sim(1:Nx, Ny) - T_piece);
            energy_loss(t) = energy_loss_sides_ligne_1 + energy_loss_sides_ligne_f + energy_loss_top_down + energy_loss_sides_colonne_1 + energy_loss_sides_colonne_f;
            
            % Sauvegarde et envoi du snapshot toutes les snapshotPeriod_local itérations (et à t==1)
            if mod(t, snapshotPeriod_local) == 0 || t == 1
                localSnapshotIndex = localSnapshotIndex + 1;
                localSnapshots(:, :, localSnapshotIndex) = T_sim;
                % Envoi d'une structure contenant le snapshot et le temps (itération) courant
                dataOut.snapshot = T_sim;
                dataOut.t = t;
                send(dq_local, dataOut);
            end
        end
    end

    %% --- Callback du DataQueue ---
    function updateSnapshot(data)
        % Mise à jour de l'affichage dès qu'une nouvelle snapshot est reçue.
        % On accède ici aux variables f1_surf, timeText et dt (définies dans le parent)
        latestSnapshot = data.snapshot;
        currentTime = data.t * dt;
        set(f1_surf, 'ZData', latestSnapshot - 273.15);
        set(timeText, 'String', ['Temps : ' num2str(currentTime, '%.2f') ' s']);
        drawnow;  % Forcer le rafraîchissement de la figure
    end
end

%% --- Fonction de chargement du JSON ---
function params = load_json_params(filename)
    fid = fopen(filename, 'r');
    raw = fread(fid, inf, 'uint8=>char')';
    fclose(fid);
    params = jsondecode(raw);
end
