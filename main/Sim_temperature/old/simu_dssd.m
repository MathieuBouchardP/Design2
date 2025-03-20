function simu_dssd(json_path)
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
    T_init             = params.conditions_initiales.T_init;
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
    
    %% Puissance déposée par l'actuateur
    P = zeros(Nx,Ny);
    nb_elts_pin = (Pin_loc_y_max - Pin_loc_y_min + 1) * (Pin_loc_x_max - Pin_loc_x_min + 1);
    P(Pin_loc_x_min:Pin_loc_x_max, Pin_loc_y_min:Pin_loc_y_max) = Pin * couplage_tec / nb_elts_pin;
    
    %% Conditions initiales
    T_piece = 273.15 + T_piece;  % conversion en Kelvin
    T_init = 273.15 + T_init;  % conversion en Kelvin
    T       = T_init .* ones(Nx, Ny);

    % Initialiser les thermistances avec la température initiale réelle
    Therm1_loc = [(fix(Therm1_loc_x/dx) + 1), (fix(Therm1_loc_y/dy) + 1)];
    Therm2_loc = [(fix(Therm2_loc_x/dx) + 1), (fix(Therm2_loc_y/dy) + 1)];
    Therm3_loc = [(fix(Therm3_loc_x/dx) + 1), (fix(Therm3_loc_y/dy) + 1)];
    
    %% Précalcul des constantes
    dt_dx2 = (alpha * dt) / dx^2;
    dt_dy2 = (alpha * dt) / dy^2;
    conv_term_top     = (2 * aire_top * h_conv * dt) / (volume * rho * cp);
    conv_term_sides_y = (aire_sides_y * h_conv * dt) / (volume * rho * cp);
    conv_term_sides_x = (aire_sides_x * h_conv * dt) / (volume * rho * cp);
    power_term        = dt / (rho * cp * volume);
    deposited         = 0;
    %% Lancement de la simulation en arrière-plan via parfeval
    
    stop(updateTimer);
    delete(updateTimer);
    
    runtime = toc;
    fprintf('Simulation complète en %.6f secondes\n', runtime);
    
    %% Nettoyage du pool à la fin
    delete(pool);
    
    %% --- Fonction de simulation (exécutée en parallèle) ---
    function simulationWith_snapshots(T0, Nt_local, snapshotPeriod_local, dq_local, )
        T = T0;
        % Tableaux locaux pour mesurer (pour le worker)
        localTherm1 = zeros(1, Nt_local);
        localTherm2 = zeros(1, Nt_local);
        localTherm3 = zeros(1, Nt_local);
        localEnergyAdded = zeros(1, Nt_local);
        localEnergyLoss = zeros(1, Nt_local);
        

        function T_update = sim_n_iter(n_iter, init)
            if init
                persistent  
            end

    









        for t = 1:Nt_local
            function simul_n_iteration(nbr_ite, T)
            end
            % Calcul de Tnew par schéma numérique (inchangé)
            Tnew = T;
            Tnew(2:Nx-1, 2:Ny-1) = T(2:Nx-1, 2:Ny-1) + dt_dx2*(T(1:Nx-2, 2:Ny-1)-2*T(2:Nx-1, 2:Ny-1)+T(3:Nx, 2:Ny-1)) + dt_dy2*(T(2:Nx-1, 1:Ny-2)-2*T(2:Nx-1, 2:Ny-1)+T(2:Nx-1, 3:Ny));
            Tnew(1, 2:Ny-1) = T(1, 2:Ny-1) + dt_dx2*(T(2, 2:Ny-1)-T(1, 2:Ny-1)) + dt_dy2*(T(1, 1:Ny-2)-2*T(1, 2:Ny-1)+T(1, 3:Ny));
            Tnew(Nx, 2:Ny-1) = T(Nx, 2:Ny-1) + dt_dx2*(T(Nx-1, 2:Ny-1)-T(Nx, 2:Ny-1)) + dt_dy2*(T(Nx, 1:Ny-2)-2*T(Nx, 2:Ny-1)+T(Nx, 3:Ny));
            Tnew(2:Nx-1, 1) = T(2:Nx-1, 1) + dt_dy2*(T(1:Nx-2, 1)-2*T(2:Nx-1, 1)+T(3:Nx, 1)) + dt_dx2*(T(2:Nx-1, 2)-T(2:Nx-1, 1));
            Tnew(2:Nx-1, Ny) = T(2:Nx-1, Ny) + dt_dy2*(T(1:Nx-2, Ny)-2*T(2:Nx-1, Ny)+T(3:Nx, Ny)) + dt_dx2*(T(2:Nx-1, Ny-1)-T(2:Nx-1, Ny));
            Tnew(1, 1) = T(1, 1) + dt_dx2*(T(2, 1)-T(1, 1)) + dt_dy2*(T(1, 2)-T(1, 1));
            Tnew(1, Ny) = T(1, Ny) + dt_dx2*(T(2, Ny)-T(1, Ny)) + dt_dy2*(T(1, Ny-1)-T(1, Ny));
            Tnew(Nx, 1) = T(Nx, 1) + dt_dx2*(T(Nx-1, 1)-T( Nx, 1)) + dt_dy2*(T(Nx, 2)-T(Nx, 1));
            Tnew(Nx, Ny) = T(Nx, Ny) + dt_dx2*(T(Nx-1, Ny)-T(Nx, Ny)) + dt_dy2*(T(Nx, Ny-1)-T(Nx, Ny));
            Tnew(1:Nx,1:Ny) = Tnew(1:Nx,1:Ny) - conv_term_top .* (T(1:Nx,1:Ny)-T_piece);
            Tnew(1, :) = Tnew(1, :) - conv_term_sides_y .* (T(1, :) - T_piece);
            Tnew(Nx, :) = Tnew(Nx, :) - conv_term_sides_y .* (T(Nx, :) - T_piece);
            Tnew(:, 1) = Tnew(:, 1) - conv_term_sides_x .* (T(:, 1) - T_piece);
            Tnew(:, Ny) = Tnew(:, Ny) - conv_term_sides_x .* (T(:, Ny) - T_piece);
            
            if pert_pow ~= 0 && deposited == 0 
                if (round(t_pert_deb/dt) == t)
                    nb_elts_pert = (pert_loc_y_max-pert_loc_y_min+1)*(pert_loc_x_max-pert_loc_x_min+1);
                    ajout = (pert_pow/nb_elts_pert);
                    P(pert_loc_x_min:pert_loc_x_max, pert_loc_y_min:pert_loc_y_max) = P(pert_loc_x_min:pert_loc_x_max, pert_loc_y_min:pert_loc_y_max) + ajout;
                    deposited = deposited + 1;
                end
            end
            if pert_pow ~= 0 && deposited == 1
                if (t == round(t_pert_fin/dt))
                    P(pert_loc_x_min:pert_loc_x_max, pert_loc_y_min:pert_loc_y_max) = P(pert_loc_x_min:pert_loc_x_max, pert_loc_y_min:pert_loc_y_max) - ajout;
                    deposited = deposited + 1;
                end 
            end
            Tnew = Tnew + power_term .* P;

            energy_added(t) = sum(P(:)) * dt;
            energy_loss_sides_ligne_1 = (h_conv*aire_sides_y*dt)* sum(T(1, 1:Ny)-T_piece);
            energy_loss_sides_ligne_f = (h_conv*aire_sides_y*dt)*sum(T(Nx, 1:Ny)-T_piece);
            energy_loss_top_down = (h_conv*2*aire_top*dt)*sum(sum(T(1:Nx, 1:Ny)-T_piece));
            energy_loss_sides_colonne_1 = (h_conv*aire_sides_x*dt)* sum(T(1:Nx, 1)-T_piece);
            energy_loss_sides_colonne_f = (h_conv*aire_sides_x*dt)* sum(T(1:Nx, Ny)-T_piece);
            energy_loss(t) = energy_loss_sides_ligne_1 + energy_loss_sides_ligne_f + energy_loss_top_down + energy_loss_sides_colonne_1 + energy_loss_sides_colonne_f;
            
            T = Tnew;
            localTherm1(t) = T(Therm1_loc(1), Therm1_loc(2));
            localTherm2(t) = T(Therm2_loc(1), Therm2_loc(2));
            localTherm3(t) = T(Therm3_loc(1), Therm3_loc(2));
            localEnergyAdded(t) = sum(P(:)) * dt;
            localEnergyLoss(t) = energy_loss_sides_ligne_1 + energy_loss_sides_ligne_f + energy_loss_top_down + energy_loss_sides_colonne_1 + energy_loss_sides_colonne_f;
            
            if mod(t, snapshotPeriod_local) == 0 || t == 1
                dataOut.snapshot = T;
                dataOut.t = t;
                dataOut.therm1 = localTherm1(t);
                dataOut.therm2 = localTherm2(t);
                dataOut.therm3 = localTherm3(t);
                dataOut.energyAdded = localEnergyAdded(t);
                dataOut.energyLoss = localEnergyLoss(t);
                current_state = appConstant.Value.state;
                dataOut.state = current_state;
                send(dq_local, dataOut);
            end
        end
    end

    %% --- Callback du DataQueue ---
    function updateSnapshot(data)
        if data.state == 1
            return;
        end
    % Mise à jour du graphique 1 (surface)
    set(f1_surf, 'ZData', data.snapshot - 273.15);
    currentTime = data.t * dt;
    set(timeText, 'String', sprintf('Temps : %.2f s', currentTime));
    
    % Déclaration persistante pour stocker les données secondaires
    persistent timeData therm1Data therm2Data therm3Data energyAddedData energyLossData updateCounter
    
    % Initialisation ou mise à jour des données persistantes
    if isempty(timeData)
        timeData = currentTime;
        therm1Data = data.therm1;
        therm2Data = data.therm2;
        therm3Data = data.therm3;
        energyAddedData = data.energyAdded;
        energyLossData = data.energyLoss;
        updateCounter = 1;
    else
        timeData(end+1) = currentTime;
        if ~isempty(data.therm1)
            therm1Data(end+1) = data.therm1;
            therm2Data(end+1) = data.therm2;
            therm3Data(end+1) = data.therm3;
            energyAddedData(end+1) = data.energyAdded;
            energyLossData(end+1) = data.energyLoss;
        end
        updateCounter = updateCounter + 1;
    end
    
    % Mise à jour des graphiques secondaires une fois sur deux
    if mod(updateCounter, 2) == 0 && ~isempty(therm1Data)
        set(f2_t1, 'XData', timeData, 'YData', therm1Data - 273.15);
        set(f2_t2, 'XData', timeData, 'YData', therm2Data - 273.15);
        set(f2_t3, 'XData', timeData, 'YData', therm3Data - 273.15);
        set(f3_add, 'XData', timeData, 'YData', energyAddedData);
        set(f3_loss, 'XData', timeData, 'YData', energyLossData);
    end
    drawnow expose;
end
end

%% --- Fonction de chargement du JSON ---
function params = load_json_params(filename)
    fid = fopen(filename, 'r');
    raw = fread(fid, inf, 'uint8=>char')';
    fclose(fid);
    params = jsondecode(raw);
end
