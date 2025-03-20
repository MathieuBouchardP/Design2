function Copy_of_sim_with_snapshots(json_path)
    %% Réinitialisation des variables persistantes et nettoyage du pool
    clear updateSnapshot; 
    pool = gcp('nocreate');
    if ~isempty(pool)
        delete(pool);
    end
    pool = parpool('local', 1);  % Utilisation d'un seul worker pour cette simulation

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
    
    % Initialiser les thermistances avec la température initiale réelle
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
    %set(gcf, 'Renderer', 'opengl');
    set(f1, 'Interruptible', 'on', 'BusyAction', 'cancel');
    sgtitle(strcat("Distribution de température sur une plaque d'", materiau));
    set(gcf, 'Units', 'normalized', 'OuterPosition', [0 0 1 1]);
    set(f1, 'Color', 'w');
    rotate3d on;
    
    % Graphique 1 : Surface de la plaque
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
    
    % Graphique 2 : Température aux thermistances
    subplot(132);
    hold on;
    initTherm1 = T(Therm1_loc(1), Therm1_loc(2));
    initTherm2 = T(Therm2_loc(1), Therm2_loc(2));
    initTherm3 = T(Therm3_loc(1), Therm3_loc(2));
    f2_t1 = plot(0, initTherm1 - 273.15, 'r', 'DisplayName', 'Thermistance 1');
    f2_t2 = plot(0, initTherm2 - 273.15, 'g', 'DisplayName', 'Thermistance 2');
    f2_t3 = plot(0, initTherm3 - 273.15, 'b', 'DisplayName', 'Thermistance 3');
    grid on; 
    ax = gca;
    ax.FontSize = 16;
    xlabel('Temps [s]', 'FontSize', 16);
    ylabel('Température [°C]', 'FontSize', 16);
    title('Température aux thermistances', 'FontSize', 16);
    legend('show', 'FontSize', 14, 'Location', 'best');
    
    % Graphique 3 : Bilans énergétiques
    subplot(133);
    hold on;
    f3_add = plot(0, 0, 'LineWidth', 1.5, 'DisplayName', 'Energie déposée');
    f3_loss = plot(0, 0, 'LineWidth', 1.5, 'DisplayName', 'Energie dissipée');
    xlabel('Temps [s]', 'FontSize', 16);
    ylabel('Énergie', 'FontSize', 16);
    legend('show', 'FontSize', 16, 'Location', 'southeast');
    grid on;
    
    %% Création d'un timer pour forcer le rafraîchissement continu
    updateTimer = timer('ExecutionMode', 'fixedRate', 'Period', 0.05, 'TimerFcn', @(~,~) drawnow);
    start(updateTimer);
    
    %% Création du DataQueue pour communication depuis le worker
    dq = parallel.pool.DataQueue;
    afterEach(dq, @updateSnapshot);
    
    %% Lancement de la simulation en arrière-plan via parfeval
    snapshotPeriod = 750;  % Fréquence de snapshot (ajustée en fonction du problème)
    f = parfeval(@simulationWithSnapshots, 0, T, Nt, snapshotPeriod, dq);
    
    while ~strcmp(f.State, 'finished')
        drawnow limitrate;   
        pause(0.05);
    end
    wait(f);  % S'assurer que la simulation est terminée
    
    stop(updateTimer);
    delete(updateTimer);
    
    runtime = toc;
    fprintf('Simulation complète en %.6f secondes\n', runtime);
    
    %% Nettoyage du pool à la fin
    delete(pool);
    
    %% --- Fonction de simulation (exécutée en parallèle) ---
    function simulationWithSnapshots(T0, Nt_local, snapshotPeriod_local, dq_local)
        T_sim = T0;
        % Tableaux locaux pour mesurer (pour le worker)
        localTherm1 = zeros(1, Nt_local);
        localTherm2 = zeros(1, Nt_local);
        localTherm3 = zeros(1, Nt_local);
        localEnergyAdded = zeros(1, Nt_local);
        localEnergyLoss = zeros(1, Nt_local);
        
        for t = 1:Nt_local
            % Calcul de Tnew par schéma numérique (exemple simplifié)
            Tnew = T_sim;
            % Conduction interne
            Tnew(2:Nx-1, 2:Ny-1) = T_sim(2:Nx-1, 2:Ny-1) ...
                + dt_dx2 * (T_sim(1:Nx-2, 2:Ny-1) - 2*T_sim(2:Nx-1, 2:Ny-1) + T_sim(3:Nx, 2:Ny-1)) ...
                + dt_dy2 * (T_sim(2:Nx-1, 1:Ny-2) - 2*T_sim(2:Nx-1, 2:Ny-1) + T_sim(2:Nx-1, 3:Ny));
            % Conduction aux bords (simplifié)
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
            
            % Mise à jour de T_sim et enregistrement des mesures locales
            T_sim = Tnew;
            localTherm1(t) = T_sim(Therm1_loc(1), Therm1_loc(2));
            localTherm2(t) = T_sim(Therm2_loc(1), Therm2_loc(2));
            localTherm3(t) = T_sim(Therm3_loc(1), Therm3_loc(2));
            localEnergyAdded(t) = sum(P(:)) * dt;
            localEnergyLoss(t) = (h_conv * aire_top * dt) * mean(T_sim(:) - T_piece);
            
            % Envoi des données toutes les snapshotPeriod_local itérations (et à t==1)
            if mod(t, snapshotPeriod_local) == 0 || t == 1
                dataOut.snapshot = T_sim;
                dataOut.t = t;
                dataOut.therm1 = localTherm1(t);
                dataOut.therm2 = localTherm2(t);
                dataOut.therm3 = localTherm3(t);
                dataOut.energyAdded = localEnergyAdded(t);
                dataOut.energyLoss = localEnergyLoss(t);
                send(dq_local, dataOut);
            end
        end
    end

    %% --- Callback du DataQueue ---
    function updateSnapshot(data)
        % Mise à jour du graphique 1 (surface) à chaque snapshot
        set(f1_surf, 'ZData', data.snapshot - 273.15);
        currentTime = data.t * dt;
        set(timeText, 'String', ['Temps : ' num2str(currentTime, '%.2f') ' s']);
        
        % Mise à jour des graphiques secondaires seulement si des données sont présentes
        persistent timeData therm1Data therm2Data therm3Data energyAddedData energyLossData updateCounter
        if isempty(timeData)
            timeData = currentTime;
            therm1Data = data.therm1;
            therm2Data = data.therm2;
            therm3Data = data.therm3;
            energyAddedData = data.energyAdded;
            energyLossData = data.energyLoss;
            updateCounter = 1;
        else
            timeData = [timeData, currentTime];
            if ~isempty(data.therm1)
                therm1Data = [therm1Data, data.therm1];
                therm2Data = [therm2Data, data.therm2];
                therm3Data = [therm3Data, data.therm3];
                energyAddedData = [energyAddedData, data.energyAdded];
                energyLossData = [energyLossData, data.energyLoss];
            end
            updateCounter = updateCounter + 1;
        end
        
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
