close all
clear
clc
tic;
%F(length(zoom)) = struct('cdata',[],'colormap',[]);
%writerObj = VideoWriter('TemperatureDistribution2D.avi');
%open(writerObj);

%% Lire le contenu du fichier JSON
    fid = fopen('param.json', 'r');
    raw = fread(fid, inf, 'uint8=>char')';
    fclose(fid);
        
    %% Décoder le JSON en structure MATLAB
    params = jsondecode(raw);

    %% Accéder aux valeurs
    TempsTotal = params.simulation.TempsTotal;     % Durée de la simulation
    Lx = params.geometrie.Lx;       % Longueur [m]
    Ly = params.geometrie.Ly;       % Largeur [m]
    epaisseur = params.geometrie.epaisseur;     % epaisseur[m]

    Nx = params.grille.Nx;      % Nombre d'élements en x
    Ny = params.grille.Ny;      % Nombre d'élements en y

    k = params.materiau.k;      % Conductivité Thermique [W/m·K]
    rho = params.materiau.rho;  % Densité [kg/m^3]
    cp = params.materiau.cp;    % Chaleur spécifique  [J/kg·K]
    materiau = params.materiau.nom; % Nom du matériau
    alpha = k/(rho*cp);    % Diffusivité Thermique [m^2/s]
    %% Paramètres caclulés

    dx = Lx / Nx;               % Pas de discrétisation en x
    dy = Ly / Ny;               % Pas de discrétisation en y
    dz = epaisseur;             % épaisseur en z

    dt = (1/(4*alpha))*(dx^2 * dy^2)/(dx^2 + dy^2);    % Pas en temps [s], Choisi pour assurer la stabilité
    Nt = round(TempsTotal/dt);  % Nombre d'itérations temporelles
    %% Paramètres de Convection
    h_conv = params.conditions_limites.h_conv;  % Coeff. de convection[W/m^2·K]

    convection_activee = params.conditions_limites.convection_activee;    %activation ou non de la convection
    
    %% Paramètres de Puissance
    Pin = params.puissance.pin ;      % Puissance [W] (example: 5 W localized input)

    Pin_loc_x_min = params.puissance.pin_loc_x_min;     %Localisation de la puissance en x - min
    Pin_loc_x_min = (fix(Pin_loc_x_min/dx)+ 1);

    Pin_loc_x_max = params.puissance.pin_loc_x_max;     %Localisation de la puissance en x - min
    Pin_loc_x_max = (fix(Pin_loc_x_max/dx)+ 1);

    Pin_loc_y_min = params.puissance.pin_loc_y_min;     %Localisation de la puissance en y - min
    Pin_loc_y_min = (fix(Pin_loc_y_min/dy) + 1);

    Pin_loc_y_max = params.puissance.pin_loc_y_max;     %Localisation de la puissance en y - max
    Pin_loc_y_max = (fix(Pin_loc_y_max/dy) + 1);

    T_piece = params.conditions_initiales.T_piece ;  % Température pièce en Celsius

    T_loc_x = params.conditions_initiales.T_loc_x ;  % Localisation d'un dirac en Température sur x [m]
    T_loc_x = (fix(T_loc_x/dx) + 1);

    T_loc_y = params.conditions_initiales.T_loc_y ;  % Localisation d'un dirac en Température sur y [m]
    T_loc_y = (fix(T_loc_y/dy) + 1);

    %% Localisation des thermistances
    Therm1_loc_x = params.conditions_initiales.Therm1_loc_x ;   % Localisation de la thermistance sur x [m]
    Therm1_loc_y = params.conditions_initiales.Therm1_loc_y ;   % Localisation de la thermistance sur y [m]
    Therm2_loc_x = params.conditions_initiales.Therm2_loc_x ;   % Localisation de la thermistance sur x [m]
    Therm2_loc_y = params.conditions_initiales.Therm2_loc_y ;   % Localisation de la thermistance sur y [m]
    Therm3_loc_x = params.conditions_initiales.Therm3_loc_x ;   % Localisation de la thermistance sur x [m]
    Therm3_loc_y = params.conditions_initiales.Therm3_loc_y ;   % Localisation de la thermistance sur y [m]
    
    %% Paramètres de perturbations 
    pert_loc_x_min = params.pertub.pert_loc_x_min;
    pert_loc_x_min = fix(pert_loc_x_min/dx)+1;
    pert_loc_x_max = params.pertub.pert_loc_x_max;
    pert_loc_x_max = fix(pert_loc_x_max/dx)+1;
    pert_loc_y_min = params.pertub.pert_loc_y_min;
    pert_loc_y_min = fix(pert_loc_y_min/dy)+1;
    pert_loc_y_max = params.pertub.pert_loc_y_max;
    pert_loc_y_max = fix(pert_loc_y_max/dy)+1;
    pert_pow = params.pertub.pert_pow;
    t_pert_deb = params.pertub.t_pert_deb;
    t_pert_fin = params.pertub.t_pert_fin;

    %% Pour la convection

    aire_sides_y = dy*dz ;      % aire des côtés sur la largeur y
    aire_sides_x = dx*dz;       % aire des côtés sur y
    aire_top = dx*dy;           % aire du dessus et du dessous

    volume = dx*dy*dz;          % Volume d'un élément

    Temps = (0: Nt-1)*dt;       % Vecteur de temps pour la simulation
    x = ((1:Nx) - 0.5).*dx;     % Vecteur de coordonnées en x
    y = ((1:Ny) - 0.5 ).*dy;    % Vecteur de coordonnées en y
    
    [Y, X] = meshgrid(y, x);    % Génération des matrices coordonnées 2D
    
    %% Puissance déposée par l'actuateur
   
    P = zeros(Nx,Ny);           % Matrice de puissance à ajouter à chaque elt

    nb_elts_pin = (Pin_loc_y_max - Pin_loc_y_min + 1)* (Pin_loc_x_max - Pin_loc_x_min + 1); %nombre d'élements 
                                                                                            % couverts par l'actuateur
    P(Pin_loc_x_min:Pin_loc_x_max, Pin_loc_y_min:Pin_loc_y_max) = Pin/nb_elts_pin;    % La puissance est mise les élements indiqués
                                         
    
    %% Conditions initiales
    T_piece = 273.15 + T_piece;         % conversion Température pièce en [K]
    T = T_piece.*ones(Nx, Ny);          % Temperature de tous les éléments

    %T(T_loc_x, T_loc_y) =          % Un élement plus chaud
     
    Therm1_loc = [(fix(Therm1_loc_x/dx) + 1) ,(fix(Therm1_loc_y/dy) + 1)] ;   % Position de la mesure de température(Thermistance1)
    Therm2_loc = [(fix(Therm2_loc_x/dx) + 1) ,(fix(Therm2_loc_y/dy) + 1)] ;   % Position de la mesure de température(Thermistance2)
    Therm3_loc = [(fix(Therm3_loc_x/dx) + 1) ,(fix(Therm3_loc_y/dy) + 1)] ;   % Position de la mesure de température(Thermistance3)

   %% Préallocation des vecteurs qui seront utilisés dans la boucle

    energy_added = zeros(1,Nt);
    energy_loss = zeros(1,Nt);

    thermistance_1 = T_piece.*ones(1,Nt);
    thermistance1 = zeros(1, Nt);
    thermistance_2 = T_piece.*ones(1,Nt);
    thermistance2 = zeros(1, Nt);
    thermistance_3 = T_piece.*ones(1,Nt);
    thermistance3 = zeros(1, Nt);
    Tnew = T;

    %% Configuration de la figure

    f1 = figure(1);
    sgtitle(strcat("Distribution de température sur une plaque d'", materiau));
    set(gcf, 'Units', 'normalized', 'OuterPosition', [0 0 1 1]); % Position et taille optimisées
    set(f1, 'Color', 'w'); % Fond blanc pour un meilleur contraste
    rotate3d on; % Active l’interaction avec la souris
    
    %% Précalcul des constantes

    dt_dx2 = (alpha * dt) / dx^2;
    dt_dy2 = (alpha * dt) / dy^2;
    conv_term_top = (2 * aire_top * h_conv * dt) / (volume * rho * cp);
    conv_term_sides_y = (aire_sides_y * h_conv * dt) / (volume * rho * cp);
    conv_term_sides_x = (aire_sides_x * h_conv * dt) / (volume * rho * cp);
    power_term = dt / (rho * cp * volume);
    deposited = 0;

%% Création de la figure
    % Initier le graphique 1
    subplot(131)

    f1_surf = meshc(1000*X, 1000*Y, T-273.15);
    shading interp;
    xlabel('Longueur de la plaque (mm)')
    ylabel('Largeur de la plaque (mm)')
    zlabel('Température en degré Celsius')
    %zlim([20 40]);
    title(['Température à t = ', num2str(69), ' s'],'FontSize',16);
    timeText = title(['Temps : ' num2str(0, '%.2f') ' s'], 'FontSize', 16);
    colorbar; % Échelle de couleurs
    colormap jet; 
    clim([10 50]);
    grid on;
    view(3); 
    zticks(floor(0):1:ceil(100)); % Forcer à ce qu'il gradu à tout les degrée
    axis manual;
    pbaspect([Lx Ly min([Lx Ly])]); % forcer le ratio
    axis 'auto z' % Libérer l'axe des z

    subplot(132); % Sélectionne le deuxième sous-graphique (1 ligne, 3 colonnes, position 2)
    % Initier le graphique 2
    hold on;
    f2_t1 = plot(Temps(1:2), thermistance1(1:2)-273.15 , 'r', 'DisplayName', 'Thermistance 1'); % Courbe verte
    f2_t2 = plot(Temps(1:2), thermistance2(1:2) - 273.15, 'g', 'DisplayName', 'Thermistance 2'); % Courbe rouge
    f2_t3 = plot(Temps(1:2), thermistance3(1:2) - 273.15, 'b', 'DisplayName', 'Thermistance 3'); % Courbe bleue
    grid on; 
    ax = gca; % Récupère l'axe actuel
    ax.FontSize = 16; % Taille de la police pour les labels
    xlabel('Temps [s]', 'FontSize', 16);
    ylabel('Température [°C]', 'FontSize', 16);
    title('Température aux thermistances', 'FontSize', 16);
    legend('show', 'FontSize', 14, 'Location', 'best'); % Affiche la légende
    
    % initier le graphique 3
    subplot(133);
    hold on;
    f3_add = plot(Temps(1:2),energy_added(1:2));
    f3_loss = plot(Temps(1:2),energy_loss(1:2));
    xlabel('Temps [s]','FontSize',16)
    ylabel('Énergie dans l''itération','FontSize',16)
    legend('Energie déposée','Energie dissipée par convection','FontSize',16,'Location','southeast')
    grid on
%% Initiation du multi-thread
    %pool = gcp('nocreate');
    %if isempty(pool)
    %    pool = parpool; % Crée un pool si aucun n'est actif
    %end


%% Update de l'affichage (obselète) ne pas supprimer

function update_display(f1_surf, surf, timeText, t, f2_t1, t1, f2_t2, t2, f2_t3, t3, f3_add, add, f3_loss, loss, Temps, dt)
        set(f1_surf , 'ZData', surf - 273.15);
        set(timeText, 'String', ['Temps : ' num2str(t * dt, '%.2f') ' s']);
    
        % Mise à jour des courbes de température
        set(f2_t1, 'XData', Temps(1:t), 'YData', t1(1:t) - 273.15);
        set(f2_t2, 'XData', Temps(1:t), 'YData', t2(1:t) - 273.15);
        set(f2_t3, 'XData', Temps(1:t), 'YData', t3(1:t) - 273.15);
    
        % Mise à jour des courbes d’énergie
        set(f3_add, 'XData', Temps(1:t), 'YData', add(1:t));
        set(f3_loss, 'XData', Temps(1:t), 'YData', loss(1:t));
        axis auto;
        drawnow limitrate;
end
% Boucle principale
for t = 1:Nt
    Tnew = T;
    
    % Conduction interne (milieu)
    Tnew(2:Nx-1, 2:Ny-1) = T(2:Nx-1, 2:Ny-1) ...
        + dt_dx2 * (T(1:Nx-2, 2:Ny-1) - 2*T(2:Nx-1, 2:Ny-1) + T(3:Nx, 2:Ny-1)) ...
        + dt_dy2 * (T(2:Nx-1, 1:Ny-2) - 2*T(2:Nx-1, 2:Ny-1) + T(2:Nx-1, 3:Ny));
    
    % Conduction aux bords
    Tnew(1, 2:Ny-1) = T(1, 2:Ny-1) ... % Première ligne
        + dt_dx2 * (T(2, 2:Ny-1) - T(1, 2:Ny-1)) ...
        + dt_dy2 * (T(1, 1:Ny-2) - 2.*T(1, 2:Ny-1) + T(1, 3:Ny));
    
    Tnew(Nx, 2:Ny-1) = T(Nx, 2:Ny-1) ... % Dernière ligne
        + dt_dx2 * (T(Nx-1, 2:Ny-1) - T(Nx, 2:Ny-1)) ...
        + dt_dy2 * (T(Nx, 1:Ny-2) - 2.*T(Nx, 2:Ny-1) + T(Nx, 3:Ny));
    
    Tnew(2:Nx-1, 1) = T(2:Nx-1, 1) ... % Première colonne
        + dt_dy2 * (T(1:Nx-2, 1) - 2.*T(2:Nx-1, 1) + T(3:Nx, 1)) ...
        + dt_dx2 * (T(2:Nx-1, 2) - T(2:Nx-1, 1));
    
    Tnew(2:Nx-1, Ny) = T(2:Nx-1, Ny) ... % Dernière colonne
        + dt_dy2 * (T(1:Nx-2, Ny) - 2.*T(2:Nx-1, Ny) + T(3:Nx, Ny)) ...
        + dt_dx2 * (T(2:Nx-1, Ny-1) - T(2:Nx-1, Ny));
    
    % Conduction aux coins
    Tnew(1, 1) = T(1, 1) ... % Coin haut-gauche
        + dt_dx2 * (T(2, 1) - T(1, 1) ) ...
        + dt_dy2 * (T(1, 2) - T(1, 1));
    
    Tnew(1, Ny) = T(1, Ny) ... % Coin haut-droite
        + dt_dx2 * (T(2, Ny) - T(1, Ny)) ...
        + dt_dy2 * (T(1, Ny-1) - T(1, Ny));
    
    Tnew(Nx, 1) = T(Nx, 1) ... % Coin bas-gauche
        + dt_dx2 * (T(Nx-1, 1) - T(Nx, 1)) ...
        + dt_dy2 * (T(Nx, 2) - T(Nx, 1));
    
    Tnew(Nx, Ny) = T(Nx, Ny) ... % Coin bas-droite
        + dt_dx2 * (T(Nx-1, Ny) - T(Nx, Ny) ) ...
        + dt_dy2 * (T(Nx, Ny-1) - T(Nx, Ny) );
    
    % Convection
    Tnew(1:Nx,1:Ny) = Tnew(1:Nx,1:Ny) - conv_term_top .* (T(1:Nx,1:Ny) - T_piece);% Convection en haut et en bas
    Tnew(1, :) = Tnew(1, :) - conv_term_sides_y .* (T(1, :) - T_piece); % Première ligne
    Tnew(Nx, :) = Tnew(Nx, :) - conv_term_sides_y .* (T(Nx, :) - T_piece); % Dernière ligne
    Tnew(:, 1) = Tnew(:, 1) - conv_term_sides_x .* (T(:, 1) - T_piece); % Première colonne
    Tnew(:, Ny) = Tnew(:, Ny) - conv_term_sides_x .* (T(:, Ny) - T_piece); % Dernière colonne
    
    % Ajout des perturbations
     
    if pert_pow ~= 0  && deposited == 0 
        if (round(t_pert_deb/dt) == t)
            nb_elts_pert = (pert_loc_y_max - pert_loc_y_min + 1)* (pert_loc_x_max - pert_loc_x_min + 1);
            ajout = (pert_pow/nb_elts_pert);
            P(pert_loc_x_min:pert_loc_x_max, pert_loc_y_min:pert_loc_y_max) = ...
                P(pert_loc_x_min:pert_loc_x_max, pert_loc_y_min:pert_loc_y_max) + ajout;
            deposited = deposited + 1;
        end
    end
  
    if pert_pow ~= 0  && deposited == 1
        if (t == round(t_pert_fin/dt))
            P(pert_loc_x_min:pert_loc_x_max, pert_loc_y_min:pert_loc_y_max) = ...
                P(pert_loc_x_min:pert_loc_x_max, pert_loc_y_min:pert_loc_y_max) - ajout;
            deposited = deposited + 1;
        end 
    end
    % Ajout de la puissance
    Tnew = Tnew + power_term .* P;

 
    
    % Mise à jour de la température
    T = Tnew;
    thermistance1(t) = T(Therm1_loc(1), Therm1_loc(2));
    thermistance2(t) = T(Therm2_loc(1), Therm2_loc(2));
    thermistance3(t) = T(Therm3_loc(1), Therm3_loc(2));

   %% Pour vérification, bilan d'énergie

    % Énergie ajoutée dans le systeme par l'actuateur
    energy_added(t) = sum(P(:)) * dt;

    % Énergie dissipée par convection
    energy_loss_sides_ligne_1  = (h_conv * aire_sides_y * dt)* sum(T(1, 1:Ny) - T_piece);
    energy_loss_sides_ligne_f  =(h_conv * aire_sides_y * dt)*sum(T(Nx, 1:Ny) - T_piece) ;

    energy_loss_top_down    = (h_conv * 2*aire_top * dt)*sum(sum(T(1:Nx, 1:Ny) - T_piece));
    
    energy_loss_sides_colonne_1  = (h_conv * aire_sides_x * dt)* sum(T(1:Nx, 1) - T_piece);
    energy_loss_sides_colonne_f = (h_conv * aire_sides_x * dt)* sum(T(1:Nx, Ny) - T_piece);

    % Énergie totale dissipée
    energy_loss(t) = energy_loss_sides_ligne_1 + energy_loss_sides_ligne_f + energy_loss_top_down + energy_loss_sides_colonne_1 + energy_loss_sides_colonne_f;
    
    if abs((energy_added(t) - energy_loss(t))/(energy_added(t)+1e-9)) < 1e-3 
        runtime = toc;
        fprintf('Temps d''exécution : %.6f secondes\n', runtime);
        tic;
        update_display(f1_surf, Tnew, timeText, t, f2_t1, thermistance1, f2_t2, thermistance2, f2_t3, thermistance3, f3_add, energy_added, f3_loss, energy_loss, Temps, dt);
        fig_time = toc;
        fprintf('Temps d''affichage : %.6f secondes\n', fig_time);
        break
    end
    
    %% Affichage
    
    if mod(t, round(Nt/10)) == 0 || t==1 % affichage en 1000 intervale
    %if mod(t, 100) == 0 || t==1  % affichage à chaque 5000 iteration
    %if t == Nt   % Mode où on affiche juste le résultat final

        runtime = toc;
        %parfeval(pool, @update_display, 0, f1_surf, Tnew, timeText, t, f2_t1, thermistance1, f2_t2, thermistance2, f2_t3, thermistance3, f3_add, energy_added, f3_loss, energy_loss, Temps, dt);
        fprintf('Temps d''exécution : %.6f secondes\n', runtime);
        update_display(f1_surf, Tnew, timeText, t, f2_t1, thermistance1, f2_t2, thermistance2, f2_t3, thermistance3, f3_add, energy_added, f3_loss, energy_loss, Temps, dt);
        %tic

        % Mise à jour de la plaque
        %set(f1_surf , 'ZData', Tnew - 273.15);
        %set(timeText, 'String', ['Temps : ' num2str(t * dt, '%.2f') ' s']);
    
        % Mise à jour des courbes de température
        %set(f2_t1, 'XData', Temps(1:t), 'YData', thermistance1(1:t) - 273.15);
        %set(f2_t2, 'XData', Temps(1:t), 'YData', thermistance2(1:t) - 273.15);
        %set(f2_t3, 'XData', Temps(1:t), 'YData', thermistance3(1:t) - 273.15);
    
        % Mise à jour des courbes d’énergie
        %set(f3_add, 'XData', Temps(1:t), 'YData', energy_added(1:t));
        %set(f3_loss, 'XData', Temps(1:t), 'YData', energy_loss(1:t));
        %axis auto;
        %drawnow limitrate;
        %F = getframe(gcf);
        %writeVideo(writerObj,F);
    end
end
%fig_time = toc;
%fprintf('Temps d''affichage : %.6f secondes\n', fig_time);

%close(writerObj);