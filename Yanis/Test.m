close all
clear all
clc

F(length(zoom)) = struct('cdata',[],'colormap',[]);
writerObj = VideoWriter('TemperatureDistribution2D.avi');
open(writerObj);

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
    alpha = k/(rho*cp);    % Diffusivité Thermique [m^2/s]
    h_conv = params.conditions_limites.h_conv;  % Coeff. de convection[W/m^2·K]

    convection_activee = params.conditions_limites.convection_activee;    %activation ou non de la convection
    Pin = params.puissance.pin ;      % Puissance [W] (example: 5 W localized input)
    Pin_loc_x = params.puissance.pin_loc_x;     %Localisation de la puissance en x
    Pin_loc_y = params.puissance.pin_loc_y;     %Localisation de la puissance en y

    T_piece = params.conditions_initiales.T_piece ;  % Température pièce en Celsius
    T_loc_x = params.conditions_initiales.T_loc_x ;  % Localisation d'un dirac en Température sur x [m]
    T_loc_y = params.conditions_initiales.T_loc_y ;  % Localisation d'un dirac en Température sur y [m]

    Therm1_loc_x = params.conditions_initiales.Therm1_loc_x ;   % Localisation de la thermistance sur x [m]
    Therm1_loc_y = params.conditions_initiales.Therm1_loc_y ;   % Localisation de la thermistance sur y [m]
    Therm2_loc_x = params.conditions_initiales.Therm2_loc_x ;   % Localisation de la thermistance sur x [m]
    Therm2_loc_y = params.conditions_initiales.Therm2_loc_y ;   % Localisation de la thermistance sur y [m]
    Therm3_loc_x = params.conditions_initiales.Therm3_loc_x ;   % Localisation de la thermistance sur x [m]
    Therm3_loc_y = params.conditions_initiales.Therm3_loc_y ;   % Localisation de la thermistance sur y [m]
    
    %% Paramètres caclulés

    dx = Lx / Nx;               % Pas de discrétisation en x
    dy = Ly / Ny;               % Pas de discrétisation en y
    dz = epaisseur;             % épaisseur en z

    dt = (1/(4*alpha))*((dx^2*dy^2)/(dx^2+dy^2));    % Pas en temps [s], Choisi pour assurer la stabilité
    Nt = round(TempsTotal/dt);  % Nombre d'itérations temporelles

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
    Pin_Loc = [(fix(Pin_loc_x/dx)+ 1), (fix(Pin_loc_y/dy) + 1)] ;   %Élément qui reçoit la puissance de l'actuateur 
    
    P = zeros(Nx,Ny);           % Matrice de puissance à ajouter à chaque elt
    P(Pin_Loc(1), Pin_Loc(2)) = Pin;         % La puissance est mise sur un seul élément
                                            % En pratique, la puissance déposée sera probablement 
                                            % répartie sur plusieurs éléments. 
    
    %% Conditions initiales
    T_piece = 273.15+ T_piece;           % conversion Température pièce en [K]
    T = T_piece .* ones(Nx, Ny);  % Temperature de tous les éléments

    T_loc= [(fix(T_loc_x/dx) + 1), (fix(T_loc_y/dy)+ 1)] ; % Element qui est plus chaud
    
    T(T_loc(1), T_loc(2)) = 273.15 ;        % Un élement plus chaud
     
    Therm1_loc = [(fix(Therm1_loc_x/dx) + 1) ,(fix(Therm1_loc_y/dy) + 1)] ;   % Position de la mesure de température(Thermistance1)
    Therm2_loc = [(fix(Therm2_loc_x/dx) + 1) ,(fix(Therm2_loc_y/dy) + 1)] ;   % Position de la mesure de température(Thermistance2)
    Therm3_loc = [(fix(Therm3_loc_x/dx) + 1) ,(fix(Therm3_loc_y/dy) + 1)] ;   % Position de la mesure de température(Thermistance3)

   %% Préallocation des vecteurs qui seront utilisés dans la boucle

    energy_added = zeros(1,Nt);
    energy_loss = zeros(1,Nt);

    thermistance_1 = T_piece.*ones(1,Nt);
    thermistance_2 = T_piece.*ones(1,Nt);
    thermistance_3 = T_piece.*ones(1,Nt);
    Tnew = 0*T;

    %% COnfiguration de la figure

    f1 = figure(1);
    set(f1, 'Position', [0, 0, 2000, 1000]); % Position et taille optimisées
    set(f1, 'Color', 'w'); % Fond blanc pour un meilleur contraste
    rotate3d on; % Active l’interaction avec la souris

    %% Précalcul des constantes

    dt_dx2 = (alpha * dt) / dx^2;
    dt_dy2 = (alpha * dt) / dy^2;
    conv_term_top = (2 * aire_top * h_conv * dt) / (volume * rho * cp);
    conv_term_sides_y = (aire_sides_y * h_conv * dt) / (volume * rho * cp);
    conv_term_sides_x = (aire_sides_x * h_conv * dt) / (volume * rho * cp);
    power_term = dt / (rho * cp * volume);

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
    
    % Ajout de la puissance
    Tnew(Pin_Loc(1), Pin_Loc(2)) = Tnew(Pin_Loc(1), Pin_Loc(2)) + power_term * Pin;
    
    % Mise à jour de la température
    T = Tnew;
    thermistance1(t) = T(Therm1_loc(1), Therm1_loc(2));
    thermistance2(t) = T(Therm2_loc(1), Therm2_loc(2));
    thermistance3(t) = T(Therm3_loc(1), Therm3_loc(2));
    % Pour vérification, bilan d'énergie

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
    
  
    
    if mod(t, round(Nt/1000)) == 0 || t==1
       clf
    % Création de la figure
    subplot(131)
    surfc(1000*X, 1000*Y, T-273.15);
    set(gca, 'Color', 'w'); % Fond blanc
    shading interp;
    xlabel('Longueur de la plaque (mm)')
    ylabel('Largeur de la plaque (mm)')
    zlabel('Température en degré Celsius')
    title('Distribution de la température sur la plaque')
    colorbar; % Échelle de couleurs
    colormap jet; 
    grid on;
    axis tight; % Ajuste les axes aux données
    view(3); 
    axis vis3d;

    subplot(132); % Sélectionne le deuxième sous-graphique (1 ligne, 3 colonnes, position 2)

    % Tracer les courbes
    plot(Temps(1:t), thermistance1(1:t)-273.15 , 'g', 'DisplayName', 'Thermistance 1'); % Courbe verte
    hold on;
    plot(Temps(1:t), thermistance2(1:t) - 273.15, 'r', 'DisplayName', 'Thermistance 2'); % Courbe rouge
    plot(Temps(1:t), thermistance3(1:t) - 273.15, 'b', 'DisplayName', 'Thermistance 3'); % Courbe bleue
    hold off; 

    % Configuration du graphe
    grid on; % Active la grille

    % Personnalisation des axes
    ax = gca; % Récupère l'axe actuel
    ax.FontSize = 16; % Taille de la police pour les labels
    xlabel('Temps [s]', 'FontSize', 16);
    ylabel('Température [°C]', 'FontSize', 16);
    title('Température aux thermistances', 'FontSize', 16);

    % Ajouter une légende
    legend('show', 'FontSize', 14, 'Location', 'best'); % Affiche la légende
          
    subplot(133)
    hold on
    plot(Temps(1:t),energy_added(1:t))
    plot(Temps(1:t),energy_loss(1:t))
    xlabel('Temps [s]','FontSize',16)
    ylabel('Énergie dans l''itération','FontSize',16)
    legend('Energie déposée','Energie dissipée par convection','FontSize',16,'Location','southeast')
    grid on
    drawnow limitrate;
    end
end
close(writerObj);