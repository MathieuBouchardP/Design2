function config(fichier)

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
    alpha = params.materiau.alpha;    % Diffusivité Thermique [m^2/s]
    h_conv = params.conditions_limites.h_conv;  % Coeff. de convection[W/m^2·K]

    convection_activee = params.conditions_limites.convection_activee;    %activation ou non de la convection
    Pin = params.puissance.pin ;      % Puissance [W] (example: 5 W localized input)
    Pin_loc_x = params.puissance.pin_loc_x;     %Localisation de la puissance en x
    Pin_loc_y = params.puissance.pin_loc_y;     %Localisation de la puissance en y

    T_piece = params.conditions_initiales.T_piece ;  % Température pièce en Celsius
    T_loc_x = params.conditions_initiales.T_loc_x ;  % Localisation d'un dirac en Température sur x [m]
    T_loc_y = params.conditions_initiales.T_loc_y ;  % Localisation d'un dirac en Température sur y [m]

    Therm_loc_x = params.conditions_initiales.Therm_loc_x ;   % Localisation de la thermistance sur x [m]
    Therm_loc_y = params.conditions_initiales.Therm_loc_y ;   % Localisation de la thermistance sur y [m]
    
    %% Paramètres caclulés

    dx = Lx / Nx;               % Pas de discrétisation en x
    dy = Ly / Ny;               % Pas de discrétisation en y
    dz = epaisseur;             % épaisseur en z

    dt = (1/(10*alpha))*((dx^2*dy^2)/(dx^2+dy^2));    % Pas en temps [s], Choisi pour assurer la stabilité
    Nt = round(TempsTotal/dt);  % Nombre d'itérations temporelles

    %% Pour la convection

    aire_sides_y = dy*dz ;      % aire des côtés sur y 
    aire_sides_x = dx*dz;       % aire des côtés sur y
    aire_top = dx*dy;           % aire du dessus et du dessous

    volume = dx*dy*dz;          % Volume d'un élément

    Temps = (0: Nt-1)*dt;       % Vecteur de temps pour la simulation
    x = ((1:Nx) - 0.5).*dx;     % Vecteur de coordonnées en x
    y = ((1:Ny) - 0.5 ).*dy;    % Vecteur de coordonnées en y
    
    [X, Y] = meshgrid(x, y);    % Génération des matrices coordonnées 2D
    Position = [X(:), Y(:)];     % Convertir X et Y en vecteurs colonne et les combiner

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

    T(T_loc(1), T_loc(2)) = 273.15 + 3;        % Un élement plus chaud
     
    Therm_loc = [(fix(Therm_loc_x/dx) + 1) ,(fix(Therm_loc_y/dy) + 1)] ;   % Position de la mesure de température(Thermistance)

   %% Préallocation des vecteurs qui seront utilisés dans la boucle

    energy_added = zeros(1,Nt);
    energy_loss = zeros(1,Nt);

    thermistance = zeros(1,Nt);
    Tnew = 0*T;
end 
