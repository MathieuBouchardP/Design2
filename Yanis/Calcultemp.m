
for t = 1:Nt            % Itération sur chaque temps
    Tnew = T;
    %Conduction extremité haut-gauche
    Tnew(1 , 1) = Tnew(1, 1) + ((k*dt)/(rho*cp*dx^2))*(T(1, 2) - 2*T(1 , 1) + T(2, 1));
    %Conduction extremité bas-gauche
    Tnew(Nx , 1) = Tnew(Nx, 1) + ((k*dt)/(rho*cp*dx^2))*(T(Nx, 2) - 2*T(Nx, 1) + T(Nx-1, 1));
    %Conduction extremité haut-droite
    Tnew(1 , Ny) = Tnew(1 , Ny) + ((k*dt)/(rho*cp*dx^2))*(T(1, Ny -1) - 2*T(1 , Ny) + T(2, Ny));
    %Conduction extremité bas-droite
    Tnew(Nx , Ny) = Tnew(Nx, Ny) + ((k*dt)/(rho*cp*dx^2))*(T(Nx-1, Ny) - 2*T(Nx, Ny) + T(Nx, Ny-1));
    %Conduction extremité  1ere ligne
    Tnew(1, 2:Ny-1) = Tnew(1, 2:Ny-1) + ((k*dt)/(rho*cp*dx^2)).*(T(1, 1:Ny-2 ) - 3.*T(1, 2:Ny-1) + T(1, 3:Ny) + T(2, 2:Ny-1) );
    %Conduction extremité dernière ligne
    Tnew(Nx, 2:Ny-1) = T(Nx, 2:Ny-1) + ((k*dt)/(rho*cp*dx^2)).*(T(Nx, 1:Ny-2 ) - 3.*T(Nx, 2:Ny-1) + T(1, 3:Ny) + T(Nx-1, 2:Ny-1) );
    %Conduction extremité 1ere colonne
    Tnew(2:Nx-1, 1) = T(2:Nx-1, 1) + ((k*dt)/(rho*cp*dx^2)).*(T(1:Nx-2, 1 ) - 3.*T(2:Nx-1, 1) + T(3:Nx, 1) + T(2:Nx-1, 2) );
    %Conduction extremité dernière colonne
    Tnew(2:Nx-1, Ny) = T(2:Nx-1, Ny) + ((k*dt)/(rho*cp*dx^2)).*(T(1:Nx-2, Ny ) - 3.*T(2:Nx-1, 1) + T(3:Nx, Ny) + T(2:Nx-1, Ny-1) );
    %Conduction des éléments au milieu
    Tnew(2:Nx-1, 2:Ny-1) = T(2:Nx-1, 2:Ny-1) + ((k*dt)/(rho*cp*dx^2)).*(T(1:Nx-2, 2:Ny-1) - 4.*T(2:Nx-1, Ny) + T(3:Nx, 2:Ny-1) + T(2:Nx-1, 1:Ny-2) + T(2:Nx-1, 3:Ny));
    %Convection en haut et en bas 
    Tnew(1:Nx, 1:Ny) = T(1:Nx, 1:Ny) + ((2*aire_top*h_conv*dt)/(volume*rho*cp)).* (T_piece - T(1:Nx, 1:Ny));
    %Convection première ligne
    Tnew(1, 1:Ny) = T(1, 1:Ny) + ((aire_sides_y*h_conv*dt)/(volume*rho*cp)).* (T_piece - T(1, 1:Ny));
    %Convection dernière ligne
    Tnew(Nx, 1:Ny) = T(Nx, 1:Ny) + ((aire_sides_y*h_conv*dt)/(volume*rho*cp)).* (T_piece - T(Nx, 1:Ny));
    %Convection première colonne
    Tnew(1:Nx,1) = T(1:Nx, 1) + ((aire_sides_x*h_conv*dt)/(volume*rho*cp)).* (T_piece - T(1:Nx, 1));
    %Convection dernière colonne
    Tnew(1:Nx, Ny) = T(1:Nx,Ny) + ((aire_sides_x*h_conv*dt)/(volume*rho*cp)).* (T_piece - T(1:Nx,Ny));
    %Rajout de la puissance
    Tnew(1:Nx, 1:Ny) = Tnew(1:Nx, 1:Ny) + (dt / (rho * cp *volume) ).* P(1:Nx, 1:Ny);
    % On a terminé la mise à jour de la température pour ce tour;
    T = Tnew;
    thermistance(t) = T(Therm_loc(1), Therm_loc(2));

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

