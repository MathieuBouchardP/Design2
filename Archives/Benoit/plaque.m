classdef plaque
    properties (Access = public)
            TempsTotal ;  % Durée de la simulation [s]
            Lx    ;       % Longueur [m]
            Ly    ;       % Largeur [m]
            epaisseur;    % Épaisseur [m]
            Nx        ;   % Nombre d'éléments en x
            Ny ;          % Nombre d'éléments en y
            k  ;        % Conductivité thermique [W/m·K]
            rho  ;       % Densité [kg/m^3]
            cp  ;        % Chaleur spécifique [J/kg·K]
            h_conv   ;   % Coefficient de convection [W/m^2·K]
                Therm1_loc_x ;
                Therm1_loc_y ;
                Therm2_loc_x ;
                Therm2_loc_y ;
                Therm3_loc_x ;
                Therm3_loc_y ;
            pert_loc_x_min ;
            pert_loc_x_max ;
            pert_loc_y_min ;
            pert_loc_y_max ;
            pert_pow   ; % Puissance de la perturbation [W]
            t_pert_deb ; % Temps début perturbation [s]
            t_pert_fin ; % Temps fin perturbation [s]
            T_piece;
            T_init;
            pin ;
            pin_start_time ; % Le temps de début d'envoie de la puissance [s]
            pin_end_time ; % Le temps de fin d'envoie de la puissance [s]
            couplage_tec = 2.1; % Coefficient de couplage avec le TEC
            pin_loc_x_min ;
            pin_loc_x_max ;
            pin_loc_y_min;
            pin_loc_y_max ;
            state ;
            f1_surf;
            timeText;
            f2_t1;
            f2_t2;
            f2_t3;
    end
end