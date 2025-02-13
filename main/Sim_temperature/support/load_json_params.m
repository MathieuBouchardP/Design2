function load_json_params(filename)
    % Lire le contenu du fichier JSON
    fid = fopen(filename, 'r');
    raw = fread(fid, inf, 'uint8=>char')';
    fclose(fid);
    
    % Décoder le JSON en structure MATLAB
    params = jsondecode(raw);
    
    % Assigner les valeurs aux variables dans l'espace de travail
    assignin('base', 'TempsTotal', params.simulation.TempsTotal);
    assignin('base', 'Lx', params.geometrie.Lx);
    assignin('base', 'Ly', params.geometrie.Ly);
    assignin('base', 'epaisseur', params.geometrie.epaisseur);
    
    assignin('base', 'Nx', params.grille.Nx);
    assignin('base', 'Ny', params.grille.Ny);
    
    assignin('base', 'k', params.materiau.k);
    assignin('base', 'rho', params.materiau.rho);
    assignin('base', 'cp', params.materiau.cp);
    assignin('base', 'materiau', params.materiau.nom);
    assignin('base', 'couplage_tec', params.puissance.couplage_TEC);
    
    
    alpha = params.materiau.k / (params.materiau.rho * params.materiau.cp);
    assignin('base', 'alpha', alpha);
    
    dx = params.geometrie.Lx / params.grille.Nx;
    dy = params.geometrie.Ly / params.grille.Ny;
    dz = params.geometrie.epaisseur;
    assignin('base', 'dx', dx);
    assignin('base', 'dy', dy);
    assignin('base', 'dz', dz);
    
    dt = (1/(4*alpha))*(dx^2 * dy^2)/(dx^2 + dy^2);
    Nt = round(params.simulation.TempsTotal / dt);
    assignin('base', 'dt', dt);
    assignin('base', 'Nt', Nt);
    
    assignin('base', 'h_conv', params.conditions_limites.h_conv);
    assignin('base', 'convection_activee', params.conditions_limites.convection_activee);
    assignin('base', 'Pin', params.puissance.pin);
    assignin('base', 'Pin_start_time', params.puissance.pin_start_time);

    Pin_end_time = params.puissance.pin_end_time;
    if isnan(Pin_end_time)
       Pin_end_time = params.simulation.TempsTotal + 1;
    end
    assignin('base', 'Pin_end_time', Pin_end_time);
    
    assignin('base', 'Pin_loc_x_min', fix(params.puissance.pin_loc_x_min/dx) + 1);
    assignin('base', 'Pin_loc_x_max', fix(params.puissance.pin_loc_x_max/dx) + 1);
    assignin('base', 'Pin_loc_y_min', fix(params.puissance.pin_loc_y_min/dy) + 1);
    assignin('base', 'Pin_loc_y_max', fix(params.puissance.pin_loc_y_max/dy) + 1);
    
    assignin('base', 'T_piece', params.conditions_initiales.T_piece);
    assignin('base', 'T_loc_x', fix(params.conditions_initiales.T_loc_x/dx) + 1);
    assignin('base', 'T_loc_y', fix(params.conditions_initiales.T_loc_y/dy) + 1);
    
    assignin('base', 'Therm1_loc_x', params.conditions_initiales.Therm1_loc_x);
    assignin('base', 'Therm1_loc_y', params.conditions_initiales.Therm1_loc_y);
    assignin('base', 'Therm2_loc_x', params.conditions_initiales.Therm2_loc_x);
    assignin('base', 'Therm2_loc_y', params.conditions_initiales.Therm2_loc_y);
    assignin('base', 'Therm3_loc_x', params.conditions_initiales.Therm3_loc_x);
    assignin('base', 'Therm3_loc_y', params.conditions_initiales.Therm3_loc_y);
    
    assignin('base', 'pert_loc_x_min', fix(params.pertub.pert_loc_x_min/dx) + 1);
    assignin('base', 'pert_loc_x_max', fix(params.pertub.pert_loc_x_max/dx) + 1);
    assignin('base', 'pert_loc_y_min', fix(params.pertub.pert_loc_y_min/dy) + 1);
    assignin('base', 'pert_loc_y_max', fix(params.pertub.pert_loc_y_max/dy) + 1);
    assignin('base', 'pert_pow', params.pertub.pert_pow);
    assignin('base', 't_pert_deb', params.pertub.t_pert_deb);
    assignin('base', 't_pert_fin', params.pertub.t_pert_fin);
end