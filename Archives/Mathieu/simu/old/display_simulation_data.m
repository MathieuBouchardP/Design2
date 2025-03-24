function display_simulation_data(filename, fps)
    % Lecture du fichier de simulation
    data = load(filename);
    T_record = data.T_record;            % Distribution de température (en K)
    therm_record = data.therm_record;      % Thermistances (en K)
    energy_record_added = data.energy_record_added;
    energy_record_loss = data.energy_record_loss;
    time_record = data.time_record;
    
    [Nx, Ny, nframes] = size(T_record);
    
    % Création d'un repère pour la distribution de température
    % Ici on utilise simplement les indices (on peut adapter si besoin)
    x = 1:Nx;
    y = 1:Ny;
    [Y, X] = meshgrid(y, x);
    
    % Création de la figure et des sous-figures
    figure;
    
    % Sous-figure 1 : Distribution de température
    h1 = subplot(1,3,1);
    hMesh = meshc(X, Y, T_record(:,:,1) - 273.15); % Conversion en °C
    xlabel('X'); ylabel('Y'); zlabel('Temp [°C]');
    title(sprintf('Température, t = %.2f s', time_record(1)));
    
    % Sous-figure 2 : Courbes des thermistances
    h2 = subplot(1,3,2);
    hold on;
    hTherm1 = plot(time_record(1), therm_record(1,1)-273.15, 'r', 'DisplayName', 'Thermistance 1');
    hTherm2 = plot(time_record(1), therm_record(2,1)-273.15, 'g', 'DisplayName', 'Thermistance 2');
    hTherm3 = plot(time_record(1), therm_record(3,1)-273.15, 'b', 'DisplayName', 'Thermistance 3');
    hold off;
    xlabel('Temps [s]'); ylabel('Temp [°C]');
    title('Thermistances');
    legend('show');
    
    % Sous-figure 3 : Bilan énergétique
    h3 = subplot(1,3,3);
    hold on;
    hEAdded = plot(time_record(1), energy_record_added(1), 'k', 'DisplayName', 'Energie ajoutée');
    hELoss  = plot(time_record(1), energy_record_loss(1), 'm', 'DisplayName', 'Energie perdue');
    hold off;
    xlabel('Temps [s]'); ylabel('Energie');
    title('Bilan énergétique');
    legend('show');
    
    % Calcul du temps de pause pour obtenir le fps désiré
    pauseTime = 1 / fps;
    
    % Boucle d'animation
    for i = 1:nframes
        % Mise à jour de la distribution de température
        set(hMesh, 'ZData', T_record(:,:,i));
        title(h1, sprintf('Température, t = %.2f s', time_record(i)));
        
        % Mise à jour des courbes de thermistances
        % set(hTherm1, 'XData', time_record(1:i), 'YData', therm_record(1,1:i)-273.15);
        % set(hTherm2, 'XData', time_record(1:i), 'YData', therm_record(2,1:i)-273.15);
        % set(hTherm3, 'XData', time_record(1:i), 'YData', therm_record(3,1:i)-273.15);
        
        % Mise à jour des courbes énergétiques
        % set(hEAdded, 'XData', time_record(1:i), 'YData', energy_record_added(1:i));
        % set(hELoss, 'XData', time_record(1:i), 'YData', energy_record_loss(1:i));
        
        drawnow;          % Forcer la mise à jour de l'affichage
        %pause(pauseTime); % Contrôle le nombre de frames par seconde
    end
end
