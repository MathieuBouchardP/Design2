function record_simulation_video(filename, fps)
    % Charger les données de simulation
    data = load(filename);
    T_record = data.T_record;            % Température en K (Nx x Ny x nframes)
    therm_record = data.therm_record;      % Températures des thermistances en K (3 x nframes)
    energy_record_added = data.energy_record_added;
    energy_record_loss = data.energy_record_loss;
    time_record = data.time_record;
    
    [Nx, Ny, nframes] = size(T_record);
    
    % Préparer les axes pour la distribution de température
    x = 1:Nx;
    y = 1:Ny;
    [Y, X] = meshgrid(y, x);
    
    % Création de la figure
    hFig = figure;
    
    % Sous-figure 1 : Distribution de température
    subplot(1,3,1)
    hMesh = meshc(X, Y, T_record(:,:,1) - 273.15); % Conversion en °C
    xlabel('X'); ylabel('Y'); zlabel('Temp [°C]');
    title(sprintf('Température, t = %.2f s', time_record(1)));
    
    % Sous-figure 2 : Courbes des thermistances
    subplot(1,3,2)
    hold on;
    hTherm1 = plot(time_record(1), therm_record(1,1)-273.15, 'r', 'DisplayName', 'Thermistance 1');
    hTherm2 = plot(time_record(1), therm_record(2,1)-273.15, 'g', 'DisplayName', 'Thermistance 2');
    hTherm3 = plot(time_record(1), therm_record(3,1)-273.15, 'b', 'DisplayName', 'Thermistance 3');
    hold off;
    xlabel('Temps [s]'); ylabel('Temp [°C]');
    title('Thermistances');
    legend('show');
    
    % Sous-figure 3 : Bilan énergétique
    subplot(1,3,3)
    hold on;
    hEAdded = plot(time_record(1), energy_record_added(1), 'k', 'DisplayName', 'Energie ajoutée');
    hELoss  = plot(time_record(1), energy_record_loss(1), 'm', 'DisplayName', 'Energie perdue');
    hold off;
    xlabel('Temps [s]'); ylabel('Energie');
    title('Bilan énergétique');
    legend('show');
    
    % Créer et configurer l'objet VideoWriter
    videoFilename = sprintf('simulation_video_%s.mp4', char(datetime('now', 'Format', 'yyyyMMdd_HHmmss')));
    v = VideoWriter(videoFilename, 'MPEG-4');
    v.FrameRate = fps;
    open(v);
    
    % Boucle d'animation
    for i = 1:nframes
        % Mise à jour de la distribution de température
        subplot(1,3,1)
        set(hMesh, 'ZData', T_record(:,:,i) - 273.15);
        title(sprintf('Température, t = %.2f s', time_record(i)));
        
        % Mise à jour des courbes des thermistances
        subplot(1,3,2)
        set(hTherm1, 'XData', time_record(1:i), 'YData', therm_record(1,1:i)-273.15);
        set(hTherm2, 'XData', time_record(1:i), 'YData', therm_record(2,1:i)-273.15);
        set(hTherm3, 'XData', time_record(1:i), 'YData', therm_record(3,1:i)-273.15);
        
        % Mise à jour des courbes d'énergie
        subplot(1,3,3)
        set(hEAdded, 'XData', time_record(1:i), 'YData', energy_record_added(1:i));
        set(hELoss, 'XData', time_record(1:i), 'YData', energy_record_loss(1:i));
        
        drawnow; % Force la mise à jour de la figure
        
        % Capturer la frame et l'écrire dans la vidéo
        frame = getframe(hFig);
        writeVideo(v, frame);
    end
    
    % Fermer l'objet VideoWriter
    close(v);
    fprintf('Vidéo sauvegardée sous : %s\n', videoFilename);
end
