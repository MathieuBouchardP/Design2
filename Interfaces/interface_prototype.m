function interface_prototype()
    % Création de la fenêtre
    fig = figure('Name', 'Interface Prototype', 'Position', [1900, 1200, 1200, 1200]);

    % Variables globales
    global serialObj consigneData tempMesureeData timeData plotConsigne plotTemp plotSimu startTime tempSimuData;

    % Initialisation des données du graphe
    consigneData = [];
    tempMesureeData = [];
    timeData = [];
    tempSimuData = [];
    startTime = tic; % Chronomètre pour le temps écoulé

    % ========= GRAPHE 1 : Consigne & Température Mesurée =========
    ax1 = axes('Parent', fig, 'Position', [0.1, 0.6, 0.4, 0.35]);
    hold(ax1, 'on');
    plotConsigne = plot(ax1, timeData, consigneData, 'r', 'LineWidth', 2, 'DisplayName', 'Consigne');
    plotTemp = plot(ax1, timeData, tempMesureeData, 'b', 'LineWidth', 2, 'DisplayName', 'Température Mesurée');
    xlabel('Temps (s)');
    ylabel('Température (°C)');
    title('Consigne vs Température Mesurée');
    legend;
    grid on;
    hold(ax1, 'off');

    % ========= GRAPHE 2 : Simulation =========
    ax2 = axes('Parent', fig, 'Position', [0.1, 0.1, 0.4, 0.35]);
    hold(ax2, 'on');
    plotSimu = plot(ax2, [], [], 'g', 'LineWidth', 2, 'DisplayName', 'Température Simulation');
    xlabel('Temps (s)');
    ylabel('Température (°C)');
    title('Simulation');
    legend;
    grid on;
    hold(ax2, 'off');

    % ========= Interface utilisateur =========
    % Champ pour entrer la consigne
    uicontrol('Style', 'text', 'Position', [700, 300, 120, 30], 'String', 'Consigne (°C)');
    tempInput = uicontrol('Style', 'edit', 'Position', [850, 300, 120, 30]);

    % Bouton pour envoyer la consigne
    uicontrol('Style', 'pushbutton', 'Position', [1000, 300, 120, 30], ...
              'String', 'Envoyer', 'Callback', @(~, ~) sendConsigne(tempInput));

    % Bouton pour démarrer la connexion série
    uicontrol('Style', 'pushbutton', 'Position', [700, 400, 120, 40], ...
              'String', 'Démarrer COM', 'Callback', @connectSerial);

    % Bouton pour arrêter la connexion série
    uicontrol('Style', 'pushbutton', 'Position', [1000, 400, 120, 40], ...
              'String', 'Arrêter COM', 'Callback', @disconnectSerial);

    % Bouton pour charger la simulation
    uicontrol('Style', 'pushbutton', 'Position', [850, 200, 120, 30], ...
              'String', 'Charger Simulation', 'Callback', @loadSimulation);

    % Indicateur d'état de la connexion
    statusText = uicontrol('Style', 'text', 'Position', [840, 100, 150, 30], ...
                           'String', 'Statut : Déconnecté', 'ForegroundColor', 'red');

    % ========= FONCTIONS =========

    % Fonction pour démarrer la communication série
    function connectSerial(~, ~)
        try
            serialObj = serialport("COM3", 9600); % Adapter le port selon l'Arduino
            pause(2); % Pause pour stabiliser la connexion
            statusText.String = 'Statut : Connecté';
            statusText.ForegroundColor = 'green';
            % Démarrer la lecture des données de l'Arduino
            configureCallback(serialObj, "terminator", @readTemperature);
        catch
            statusText.String = 'Erreur de connexion';
            statusText.ForegroundColor = 'red';
        end
    end

    % Fonction pour arrêter la communication série
    function disconnectSerial(~, ~)
        if exist('serialObj', 'var') && ~isempty(serialObj)
            clear serialObj;
            statusText.String = 'Statut : Déconnecté';
            statusText.ForegroundColor = 'red';
        else
            statusText.String = 'Statut : Aucune connexion active';
            statusText.ForegroundColor = 'red';
        end
    end

    % Fonction pour lire la température envoyée par Arduino
    function readTemperature(src, ~)
        tempMesuree = str2double(readline(src)); % Lire et convertir la température
        if ~isnan(tempMesuree)
            % Ajouter les nouvelles données à la liste
            timeData = [timeData, toc(startTime)];
            tempMesureeData = [tempMesureeData, tempMesuree];
            
            % Simuler une température (par exemple, pour la démonstration)
            tempSimuData = [tempSimuData, tempMesuree + randn()*0.5]; % Température simulée avec bruit
            
            % Mettre à jour les graphiques
            updateGraphs();
        end
    end

    % Fonction pour mettre à jour les graphiques
    function updateGraphs()
        plotConsigne.YData = consigneData;
        plotTemp.YData = tempMesureeData;
        plotSimu.YData = tempSimuData;
        drawnow;
    end

    % Fonction pour envoyer la consigne
    function sendConsigne(tempInput)
        consigne = str2double(tempInput.String); % Récupérer la consigne
        if ~isnan(consigne)
            consigneData = [consigneData, consigne];
            % Mettre à jour les graphiques
            updateGraphs();
        end
    end

    % Fonction pour charger la simulation (optionnelle)
    function loadSimulation(~, ~)
        % Logique de chargement de simulation ici (exemple d'une fonction à compléter)
    end
end
