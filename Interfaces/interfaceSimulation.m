function interfaceSimulation()
    % Créer une fenêtre graphique principale avec une taille agrandie
    fig = figure('Name', 'Simulation Thermique', ...
                 'Position', [100, 100, 1000, 700], 'Color', [0.9 0.9 0.9]);

    % Titre
    uicontrol('Style', 'text', 'String', 'Simulation Thermique', ...
              'FontSize', 14, 'FontWeight', 'bold', ...
              'Position', [400, 650, 200, 30]);

    uicontrol('Style', 'text', 'String', 'Paramètres', ...
              'FontSize', 14, 'FontWeight', 'bold', ...
              'Position', [200, 470, 200, 40]);

    % === Section : Charger JSON ===
    uicontrol('Style', 'pushbutton', 'String', 'Charger JSON', ...
              'Position', [420, 100, 100, 40], 'FontSize', 10, ...
              'Callback', @chargerFichierJSON);

    % === Section : Paramètres ===
    uicontrol('Style', 'text', 'String', 'Longueur (m) :', ...
              'Position', [50, 400, 120, 35], 'FontSize', 10);
    handles.editLongueur = uicontrol('Style', 'edit', ...
                                     'Position', [180, 400, 90, 36]);

    uicontrol('Style', 'text', 'String', 'Largeur (m) :', ...
              'Position', [290, 400, 120, 35], 'FontSize', 10);
    handles.editLargeur = uicontrol('Style', 'edit', ...
                                    'Position', [420, 400,90, 36]);

    uicontrol('Style', 'text', 'String', 'Épaisseur (m) :', ...
              'Position', [50, 350, 120, 35], 'FontSize', 10);
    handles.editEpaisseur = uicontrol('Style', 'edit', ...
                                      'Position', [180, 350, 90, 36]);

    uicontrol('Style', 'text', 'String', 'Durée (s) :', ...
              'Position', [50, 250, 120, 35], 'FontSize', 10);
    handles.editDuree = uicontrol('Style', 'edit', ...
                                  'Position', [180, 250, 90, 36]);

    uicontrol('Style', 'text', 'String', 'Conductivité :', ...
              'Position', [50, 300, 120, 35], 'FontSize', 10);
    handles.editConductivite = uicontrol('Style', 'edit', ...
                                         'Position', [180, 300, 90, 36]);

    uicontrol('Style', 'text', 'String', 'Masse volumique :', ...
              'Position', [290, 350, 120, 35], 'FontSize', 10);
    handles.editmasse_volumique = uicontrol('Style', 'edit', ...
                                            'Position', [420, 350,90, 36]);

    uicontrol('Style', 'text', 'String', 'Capacité calorifique :', ...
              'Position', [290, 300, 120, 35], 'FontSize', 10);
    handles.editcapacite_calorifique = uicontrol('Style', 'edit', ...
                                                 'Position', [420, 300, 90, 36]);

    uicontrol('Style', 'text', 'String', 'Coefficient de convection :', ...
              'Position', [290, 250, 120, 35], 'FontSize', 10);
    handles.editcoefconvection = uicontrol('Style', 'edit', ...
                                           'Position', [420, 250, 90, 36]);

    % === Bouton : Démarrer la simulation ===
    uicontrol('Style', 'pushbutton', 'String', 'Simulation', ...
              'Position', [50, 100, 100, 40], 'FontSize', 10, 'FontWeight', 'bold', ...
              'Callback', @lancerSimulation);

    % === Zone pour afficher la courbe de température ===
    handles.axesGraph = axes('Parent', fig, 'Position', [0.6, 0.3, 0.3, 0.6]);
    title('Évolution de la température');
    xlabel('Position (m)');
    ylabel('Température (°C)');

    % Sauvegarder les éléments dans handles pour les récupérer plus tard
    guidata(fig, handles);
end

function chargerFichierJSON(~, ~)
    % Récupérer les handles de l'interface
    handles = guidata(gcf);

    % Demander à l'utilisateur de sélectionner un fichier JSON
    [fichier, chemin] = uigetfile('*.json', 'Sélectionnez un fichier JSON');

    % Vérifier si l'utilisateur a choisi un fichier
    if fichier ~= 0
        % Construire le chemin complet du fichier
        cheminFichier = fullfile(chemin, fichier);
        
        try
            % Lire et décoder le fichier JSON
            data = jsondecode(fileread(cheminFichier));

            % Vérifier si les champs nécessaires sont présents dans le JSON
            if isfield(data, 'plaque') && isfield(data.plaque, 'longueur') && ...
               isfield(data.plaque, 'largeur') && isfield(data.plaque, 'epaisseur') && ...
               isfield(data.plaque, 'conductivite') && isfield(data.plaque, 'masse_volumique') && ...
               isfield(data.plaque, 'capacite_calorifique') && isfield(data.simulation, 'duree') && ...
               isfield(data.simulation, 'coefficient_convection')

                % Mettre à jour les champs de l'interface avec les valeurs du JSON
                set(handles.editLongueur, 'String', num2str(data.plaque.longueur));
                set(handles.editLargeur, 'String', num2str(data.plaque.largeur));
                set(handles.editEpaisseur, 'String', num2str(data.plaque.epaisseur));
                set(handles.editConductivite, 'String', num2str(data.plaque.conductivite));
                set(handles.editmasse_volumique, 'String', num2str(data.plaque.masse_volumique));
                set(handles.editcapacite_calorifique, 'String', num2str(data.plaque.capacite_calorifique));
                set(handles.editDuree, 'String', num2str(data.simulation.duree));
                set(handles.editcoefconvection, 'String', num2str(data.simulation.coefficient_convection));

                disp('Les paramètres ont été chargés avec succès.');
            else
                disp('Erreur : Le fichier JSON est mal formaté.');
            end
        catch ME
            % En cas d'erreur dans la lecture du JSON
            disp(['Erreur lors du chargement du fichier JSON : ', ME.message]);
        end
    else
        disp('Aucun fichier sélectionné.');
    end
end

function lancerSimulation(~, ~)
    % Récupérer les handles de l'interface
    handles = guidata(gcf);

    % Récupérer les valeurs saisies dans les champs de l'interface
    longueur = str2double(get(handles.editLongueur, 'String'));
    largeur = str2double(get(handles.editLargeur, 'String'));
    epaisseur = str2double(get(handles.editEpaisseur, 'String'));
    duree = str2double(get(handles.editDuree, 'String'));
    conductivite = str2double(get(handles.editConductivite, 'String'));
    masse_volumique = str2double(get(handles.editmasse_volumique, 'String'));
    capacite_calorifique = str2double(get(handles.editcapacite_calorifique, 'String'));
    coefconvection = str2double(get(handles.editcoefconvection, 'String'));

    % Calculs pour la simulation (exemple, à adapter selon la physique du problème)
    % Ici, nous créons des valeurs fictives pour la température en fonction de la position
    x = linspace(0, longueur, 100);
    temperature = sin(pi * x / longueur) * 100; % Exemple de distribution de température

    % Afficher les résultats dans le graphique
    axes(handles.axesGraph);
    plot(x, temperature);
    title('Température en fonction de la position');
    xlabel('Position (m)');
    ylabel('Température (°C)');
end
