function Affichage_3D(data, T)
    % [input] data
    % colonne 1 : abscisses x
    % colonne 2 : ordonnées y
    % colonne 3 : températures

    % Extraction des valeurs x, y, et z
    x = data(:, 1) * 1000;
    y = data(:, 2)* 1000;
    z = T(:);

    % Déterminer la taille de la grille
    unique_x = unique(x); % Valeurs uniques de x
    unique_y = unique(y); % Valeurs uniques de y
    nx = length(unique_x); % Nombre de valeurs uniques en x
    ny = length(unique_y); % Nombre de valeurs uniques en y

    % Restructuration des données en matrices pour surf
    X = reshape(x, ny, nx);
    Y = reshape(y, ny, nx);
    Z = reshape(z, ny, nx);

    % Création de la figure
    f1 = figure(1);
    f1.Position = [1200 1000 1400 460];
    f1.Position =[600 1000 2000 460];
    surfc(X, Y, Z);
    xlabel('Longueur de la plaque (mm)')
    ylabel('Largeur de la plaque (mm)')
    zlabel('Température en degré Celsius')
    title('Distribution de la température sur la plaque')
    colorbar; % Échelle de couleurs
    colormap jet; 
    grid on;
end

