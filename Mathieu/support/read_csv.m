function dataStruct = read_csv(filename)
    if nargin < 1
        filename = "data.csv";
    end
    % Lire le fichier CSV
    data = readtable(filename);
    
    % Convertir chaque colonne en un vecteur et stocker dans une structure
    dataStruct = struct();
    columnNames = data.Properties.VariableNames;
    
    for i = 1:length(columnNames)
        dataStruct.(columnNames{i}) = data{:, i};
        %assignin('base', columnNames{i}, data{:, i}); % Assigner chaque colonne à une variable dans l'espace de travail
    end
end