function save_model(model, base_file_name,  folder_name)
%% Enregistre le modèle dans un sous dossier

% Rend l'argument base_file_name facultatif
if nargin < 3
    folder_name = "Identified_models"; % Valeur par défaut pour b
end

if nargin < 2
    base_file_name = "model"; % Valeur par défaut pour c
end


if ~exist(folder_name, 'dir')
    mkdir(folder_name);
end

counter = 1;
file_name = sprintf('%s_%d.mat', base_file_name, counter);
file_path = fullfile(folder_name, file_name);

while exist(file_path, 'file')
    counter = counter + 1;
    file_name = sprintf('%s_%d.mat', base_file_name, counter);
    file_path = fullfile(folder_name, file_name);
end
save(file_path, 'model');
fprintf('Fichier enregistré sous : %s\n', file_path);

end

