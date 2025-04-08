function indent_pertub_V2(csvfile, pert_power, start_time, end_time, list_of_options)
addpath("support")
%% Sortir les données 
res_exp = read_csv(csvfile);
y3  = res_exp.T3(1:end, 1);
y2  = res_exp.T2(1:end, 1);
y1  = res_exp.T1(1:end, 1);
t   = res_exp.Temps(1:end, 1);

if ~isnan(start_time)
    id_start = find(t >= start_time, 1);
    y3  = y3(1,id_start:end);
    y2  = y2(1,id_start:end);
    y1  = y1(1,id_start:end);
    t   = t(1,id_start:end);
end
if ~isnan(end_time)
    id_end = find(t >= end_time, 1);
    y3  = y3(1, 1:id_end-1);
    y2  = y2(1, 1:id_end-1);
    y1  = y1(1, 1:id_end-1);
    t   = t(1, 1:id_end-1);
end

%initialiser le vecteur de perturbation
p = ones(size(t))*pert_power;

% Retirer point opération
y1      = y1 - y1(1, 1); 
y2      = y2 - y2(1, 1);
y3      = y3 - y3(1, 1); 
t       = t - t(1,1);     
p(1,1)  = 0;

% List of options:
% 1: (bool) ident. de T1 
% 2: (bool) ident. de T2
% 3: (bool) ident. de T3 
% 4: (bool) ident. de perturb_T1
% 5: (bool) ident. de perturb_T2
% 6: (bool) ident. de perturb_T3
% 7: (bool) calcul du régulateur

if list_of_options(4)
    gpp_t1 = identify(y1, p, t, 1, 0, false, NaN);
        assignin('base', 'gpp_t1', gpp_t1);
end

if list_of_options(5)
    gpp_t2 = identify(y2, p, t, 1, 0, false, NaN);
        assignin('base', 'gpp_t2', gpp_t2);
end

if list_of_options(6)
    gpp_t3 = identify(y3, p, t, 1, 0, true, NaN);
        assignin('base', 'gpp_t3', gpp_t3);
end

end
