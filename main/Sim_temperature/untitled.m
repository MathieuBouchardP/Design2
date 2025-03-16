% Créer l'instance de l'application
%app = AppClass;
%app.state = 0;
%appConstant = parallel.pool.Constant(app);
%fct_simuler_plaque("param.json");
%Copy_of_sim_with_snapshots("param.json");
%sim_with_snapshots_old("param.json");
clear functions;
clear all;
slkdlajsdkad("param.json");
x = 3;
disp(x);
pause(20);


% Simulation de modifications de l'état depuis l'application
%pause(10);
%disp('Mise en pause de la simulation...');
%app.state = 1;

%pause(10);
%disp('Reprise de la simulation...');
%app.state = 0;

%pause(5);
%disp('Arrêt de la simulation...');
%app.state = 2;