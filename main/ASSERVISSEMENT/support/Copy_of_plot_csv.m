%plotter
files = {
"Essai1_0.2V.csv",
"Essai1_0.6V.csv",
"Essai1_-0.2V.csv",
"Essai1_-0.8V.csv",
"Essai1_0.4V.csv",
"Essai1_-0.4V.csv",
"Essai1_-0.6V.csv"
};
names = {
    "E_0.2",
    "E_0.6",
    "E_-0.2",
    "E_-0.8",
    "E_0.4",
    "E_-0.4",
    "E_-0.6"
};
trucs = {
    0.2,
    0.6,
    -0.2,
    -0.8,
    0.4,
    -0.4,
    -0.6
    };


addpath("../Data")
addpath("../")

    % Sortir les data%
file = 2;
    res = read_csv(files{file});
    y3 = res.T3(:,1);
    y2 = res.T2(:,1);
    y1 = res.T1(:,1);  
    t = res.Temps(:,1);

    fig = figure('Position', [100, 100, 800, 450]);  
    hold on;
    % T3
        plot(t, y3, 'o');
    % T2
        plot(t, y2, 's');
    % T1

        plot(t, y1, '^');
    legend('Température T3 (prototype)',...
        'Température T2 (prototype)',...
        'Température T1 (prototype)', 'Location', 'best')
    xlabel('Temps [s]')
    ylabel('Température [°C]')
    ax.FontSize = 14;
    axis("auto y");
    xlim([0 810]);
    set(gca, 'LooseInset', get(gca, 'TightInset'))
    grid off;
    
    len = length(t);
    factor_stuff = 0.2;
    cuttof_index = round((1-factor_stuff)*len);
    
    f_t_1 = mean(y1(cuttof_index:end,1));
    disp(f_t_1);
    f_t_2 = mean(y2(cuttof_index:end,1));
    disp(f_t_2);
    f_t_3 = mean(y3(cuttof_index:end,1));
    disp(f_t_3);
    %print(gcf, sprintf('%s.png', names{file}), '-dpng', '-r600');  % Résolution de 600 dpi

