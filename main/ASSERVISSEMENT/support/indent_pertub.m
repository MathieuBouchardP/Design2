function indent_pertub
addpath("support")

res_exp = read_csv("Essai7(Sheet1).csv");

cut = 1;
fin = 158;

y3 = res_exp.T3(cut:fin, 1); 
y2 = res_exp.T2(cut:fin, 1); 
y1 = res_exp.T1(cut:fin, 1); 
t = res_exp.Temps(cut:fin, 1);

y3 = y3(:) - y3(1, 1);
y2 = y2(:) - y2(1, 1);
y1 = y1(:) - y1(1, 1);

p_pert = 5^2/(25.12);
p_pert = ones(size(y1))*p_pert;
p_pert(1:9) = 0;

gpp_t1 = identify(y1, p_pert, t, 1, 0, false, NaN);
    assignin('base', 'gpp_t1', gpp_t1);

gpp_t2 = identify(y2, p_pert, t, 1, 0, false, NaN);
    assignin('base', 'gpp_t2', gpp_t2);

gpp_t3 = identify(y3, p_pert, t, 1, 0, true, 10);
    assignin('base', 'gpp_t3', gpp_t3);
end
