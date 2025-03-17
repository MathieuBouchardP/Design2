function info = Reg_generation(Kc, Ti, Td, Tf)
   % Définir la période d'échantillonnage
    Ts = 5;
    % Fonction de transferts du régulateur
     DF = tf([Td, 1],[Tf, 1]);

     PI = tf(1, [Ti, 1]);
    % Conversion en transformée en Z (Tustin)
    DFz = c2d(DF, Ts, 'tustin');
    PIz = c2d(PI, Ts, 'tustin');
    % Extraction des coefficients (SISO uniquement)
    DFrec = [cell2mat(DFz.Numerator); -cell2mat(DFz.Denominator)];
    PIrec = [cell2mat(PIz.Numerator); -cell2mat(PIz.Denominator)];

   % Mettre 0 dans la première valeur du dénominateur de DFrec et PIrec
    DFrec(2,1) = 0;
    PIrec(2,1) = 0;

    info = {Kc, DFrec, PIrec};
end

