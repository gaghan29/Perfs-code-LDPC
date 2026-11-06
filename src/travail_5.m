clear; clc; close all;
addpath('src');
rng(0);                                 

n = 128;
lambda = [0 0 1];                         
rho    = [0 0 0 0 0 1];                   

%% Générer une matrice H avec l'algorithme PEG
H_peg = peg_ldpc(n, lambda, rho);
fprintf('H PEG : %d x %d, degrés des variables %s, degrés des parités %s\n', ...
        size(H_peg, 1), size(H_peg, 2), ...
        mat2str(unique(sum(H_peg, 1))), mat2str(unique(sum(H_peg, 2))'));
save('H_travail5.mat', 'H_peg');

%% Question 1 : maille du code PEG et des 5 matrices du Travail 4

load('H_travail4.mat', 'H_codes'); % les 5 matrices tirées au Travail 4
matrices = [{H_peg}, H_codes];
fprintf('\nMaille du graphe de Tanner :\n');

for k = 1:numel(matrices)
    H = matrices{k};
    cycle_4 = false;
    cycle_6 = false;
    for a = 1:n
        for b = a+1:n
            p_ab = find(H(:, a) & H(:, b));        % parités reliées à x_a et à x_b
            if numel(p_ab) >= 2
                cycle_4 = true;
            end
            if numel(p_ab) == 1
                for c = b+1:n
                    p_bc = find(H(:, b) & H(:, c));    % parité entre x_b et x_c
                    p_ca = find(H(:, c) & H(:, a));    % parité entre x_c et x_a
                    if numel(p_bc) == 1 && numel(p_ca) == 1 && ...
                       ~(p_ab == p_bc && p_bc == p_ca)
                        cycle_6 = true;
                    end
                end
            end
        end
    end
    if cycle_4
        maille = '4';
    elseif cycle_6
        maille = '6';
    else
        maille = '8 ou plus';
    end
    if k == 1
        fprintf('  H PEG : maille = %s\n', maille);
    else
        fprintf('  H n°%d (Algorithme 1) : maille = %s\n', k-1, maille);
    end
end

%% Question 2 : performances du code PEG, mêmes paramètres qu'au Travail 4

load('resultats_travail4.mat', 'liste_epsilon', 'taux', 'eps_BP');
nb_trames_err_min = 100;                 
nb_trames_max = 2e4;                      
nb_iter_max = 100;                    

taux_peg = zeros(1, numel(liste_epsilon));
fprintf('\nCode PEG\n');

for i = 1:numel(liste_epsilon)
    epsilon = liste_epsilon(i);
    nb_bits_effaces = 0;  nb_trames_err = 0;  nb_trames = 0;
    while nb_trames_err < nb_trames_err_min && nb_trames < nb_trames_max
        Lc = Inf(n, 1);                         % mot tout à zéro
        Lc(rand(n, 1) < epsilon) = 0;           % canal BEC
        dec = bp_bec_travail4(H_peg, Lc, nb_iter_max);
        nb_eff = sum(dec(end, :) == 0);         % bits encore effacés à la fin
        nb_bits_effaces = nb_bits_effaces + nb_eff;
        nb_trames_err = nb_trames_err + (nb_eff > 0);
        nb_trames = nb_trames + 1;
    end
    taux_peg(i) = nb_bits_effaces / (n * nb_trames);
    fprintf('  eps = %.3f : taux = %.2e (%d trames, %d erronées)\n', ...
            epsilon, taux_peg(i), nb_trames, nb_trames_err);
end
save('resultats_travail5.mat', 'liste_epsilon', 'taux_peg');

%PEG supprime les arêtes multiples et augmente la maille (6 au lieu de 4), mais
%près du seuil le gain est faible, car c'est la longueur finie qui limite ; il apparaît à faible ε.

% Figure du Travail 4, avec la courbe PEG en plus
figure;
semilogy(liste_epsilon, taux', '-o'); hold on; grid on;
semilogy(liste_epsilon, taux_peg, 'k-s', 'LineWidth', 2);
xline(eps_BP, 'k--', '\epsilon_{BP}');
xline(0.5, 'k:', '\epsilon_{Sh} = 1 - R');
xlim([0.3 0.51]);
xlabel('\epsilon (probabilité d''effacement du canal)');
ylabel('Taux d''effacement binaire après décodage');
legend([arrayfun(@(k) sprintf('H n°%d (Algorithme 1)', k), 1:size(taux, 1), ...
        'UniformOutput', false), {'H PEG (Algorithme 2)'}], 'Location', 'southeast');
title('LDPC(128, X^2, X^5) : tirage aléatoire et PEG');

%sum(H,1) / sum(H,2) poids des colonnes (degrés des x) / des lignes (degrés des c)
% nnz(H) nombre de 1 = nombre d'arêtes E
% nnz(H)/numel(H) densité de H
% [i, j] = find(H) liste des arêtes (ci, xj), parcourue colonne par colonne
% find(H(:,j)) parités reliées à xj
% H(:,a) & H(:,b) parités reliées à xa et à xb
% any(v) vrai si au moins un élément est non nul
% setdiff(A,B) / union(A,B) A privé de B / A ∪ B
% mod(A,2) réduction dans F2
% randperm(E) les entiers 1…E dans un ordre au hasard
% randi(k) un entier au hasard entre 1 et k
% rand(n,1) < epsilon chaque bit effacé avec la probabilité ε
% Inf(n,1) , zeros , ones vecteurs de +∞, de 0, de 1