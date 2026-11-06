
clear;
clc;
close all;

%% Parametres numeriques
tol = 1e-10;       % Seuil pour considerer x comme nul
maxIter = 5000;    % Nombre maximal d'iterations
precision = 1e-4;  % Precision demandee pour le seuil

%% Dichotomie
epsilon_min = 0;
epsilon_max = 1;

while (epsilon_max - epsilon_min) > precision

    epsilon_test = (epsilon_min + epsilon_max)/2;

    % Test de convergence de la suite x_l
    succes = evolution_densite(epsilon_test, tol, maxIter);

    if succes
        % La suite converge vers zero :
        % le seuil est au-dessus de epsilon_test
        epsilon_min = epsilon_test;
    else
        % La suite ne converge pas vers zero :
        % le seuil est en dessous de epsilon_test
        epsilon_max = epsilon_test;
    end

end

%% Estimation finale du seuil
epsilon_BP = (epsilon_min + epsilon_max)/2;

fprintf('Seuil BP estime : %.6f\n', epsilon_BP);
fprintf('Intervalle final : [%.6f ; %.6f]\n', ...
        epsilon_min, epsilon_max);


%% Fonction d'evolution de densite
function succes = evolution_densite(epsilon, tol, maxIter)

    % Initialisation
    x = epsilon;

    succes = false;

    % Iterations de l'evolution de densite
    for ell = 1:maxIter

        % x_l = epsilon * lambda(1-rho(1-x_(l-1)))
        % Avec lambda(X)=X^2 et rho(X)=X^5
        x_nouveau = epsilon * (1-(1-x)^5)^2;

        % Critere de convergence vers zero
        if x_nouveau < tol
            succes = true;
            return;
        end

        % Detection d'une valeur fixe non nulle
        if abs(x_nouveau-x) < 1e-14
            succes = false;
            return;
        end

        % Passage a l'iteration suivante
        x = x_nouveau;

    end

    % Si maxIter est atteint sans atteindre tol,
    % on considere que la convergence vers zero
    % n'a pas ete demontree numeriquement.
    succes = false;

end



% Traces des fonctions f(epsilon,x)

eps1 = 0.35;
eps2 = epsilon_BP;
eps3 = 0.80;

figure;
hold on;
grid on;

% Droite y = x
plot([0 1], [0 1], 'k--', 'LineWidth', 1.5);

% Domaine de x
x = linspace(0, 1, 1000);

% Courbe epsilon = 0.35
x1 = x(x <= eps1);
plot(x1, eps1*(1-(1-x1).^5).^2, ...
    'b', 'LineWidth', 1.8);

% Courbe epsilon = epsilon_BP
x2 = x(x <= eps2);
plot(x2, eps2*(1-(1-x2).^5).^2, ...
    'r', 'LineWidth', 1.8);

% Courbe epsilon = 0.47
x3 = x(x <= eps3);
plot(x3, eps3*(1-(1-x3).^5).^2, ...
    'g', 'LineWidth', 1.8);

xlabel('x');
ylabel('f(\epsilon,x)');
title('Evolution de densite et seuil BP');

legend('y = x', ...
       '\epsilon = 0.35', ...
       '\epsilon = \epsilon_{BP}', ...
       '\epsilon = 0.47', ...
       'Location', 'northwest');

hold off;
