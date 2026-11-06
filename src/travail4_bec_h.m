clc; close all; clear;

%% Code LDPC(n, lambda, rho)  -- seule zone à modifier pour changer de code
n = 128;
lam_e = 2;   lam_w = 1;      % lambda(X) = X^3        (exposants, poids)
rho_e = 5;   rho_w = 1;      % rho(X)    = X^5
% Exemple irrégulier : lam_e = [1 2]; lam_w = [0.5 0.5];  -> lambda(X) = 0.5 X + 0.5 X^2

nH        = 5;        % nombre de tirages de H (une courbe par tirage)
nitMax    = 100;      % itérations max du décodeur BP
nErrCible = 100;      % trames erronées à observer pour chaque epsilon
maxFrames = 2e5;      % nombre max de trames par epsilon
B         = 1000;     % trames décodées en parallèle

%% Conversion en vecteurs de coefficients : lam(i) = coef de X^(i-1)
lam = zeros(1, max(lam_e) + 1);  lam(lam_e + 1) = lam_w / sum(lam_w);
rho = zeros(1, max(rho_e) + 1);  rho(rho_e + 1) = rho_w / sum(rho_w);
dl = 1:numel(lam);  dr = 1:numel(rho);

R = 1 - sum(rho ./ dr) / sum(lam ./ dl);             % rendement de conception
fprintf('Rendement de conception R = %.4f\n', R);

% n doit donner un nombre entier de noeuds de chaque degré (variable et parité)
n0 = n;  n = n_valide(n0, lam, rho);
if n ~= n0, fprintf('n ajusté : %d -> %d (degrés entiers)\n', n0, n); end

%% Seuil BP par évolution de densité : x <- eps * lambda(1 - rho(1 - x))
lamf = @(x) polyval(fliplr(lam), x);
rhof = @(x) polyval(fliplr(rho), x);
lo = 0; hi = 1;
for k = 1:40
    e = (lo + hi) / 2;  x = e;
    for it = 1:20000
        xn = e * lamf(1 - rhof(1 - x));
        if abs(xn - x) < 1e-13, x = xn; break; end
        x = xn;
    end
    if x < 1e-6, lo = e; else, hi = e; end
end
eps_BP = lo;
fprintf('Seuil BP : %.4f\n', eps_BP);

%% Plage de epsilon centrée sur le seuil (pour (3,6) : environ [0.3, 0.5])
% pour imposer une plage : eps_vec = linspace(0.3, 0.5, 20);
eps_vec = linspace(max(0.02, eps_BP - 0.13), min(0.98, eps_BP + 0.07), 20);

%% Simulation
taux = nan(nH, numel(eps_vec));
for h = 1:nH
    Hd = double(tirage_H(n, lam, rho));              % m x n
    for ie = 1:numel(eps_vec)
        nErrFr = 0;  nBitsEff = 0;  nFr = 0;
        while nErrFr < nErrCible && nFr < maxFrames
            E = rand(B, n) < eps_vec(ie);            % effacements (mot nul émis)
            for it = 1:nitMax
                C   = double(E) * Hd.';              % effacés par parité (B x m)
                Rcv = (((C == 1) * Hd) > 0) & E;     % bits récupérés par une parité à 1 effacement
                if ~any(Rcv(:)), break; end          % plus d'effacement ou plus rien récupéré
                E = E & ~Rcv;
            end
            nErrFr   = nErrFr + sum(any(E, 2));
            nBitsEff = nBitsEff + sum(E(:));
            nFr      = nFr + B;
        end
        taux(h, ie) = nBitsEff / (nFr * n);
    end
    fprintf('tirage %d/%d terminé\n', h, nH);
end

%% Tracé
taux(taux == 0) = NaN;                               % rien observé : point non tracé en échelle log
figure; hold on; grid on;
semilogy(eps_vec, taux.', '-o', 'LineWidth', 1.3);
set(gca, 'YScale', 'log');
if all(isnan(taux(:))), ylim([1e-7 1]); end
yl = ylim;
plot([eps_BP eps_BP], yl, 'k--', 'LineWidth', 1.5);
xlim([min(eps_vec) max(max(eps_vec), eps_BP)]);
xlabel('\epsilon'); ylabel('Taux d''effacement binaire');
title(sprintf('LDPC, n = %d, R = %.3f', n, R));
legend([arrayfun(@(k) sprintf('H n°%d', k), 1:nH, 'UniformOutput', false), ...
        {sprintf('\\epsilon_{BP} = %.4f', eps_BP)}], 'Location', 'southeast');

%% Plus petit n >= n0 donnant des nombres entiers de noeuds de chaque degré
function n2 = n_valide(n0, lam, rho)
    dl = 1:numel(lam);  dr = 1:numel(rho);
    entier = @(v) all(abs(v - round(v)) < 1e-9);
    for k = n0:n0 + 2000
        E  = k / sum(lam ./ dl);                     % nombre de connecteurs
        nv = E * lam ./ dl;                          % noeuds de variable par degré
        if ~entier(nv), continue; end
        E  = sum(dl .* round(nv));
        if entier(E * rho ./ dr), n2 = k; return; end
    end
    warning('Aucun n valide trouvé près de %d : arrondi des degrés par tirage_H.', n0);
    n2 = n0;
end
