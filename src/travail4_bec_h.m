clc; close all; clear;

%% Paramètres
n   = 128;  lam = [0 0 1];  rho = [0 0 0 0 0 1];    % LDPC(128, X^2, X^5)
eps_vec = linspace(0.3, 0.5, 20);
nH      = 5;           % nombre de tirages de H (une courbe par tirage)
nitMax  = 100;         % itérations max du décodeur BP
nErrCible = 100;       % trames erronées à observer pour chaque epsilon
maxFrames = 2e5;       % nombre max de trames par epsilon
B = 1000;              % trames décodées en parallèle

%% Seuil BP par évolution de densité : x <- eps * lambda(1 - rho(1 - x))
lamf = @(x) polyval(fliplr(lam), x);
rhof = @(x) polyval(fliplr(rho), x);
lo = 0; hi = 1;
for k = 1:40
    e = (lo + hi) / 2;  x = e;
    for it = 1:5000
        x = e * lamf(1 - rhof(1 - x));
    end
    if x < 1e-6, lo = e; else, hi = e; end
end
eps_BP = lo;
fprintf('Seuil BP : %.4f\n', eps_BP);

%% Simulation
taux = nan(nH, numel(eps_vec));
for h = 1:nH
    Hd = double(tirage_H(n, lam, rho));          % m x n
    for ie = 1:numel(eps_vec)
        nErrFr = 0;  nBitsEff = 0;  nFr = 0;
        while nErrFr < nErrCible && nFr < maxFrames
            E = rand(B, n) < eps_vec(ie);        % effacements (mot nul émis)
            for it = 1:nitMax
                C   = double(E) * Hd.';          % effacés par parité (B x m)
                Rcv = (((C == 1) * Hd) > 0) & E; % bits récupérés par une parité à 1 effacement
                if ~any(Rcv(:)), break; end      % plus d'effacement ou plus rien récupéré
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
taux(taux == 0) = NaN;            % rien observé : point non tracé en échelle log
figure; hold on; grid on;
semilogy(eps_vec, taux.', '-o', 'LineWidth', 1.3);
set(gca, 'YScale', 'log');
yl = ylim;
plot([eps_BP eps_BP], yl, 'k--', 'LineWidth', 1.5);
xlabel('\epsilon'); ylabel('Taux d''effacement binaire');
legend([arrayfun(@(k) sprintf('H n°%d', k), 1:nH, 'UniformOutput', false), ...
        {sprintf('\\epsilon_{BP} = %.4f', eps_BP)}], 'Location', 'southeast');
