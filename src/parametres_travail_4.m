clc; close all; clear;

%% Code LDPC(n, lambda, rho)  -- seule zone à modifier pour changer de code
n = 128;
lam_e = 2;   lam_w = 1;      % lambda(X) = X^2        (exposants, poids)
rho_e = 5;   rho_w = 1;      % rho(X)    = X^5
% Exemple irrégulier : lam_e = [1 2]; lam_w = [0.5 0.5];  -> lambda(X) = 0.5 X + 0.5 X^2

Ntir = 10000;                % nombre de tirages de H

%% Conversion en vecteurs de coefficients : lam(i) = coef de X^(i-1)
lam = zeros(1, max(lam_e) + 1);  lam(lam_e + 1) = lam_w / sum(lam_w);
rho = zeros(1, max(rho_e) + 1);  rho(rho_e + 1) = rho_w / sum(rho_w);
dl = 1:numel(lam);  dr = 1:numel(rho);

% n doit donner un nombre entier de noeuds de chaque degré (variable et parité)
n0 = n;  n = n_valide(n0, lam, rho);
if n ~= n0, fprintf('n ajusté : %d -> %d (degrés entiers)\n', n0, n); end

%% Degrés, m, E et rendement de conception
nv = round(n / sum(lam ./ dl) * lam ./ dl);          % noeuds de variable par degré
E  = sum(dl .* nv);                                  % nombre de connecteurs
nc = round(E * rho ./ dr);                           % noeuds de parité par degré
m  = sum(nc);
R  = 1 - m / n;
for d = dl(nv > 0), fprintf('variable : %d noeuds de degré %d\n', nv(d), d); end
for d = dr(nc > 0), fprintf('parité   : %d noeuds de degré %d\n', nc(d), d); end
fprintf('n = %d, m = %d, E = %d, R = %.4f\n', n, m, E, R);

%% Statistiques sur Ntir tirages
multi  = zeros(Ntir, 1);     % au moins une arête multiple ?
nmulti = zeros(Ntir, 1);     % nombre de couples (i,j) avec arête multiple
Rreel  = zeros(Ntir, 1);     % rendement réel
for t = 1:Ntir
    [H, A] = tirage_H(n, lam, rho);
    multi(t)  = any(A(:) >= 2);
    nmulti(t) = sum(A(:) >= 2);
    Rreel(t)  = 1 - rang_gf2(H) / n;
end
fprintf('Proportion avec arête multiple : %.4f\n', mean(multi));
fprintf('Nb moyen d''arêtes multiples    : %.3f\n', mean(nmulti));
fprintf('Rendement réel moyen           : %.6f\n', mean(Rreel));

%% Fonctions locales
function r = rang_gf2(H)
% rang d'une matrice binaire sur F2 (élimination de Gauss)
r = 0;  [m, n] = size(H);  H = logical(H);
for j = 1:n
    p = find(H(r+1:m, j), 1);
    if isempty(p), continue; end
    p = p + r;
    H([r+1 p], :) = H([p r+1], :);
    lignes = find(H(:, j));  lignes(lignes == r+1) = [];
    H(lignes, :) = xor(H(lignes, :), H(r+1, :));
    r = r + 1;
    if r == m, break; end
end
end

function n2 = n_valide(n0, lam, rho)
% plus petit n >= n0 donnant des nombres entiers de noeuds de chaque degré
dl = 1:numel(lam);  dr = 1:numel(rho);
entier = @(v) all(abs(v - round(v)) < 1e-9);
for k = n0:n0 + 2000
    E  = k / sum(lam ./ dl);
    nv = E * lam ./ dl;
    if ~entier(nv), continue; end
    E  = sum(dl .* round(nv));
    if entier(E * rho ./ dr), n2 = k; return; end
end
warning('Aucun n valide trouvé près de %d : arrondi des degrés par tirage_H.', n0);
n2 = n0;
end