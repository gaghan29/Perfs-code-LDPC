clc; close all; clear;

n = 128;  lam = [0 0 1];  rho = [0 0 0 0 0 1];   % lambda = X^2, rho = X^5
Ntir = 10000;

% Degrés, m, E, rendement de conception
dv = 3;  dc = 6;
E = n*dv;  m = E/dc;  R = 1 - m/n;
fprintf('dv = %d, dc = %d, E = %d, m = %d, R = %.3f\n', dv, dc, E, m, R);

multi = zeros(Ntir, 1);
Rreel = zeros(Ntir, 1);
for t = 1:Ntir
    [H, A] = tirage_H(n, lam, rho);
    multi(t) = any(A(:) >= 2);              % au moins une arête multiple
    Rreel(t) = 1 - rang_gf2(H) / n;         % rendement réel
end
fprintf('Proportion avec arête multiple : %.4f\n', mean(multi));
fprintf('Rendement réel moyen           : %.4f\n', mean(Rreel));

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