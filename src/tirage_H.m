function [H, A] = tirage_H(n, lam, rho)
    % Tirage d'une matrice H dans l'ensemble LDPC(n, lambda, rho)
    % lam(i) : coefficient de X^(i-1) dans lambda(X), fraction d'arêtes reliées à un noeud de variable de degré i
    % rho(i) : idem pour les noeuds de parité
    % Exemple (3,6)-régulier : H = tirage_H(1000, [0 0 1], [0 0 0 0 0 1]);
    
    lam = lam(:).';  rho = rho(:).';
    dl = 1:numel(lam);  dr = 1:numel(rho);
    
    % 1: calcul de m, E, d_v et d_c
    E  = n / sum(lam ./ dl);                 % n = E * int(lambda)
    nv = round(E * lam ./ dl);               % nombre de noeuds de variable de chaque degré
    [~, i] = max(nv);  nv(i) = nv(i) + n - sum(nv);   % corrige l'arrondi pour avoir n noeuds
    dv = repelem(dl, nv);                    % degré de chaque noeud de variable
    E  = sum(dv);                            % nombre réel de connecteurs
    
    nc = round(E * rho ./ dr);               % nombre de noeuds de parité de chaque degré
    dc = repelem(dr, nc);                    % degré de chaque noeud de parité
    m  = numel(dc);
    if sum(dc) ~= E
        warning('Degrés non entiers : ajustement du dernier noeud de parité.');
        dc(end) = dc(end) + E - sum(dc);
    end
    
    % 2: noeud associé à chaque connecteur
    sv = repelem(1:n, dv);                   % s_v(e) : noeud de variable du connecteur e
    sc = repelem(1:m, dc);                   % s_c(e) : noeud de parité du connecteur e
    
    % 3: permutation uniforme
    p = randperm(E);
    
    % 4-8: A(i,j) = nombre d'arêtes entre c_i et x_j
    A = accumarray([sc(p).', sv.'], 1, [m n]);
    
    % 9
    H = mod(A, 2);
end