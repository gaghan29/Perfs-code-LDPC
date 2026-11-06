function H = peg_ldpc(n, lambda, rho)
% Construction PEG d'une matrice H (Algorithme 2 de l'énoncé).
%   lambda, rho : coefficients des polynômes, comme dans tirage_ldpc
%                 ex. lambda = [0 0 1] pour X^2, rho = [0 0 0 0 0 1] pour X^5
%   H : matrice de parité m x n

int_lambda = sum(lambda ./ (1:numel(lambda)));
int_rho    = sum(rho    ./ (1:numel(rho)));
m = round(n * int_rho / int_lambda);

dv = [];
for d = 1:numel(lambda)
    nb = round(n * (lambda(d)/d) / int_lambda);   
    dv = [dv, d * ones(1, nb)];
end

dc = [];
for d = 1:numel(rho)
    nb = round(m * (rho(d)/d) / int_rho);         
    dc = [dc, d * ones(1, nb)];
end
dc(end) = dc(end) + sum(dv) - sum(dc);           

H = zeros(m, n);

for j = 1:n
    for k = 1:dv(j)                               % ajout de la k-ième arête de x_j
        degre_courant = sum(H, 2)'; % degré courant de chaque parité
        S = find(degre_courant < dc); % parités non saturées

        if k == 1
            C = S; 

        else
            % Étage 0 : x_j et les parités déjà reliées à x_j
            Nl = find(H(:, j))'; % N^l : parités des étages 0 à l
            x_arbre = j; % variables déjà dans l'arbre

            while true
                % Étage l+1 : variables reliées aux parités de N, absentes de l'arbre,
                % puis parités reliées à ces variables
                x_nouvelles = setdiff(find(any(H(Nl, :), 1)), x_arbre);
                x_arbre = [x_arbre, x_nouvelles];
                Nl_suivant = union(Nl, find(any(H(:, x_nouvelles), 2))'); % N^{l+1}

                Nl_barre = setdiff(S, Nl); % N_barre^l = S \ N^l
                Nl_barre_suivant = setdiff(S, Nl_suivant);  % N_barre^{l+1} = S \ N^{l+1}

                cas_a = isequal(Nl_suivant, Nl) && ~isempty(Nl_barre);
                cas_b = ~isempty(Nl_barre) && isempty(Nl_barre_suivant);

                if cas_a || cas_b
                    break        
                end

                Nl = Nl_suivant; % l <- l + 1
            end

            C = Nl_barre;             
        end

        d_min = min(degre_courant(C));
        candidats = C(degre_courant(C) == d_min);
        i = candidats(randi(numel(candidats)));

        H(i, j) = 1;                            
    end
end
end