clc;
close all;
clear;
disp('debut tp')

H = [1 0 1 0 0 0
     0 1 0 1 1 0
     0 0 1 1 0 1];
s = size(H);

nit = 5;                       % nombre d'itérations
eps_vec = 0.1:0.05:1;          % probabilités d'effacement
N = 2000;                      % nombre de mots simulés par epsilon

taux = zeros(nit, numel(eps_vec));

for ie = 1:numel(eps_vec)
    for essai = 1:N

        % Canal BEC, mot de code nul : LLR = +Inf si bit reçu, 0 si effacé
        Lx = Inf(1, s(2));
        Lx(rand(1, s(2)) < eps_vec(ie)) = 0;

        % Initialisation des messages
        p2v = H;
        p2v(H == 1) = 0;
        p2v(H == 0) = NaN;
        v2p = p2v;

        for j = 1:nit
            % Mise à jour variable -> parité
            for k = 1:s(2)
                for h = 1:s(1)
                    if H(h,k) == 1
                        v2p(h,k) = Lx(k);
                        for i = 1:s(1)
                            if (H(i,k) == 1) && (i ~= h)
                                v2p(h,k) = v2p(h,k) + p2v(i,k);
                            end
                        end
                    end
                end
            end

            % Mise à jour parité -> variable
            for a = 1:s(1)
                for b = 1:s(2)
                    if H(a,b) == 1
                        produit = 1;
                        for i = 1:s(2)
                            if (H(a,i) == 1) && (i ~= b)
                                produit = produit * tanh(v2p(a,i)/2);
                            end
                        end
                        p2v(a,b) = 2*atanh(produit);
                    end
                end
            end

            % LLR a posteriori après cette itération
            Lpost = Lx;
            for k = 1:s(2)
                for i = 1:s(1)
                    if H(i,k) == 1
                        Lpost(k) = Lpost(k) + p2v(i,k);
                    end
                end
            end

            % bits encore effacés : LLR a posteriori nul
            taux(j, ie) = taux(j, ie) + sum(Lpost == 0);
        end
    end
end

taux = taux / (N * s(2));

figure;
semilogy(eps_vec, taux, '-o', 'LineWidth', 1.5);
grid on;
xlabel('\epsilon');
ylabel('Taux d''effacement binaire après décodage');
legend(arrayfun(@(k) sprintf('%d itération(s)', k), 1:nit, 'UniformOutput', false), 'Location', 'southeast');