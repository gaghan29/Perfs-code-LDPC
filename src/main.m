clc;
close all;
clear;

disp('debut tp')

nit = 50;

H = [[1 0 1 0 0 0]
    [0 1 0 1 1 0]
    [0 0 1 1 0 1]];

s = size(H);

p2v = H;

for p = 1:s(1)
    for l = 1:s(2)
        if H(p,l) == 1
            p2v(p,l) = 0;
        else
            p2v(p,l) = NaN;
        end
        
    end    
end

v2p = p2v;

Lx = [-1,1,2,-1,1,-1];

disp(v2p)
disp(p2v)
disp(Lx)

for j = 1:nit
    % Mise à jour de la matrice de noeuds de variable vers noeuds de parité
    for k = 1:s(2)
        for h = 1:s(1)
            if ~isnan(v2p(h,k))
                v2p(h,k) = Lx(k);
                for i = 1:s(1)
                    if (H(i,k) == 1) && (i ~= h) 
                        v2p(h,k) = v2p(h,k) + p2v(i,k);
                    end
                end
            end
        end    
    end
    
    % Mise à jour de la matrice de noeuds de parité vers noeuds de parité
    for a = 1:s(1)
        for b = 1:s(2)
            if ~isnan(p2v(a,b))
                produit = 1;
                for i = 1:s(2)
                    if (H(a,i) == 1) && (i ~= b)
                        produit = produit* tanh(v2p(a,i)/2);
                    end
                end
                p2v(a,b) = 2*atanh(produit);
            end
        end
    end
              
    %disp('fin tp')
    if  j == nit/2
        disp(v2p)
        disp(p2v)
    end
end

disp(v2p)
disp(p2v)
