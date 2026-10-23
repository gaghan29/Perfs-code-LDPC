clc;
close all;
clear;

disp('debut tp')

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

Lx = zeros(1,s(2));

disp(v2p)
disp(Lx)

% Mise à jour de la matrice de noeuds de variable vers noeuds de parité
for k = 1:length(Lx)
    for h = 1:s(1)
        if v2p(h,k) == 0
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

disp(v2p)