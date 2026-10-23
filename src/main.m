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

for k = 1:length(Lx)
    
end

