%SUBROUTINE TO ASSEMBLE ELEMENT STIFFNESS MATRIX INTO GLOBAL STIFFNESS
%MATRIX FOR A 4-Node QUAD
function [K] = Assemble(ie, coords, nodesinele, K, Ke)
% nodesinele: connectivity matrix
% K: Global stiffness matrix
% Kstiff: Element stiffnes matrix
nodes = zeros([4,1]);
for ii = 1:4
    nodes(ii,1) = nodesinele(ie,ii);
end %end for(i)

for i = 1:4
    p = nodes(i,1);
    for j = 1:4
        l = nodes(j,1);
        K(p,l) = Ke(i,j)+K(p,l);
    end%end for(j)
end%end for(i)
%-------------------END SUBROUTINE---------------------------------
