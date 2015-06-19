%SUBROUTINE TO ASSEMBLE ELEMENT STIFFNESS MATRIX INTO GLOBAL STIFFNESS
%MATRIX FOR A 4-Node QUAD
function [M] = AssembleMassVec(ie, nodesinele, M, Me)
% nodesinele: connectivity matrix
% K: Global stiffness matrix
% Kstiff: Element stiffnes matrix
nodes = zeros([4,1]);
for ii = 1:4
    nodes(ii,1) = nodesinele(ie,ii);
    %nodes(ii-1,1) = nodesinele(ie,ii);
end %end for(i)

for i = 1:4
    p = 2*nodes(i,1);
    pmo = p-1;
    M(p) = sum(Me(2*i,:))+M(p);
    M(pmo) = sum(Me(2*i-1,:)) + M(pmo);
end%end for(i)
%-------------------END SUBROUTINE---------------------------------
