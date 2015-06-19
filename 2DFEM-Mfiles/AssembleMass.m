%SUBROUTINE TO ASSEMBLE ELEMENT STIFFNESS MATRIX INTO GLOBAL STIFFNESS
%MATRIX FOR A 4-Node QUAD
function [M] = AssembleMass(ie, nodesinele, M, Me)
% nodesinele: connectivity matrix
% K: Global stiffness matrix
% Kstiff: Element stiffnes matrix
nodes = zeros([4,1]);
for ii = 1:4
    nodes(ii,1) = nodesinele(ie,ii);
end %end for(i)

for i = 1:4
    p = nodes(i,1);
    M(p) = sum(Me(i,:))+M(p);
end%end for(i)
%-------------------END SUBROUTINE---------------------------------
