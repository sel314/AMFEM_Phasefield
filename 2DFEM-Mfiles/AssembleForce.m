%SUBROUTINE TO ASSEMBLE ELEMENT FORCE VECTOR INTO GLOBAL FORCE
%FOR A 4-NODE QUAD
function [F] = AssembleForce(ie, coords, nodesinele, F, Fe)
% nodesinele: connectivity matrix
% F: GLOBAL FORCE
% Fe: LOCAL FORCE
nodes = zeros([4,1]);
for ii = 1:4
    nodes(ii,1) = nodesinele(ie,ii);
end %end for(i)

for i = 1:4
    p = nodes(i,1);
    F(p,1) = Fe(i,1)+F(p,1);
end%end for(i)
%-------------------END SUBROUTINE---------------------------------
