%SUBROUTINE TO ASSEMBLE ELEMENT FORCE VECTOR INTO GLOBAL FORCE
%FOR A VECTOR SOLUTION
%FOR A 4-NODE QUAD
function [F] = AssembleForceVec(ie, coords, nodesinele, F, Fe)
% nodesinele: connectivity matrix
% F: GLOBAL FORCE
% Fe: LOCAL FORCE
nodes = zeros([4,1]);
for ii = 1:4
    nodes(ii,1) = nodesinele(ie,ii);
end %end for(i)

for i = 1:4
    p = 2*nodes(i,1);
    pmo = 2*nodes(i,1)-1;
    F(p,1) = Fe(2*i,1) + F(p,1);
    F(pmo,1) = Fe(2*i-1,1) + F(pmo,1);
    
end%end for(i)
%-------------------END SUBROUTINE---------------------------------
