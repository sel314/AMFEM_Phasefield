%SUBROUTINE TO CALCULATE PARTIAL DERIVATIVES OF SHAPE FUNCTION W.R.T GLOBAL
%COORDINATES
%           e
%           |
%           |
%           |
%   ----------------->z
%           |
%           |
%           |
function [B, N, jacobian] = CalcShapeFunc(ie, z, e, coords, nodesinele)
dcoords = zeros([4,2]);
jacobian = zeros([2,2]);
for i = 1:4
    dcoords(i,1) = coords(nodesinele(ie, i),1);
    dcoords(i,2) = coords(nodesinele(ie, i),2);
end%end for(i)

%declare local shape functions
Shapefunc = zeros([1,4]);
Shapefunc(1) = 0.25*(1-z)*(1-e);
Shapefunc(2) = 0.25*(1+z)*(1-e);
Shapefunc(3) = 0.25*(1+z)*(1+e);
Shapefunc(4) = 0.25*(1-z)*(1+e);
N = Shapefunc;


%declare partial derivatives of local shape functions
dShapefunc = zeros([2,4]);
dShapefunc(1,1) = -0.25*(1-e);
dShapefunc(2,1) = -0.25*(1-z);
dShapefunc(1,2) = 0.25*(1-e);
dShapefunc(2,2) = -0.25*(1+z);
dShapefunc(1,3) = 0.25*(1+e);
dShapefunc(2,3) = 0.25*(1+z);
dShapefunc(1,4) = -0.25*(1+e);
dShapefunc(2,4) = 0.25*(1-z);

%calculate the jacobian:
jacobian = dShapefunc*dcoords;
Jinv = inv(jacobian);
B = Jinv*dShapefunc;
%---------------------END SUBROUTINE-------------------------------
