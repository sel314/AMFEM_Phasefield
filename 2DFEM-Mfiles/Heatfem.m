function Heatfem(Nx, Ny, X, Y)
%implement convection BC
%finite element heat code for a 2D square
%Nx = number of elements in X direction
%Ny = number of elements in Y direction
%X = length of mesh in X direction
%Y = length of mesh in Y direction
clc
kxx = 1;
kyy = 1;
Temp_L = 300; %ebc of lower boundary
Temp_U = 500; %ebc of upper boundary
ne = Nx*Ny; %number of elements
np = (Nx+1)*(Ny+1); %number of nodal points
coords = zeros([np,2]); %initialize coordinates array
K = zeros([np, np]);%initialize Global stiffness array
M = zeros([np,np]);
Stiffness = zeros([4,4]);%initialize Local Stiffness Array
Load = zeros([np, 1]); %intialize Load vector array
nodesinele = zeros([Nx*Ny*4, 1]); %initialize element node numbering
incx = X/Nx;
incy = Y/Ny;
coordcount = 0;
y = 0;
numboundnodes = (Ny+1)*2;
boundarynodes = zeros([numboundnodes,1]);
numdofnodes = np - numboundnodes;
dofnodes = zeros([numdofnodes,1]);
%K = [kxx 0; 0 kyy];
%generate the mesh
for i = 1:Nx+1
    x = 0;
    for j = 1:Ny+1
        coords(coordcount*(Nx+1)+j,1) = x;
        coords(coordcount*(Ny+1)+j,2) = y;
        x = x+incx;
    end%end for(j)
    y = y + incy;
    coordcount = coordcount + 1; 
end %end for(i)

%find the number of essential boundary nodes
boundcount = 1;
%counter to find dof nodes
dofcount = 1;
for i = 1:np
    if coords(i,2) == 0
        boundarynodes(boundcount) = i;
        boundcount = boundcount + 1;
    elseif abs(coords(i,2)-Y) < 1e-14
        boundarynodes(boundcount) = i;
        boundcount = boundcount + 1;
    else
        dofnodes(dofcount) = i;
        dofcount = dofcount + 1;
    end%end if
end%end for(i)
boundcount = boundcount - 1;

%generate the nodal number of each finite element
count = 1;
inc = 0;
for i = 1:Ny
    for j = 1:Nx
        nodesinele(count,1) = j+inc*Ny+inc;
        count = count +1;
        nodesinele(count,1) = j+1+inc*Ny+inc;
        count = count +1;
        nodesinele(count,1) = j+Ny+2+inc*Ny+inc;
        count = count +1;
        nodesinele(count,1) = j+Ny+1+inc*Ny+inc;
        count = count +1;
    end%end for(j)
    inc = inc+1;
end%end for(i)
%Load = sparse(np,1)
%implement Boundary conditions
for i = 1:boundcount
    a = boundarynodes(i);
    if coords(a,2) == 0
        Load(a) = Temp_L;
    elseif abs(coords(a,2)-Y) < 1e-14
        Load(a) = Temp_U;
    end%end if
end%end for(i)

%K = sparse(np,np)

%Perform finite element routine
%Loop over each element, calculate stiffness matrix and assemble
for i = 1:ne
    [Kstiff, Mstiff] = elementstiff(i, np, coords, nodesinele, K);
    [K, M] = assemble(i, np, coords, nodesinele, K, Kstiff, M, Mstiff);
end%end for(i)

%PERFORM ROW AND COLUMN ADJUSTMENT METHOD
%MODIFY THE FORCE VECTOR
for i = 1:numdofnodes
    for j = 1:boundcount
        Load(dofnodes(i)) = Load(dofnodes(i))-K(boundarynodes(j),dofnodes(i))*Load(boundarynodes(j));
    end%end for(j)
end%end for(i)

%ZERO OUT ESSENTIAL BC IN Stiffness matrix
for i = 1: boundcount
    k = boundarynodes(i);
    for j = 1:np
        K(k,j) = 0;
        K(j,k) = 0;
    end %end for(j)
    K(k,k) = 1;
end%end for(i)
M;
Temp = inv(K)*Load



%%%%%%%---------------------SUBROUTINE DECLARATIONS------------%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%SUBROUTINE TO CALCULATE PARTIAL DERIVATIVES OF SHAPE FUNCTION W.R.T GLOBAL
%COORDINATES
function [B, N, jacobian] = dN(ie, z, e, np, coords, nodesinele)
dcoords = zeros([4,2]);
jacobian = zeros([2,2]);
start = ie*4 - 3;
for i = 1:4
    dcoords(i,1) = coords(nodesinele(start, 1),1);
    dcoords(i,2) = coords(nodesinele(start, 1),2);
    start = start + 1;
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

%SUBROUTINE TO CALCULATE THE ELEMENT STIFFNESS MATRIX W.R.T GLOBAL
%COORDINATES
function [Kstiff, Mstiff] = elementstiff(ie, np, coords, nodesinele, K)
Integration = [-0.57735026, 0.57735026];
Weights = [1.0, 1.0];
it = size(Integration);
Kstiff = 0;
Mstiff = 0;
for i = 1:it(2)
    I1 = Integration(i);
    W1 = Weights(i);
    for j = 1: it(2)
        I2 = Integration(i);
        W2 = Weights(i);
        [B, N, jacobian] = dN(ie, I1, I2, np, coords, nodesinele);
        Kstiff = Kstiff + transpose(B)*B*W1*W2*det(jacobian);
        Mstiff = Mstiff + transpose(N)*N*W1*W2*det(jacobian);
    end%end for(j)
end%End for(i)

%SUBROUTINE TO ASSEMBLE ELEMENT STIFFNESS MATRIX INTO GLOBAL STIFFNESS
%MATRIX
function [K, M] = assemble(ie, np, coords, nodesinele, K, Kstiff, M, Mstiff)
start = ie*4-3;
nodes = zeros([4,1]);
for i = 1:4
    nodes(i,1) = nodesinele(start,1);
    start = start + 1;
end %end for(i)

for i = 1:4
    p = nodes(i,1);
    for j = 1:4
        l = nodes(j,1);
        K(p,l) = Kstiff(i,j)+K(p,l);
        M(p,l) = Mstiff(i,j)+M(p,l);
    end%end for(j)
end%end for(i)
%-------------------END SUBROUTINE---------------------------------
