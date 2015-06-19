clc
close all
addpath('/Users/Stephen Lin/Documents/2DFEM-Mfiles')
%2-D implementation of Conservative level set method

%Variables defining domain of problem
X = 1;                      %length in X direction
Y = 1;                      %length in Y direction
Nx = 100;                     %# of elements in x direction
Ny = 100;                     %# of elements in y direction
dx = X/(Nx);                %grid space in x
dy = Y/(Ny);                %grid space in y

%FEM parameters and matrices
Integration = [-0.57735026, 0.57735026];   %Integration points
Weights = [1.0, 1.0];        %Weights for integration
ngp = 2;                     %number of integration points in each direction
ne = Nx*Ny;                 %Total # of elements
np = (Nx+1)*(Ny+1);         %Total # of nodes
coords = zeros(np,2);       %Coordinates of the node point
elenorms = zeros(ne, 2);    %element norm vectors
Fe = zeros(4,1);
M = zeros(np,1);           %Lump mass matrix
elenodes = zeros(ne,4);     %Connectivity matrix
phi = zeros(np,1);          %dof vector
Kxx = 10;                     %Thermal Conductivity
Tbc = 1.0;                  %BC for temperature
rho = 1.0;                  %material density
cp = 1.0;                   %specific heat

%Simulation parameters;
dt = 1/2*min(dx*dx,dy*dy)/Kxx;
TotalTime = 5.0;
NumSteps = ceil(TotalTime/dt);

%contour meshing
xx = 0:dx:1;
[X,Y] = meshgrid(xx);

%Generate a structured mesh

%Create coordinates with lexographic ordering
for jj = 1:Ny+1
    for ii = 1:Nx+1
        coords(ii+(jj-1)*(Nx+1),1) = (ii-1)*dx;
        coords(ii+(jj-1)*(Nx+1),2) = (jj-1)*dy;
    end%end for(ii)
end%end for(jj)

%Generate Connectivity matrix
for jj = 1:Ny
    for ii = 1:Nx
        elenodes(ii+(jj-1)*Nx, 1) = ii + (jj-1)*(Nx+1);
        elenodes(ii+(jj-1)*Nx, 2) = ii + (jj-1)*(Nx+1) + 1;
        elenodes(ii+(jj-1)*Nx, 3) = ii + jj*(Nx+1) + 1;
        elenodes(ii+(jj-1)*Nx, 4) = ii + jj*(Nx+1);
    end%end for(ii)
end%end for(jj)

phinew = zeros(np,1);
phiold = phi;
%Generate Mass Matrix
display('Assembling Mass Matrix')
for ie = 1:ne
    %integration loops
    Me = zeros(4,4);
    for jj = 1:ngp
        I1 = Integration(jj);
        W1 = Weights(jj);
        for kk = 1:ngp
            I2 = Integration(kk);
            W2 = Weights(kk);
            [B,N,jacobian] = CalcShapeFunc(ie, I1, I2, coords, elenodes);
            Me = Me + rho*cp*transpose(N)*N*W1*W2*det(jacobian);
            [Me] = diagM(Me, 4);
            [M] = AssembleMass(ie, elenodes, M, Me);
        end%end for(kk)
    end%end for(jj)
end%end for(ie)
M = sparse(M);
addpath('/Users/Stephen Lin/Documents/2DFEM-Mfiles')
phivec = zeros(4,1);
[phi] = SetBc(phi, Nx, Ny, Tbc);
phinew = zeros(np,1);
phiold = phi;
phiMat = zeros(Nx+1,Ny+1);
Kxx = 1.0;

%Precalculate shape functions
Ke = zeros(4,4);
for jj = 1:ngp
    I1 = Integration(jj);
    W1 = Weights(jj);
    for kk = 1:ngp
        I2 = Integration(kk);
        W2 = Weights(kk);
        [B,N, jacobian] = CalcShapeFunc(1,I1, I2, coords, elenodes);
        Ke = Ke + transpose(B)*B*Kxx*W1*W2*det(jacobian);
    end%end for(kk)
end%end for(jj)

%Time simulation loop
for t = 1:NumSteps
    tic
    %assemble stiffness matrix
    F = zeros(np,1);
    for ie = 1:ne
        Fe = zeros(4,1);
        for nn = 1:4
            phivec(nn,1) = phiold(elenodes(ie,nn),1);
        end%end for(nn)
        I2 = Integration(kk);
        W2 = Weights(kk);    
        Fe = Fe + Ke*phivec;
        [F] = AssembleForce(ie, coords, elenodes, F, Fe);
    end%end for(ie)
    toc
    for ii = 1:np
        phinew(ii,1) = phiold(ii,1) - dt*F(ii,1)/M(ii);
    end%end for(ii)
    %phinew = phiold - dt*(M\(F));
    [phinew] = SetBc(phinew, Nx, Ny, Tbc);
    
    %stores results in matrix
    phiold = phinew;
    for jj = 1:Ny+1
        for ii = 1:Nx+1
            phiMat(ii,jj) = phinew(ii+(jj-1)*(Nx+1));
        end%end for(ii)
    end%end for(ii)
    image(phiMat*100)
    colormap(jet)
    axis equal
    hold on
    colorbar
    pause(0.1)
end%end for(t)
