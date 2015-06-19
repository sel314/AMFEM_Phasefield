clc
close all
addpath('~/2DFEM-Mfiles')
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
Weights = [1.0, 1.0];       %Weights for integration
ngp = 2;                    %number of integration points in each direction
ne = Nx*Ny;                 %Total # of elements
np = (Nx+1)*(Ny+1);         %Total # of nodes
coords = zeros(np,2);       %Coordinates of the node point
elenorms = zeros(ne, 2);    %element norm vectors
Fe = zeros(4,1);
M = zeros(np,1);            %Lump mass matrix
RB = zeros(Ny,1);           %Holds elements on the very right boundary
BB = zeros(Nx,1);
elenodes = zeros(ne,4);     %Connectivity matrix
phi = zeros(np,1);          %dof vector
Kxx = 10;                   %Thermal Conductivity
Tbc = 1.0;                  %BC for temperature
rho = 1.0;                  %material density
cp = 1.0;                   %specific heat
qfluxr = -2.0;                 %Flux boundary condition
qfluxb = 5.0;

%Simulation parameters;
dt = 1/2*min(dx*dx,dy*dy)/Kxx;
TotalTime = 0.5;
NumSteps = ceil(TotalTime/dt);


%Sparse storage
Ip = zeros(ne*16,1);
Jp = zeros(ne*16,1);
S = zeros(ne*16,1);
globnodes = zeros(4,1);


%contour meshing
xx = 0:dx:1;
[X,Y] = meshgrid(xx);

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

phivec = zeros(4,1);
[phi] = SetLeftBc(phi, Nx, Ny, Tbc);
phinew = zeros(np,1);
phiold = phi;
phiMat = zeros(Nx+1,Ny+1);
Kxx = 1.0;

%--------------------------------------------------------------------------
%               Loop through and find boundary elements
%--------------------------------------------------------------------------

%Right edge boundary
for ii = 1:Ny
    RB(ii) = ii*Nx;
end%end for(ie)

%Bottom Edge Boundary
for ii = 1:Nx
    BB(ii) = ii;
end%end for(ii)

%--------------------------------------------------------------------------
%               Precalculate stiffness matrix
%--------------------------------------------------------------------------
jcount = 0;
for ie = 1:ne
    %Reset local stiffness matrix
    Ke = zeros(4,4);
    
    %Grab global node information
    for nn = 1:4
        globnodes(nn) = elenodes(ie,nn);
    end%end for(nn)
    
    %Loop through integration points
    for jj = 1:ngp
        I1 = Integration(jj);
        W1 = Weights(jj);
        for kk = 1:ngp
            I2 = Integration(kk);
            W2 = Weights(kk);
            
            %Calculate shape functions
            [B,N,jacobian] = CalcShapeFunc(ie, I1, I2, coords, elenodes);
            
            %Assemble stiffness matrix
            Ke = Ke + B'*B*Kxx*W1*W2*det(jacobian);
            
        end%end for(kk)
    end%end for(jj)
    
    %Sparse storage for stiffness matrix
    for m = 1:4
        for ii = 1:4
            jcount = jcount + 1;
            Ip(jcount) = globnodes(ii);
            Jp(jcount) = globnodes(m);
            S(jcount) = Ke(ii,m);
        end%end for(ii)
    end%end for(m)
end%end for(ie)

%Create sparse stiffness matrix
K = sparse(Ip,Jp, S, np, np);

%----------------------------------------------------------------------
%           Create flux boundary force
%----------------------------------------------------------------------

%Right edge flux
for jj = 1:Ny
    Fe = zeros(4,1);                %Element force vector
    ie = RB(jj);                    %Element on right boundary
    %Loop through integration points
    for ii = 1:ngp
        I1 = Integration(ii);
        W1 = Weights(ii);
        [B,N,jacobian] = CalcShapeFunc(ie, 1, I2, coords, elenodes);
        
        J = jacobian(2,2);          %dy/de
        
        Fe = Fe + N'*qfluxr*W1*J;    %line integral for flux
    end%end for(ii)
    %Scatter
    [Fflux] = AssembleForce(ie, coords, elenodes, Fflux, Fe);
end%end for(ie)

%Bottom edge flux
for jj = 1:Ny
    Fe = zeros(4,1);                %Element force vector
    ie = BB(jj);                    %Element on right boundary
    %Loop through integration points
    for ii = 1:ngp
        I1 = Integration(ii);
        W1 = Weights(ii);
        [B,N,jacobian] = CalcShapeFunc(ie, I1, -1, coords, elenodes);
        
        J = jacobian(1,1);          %dy/de
        
        Fe = Fe + N'*qfluxb*W1*J;    %line integral for flux
    end%end for(ii)
    %Scatter
    [Fflux] = AssembleForce(ie, coords, elenodes, Fflux, Fe);
end%end for(ie)

%Time simulation loop
for t = 1:NumSteps
    tic
    
    Fint = K*phiold;                    %Internal Force vector
    
    
    %---------------------------------------------------------------------
    %           Solve for new temperatures
    %---------------------------------------------------------------------
    
    %Solve for new temperatures
    phinew = phiold - dt*(Fint+Fflux)./M;
    
    toc
    %Set Essential boundary condition
    [phinew] = SetLeftBc(phinew, Nx, Ny, Tbc);
    
    %update
    phiold = phinew;
    
    %Store for visualization
    for jj = 1:Ny+1
        for ii = 1:Nx+1
            phiMat(jj,ii) = phinew(ii+(jj-1)*(Nx+1));
        end%end for(ii)
    end%end for(ii)
    
    if mod(t, 200) == 0
        image(phiMat*100)
        colormap(jet)
        hold on
        %colorbar
        pause(0.01)
    end%end if
    
    
end%end for(t)
