clc
close all
addpath('~/2DFEM-Mfiles')
%2-D implementation of Conservative level set method

%Variables defining domain of problem
X = 1;                      %length in X direction
Y = 1;                      %length in Y direction
Nx = 50;                     %# of elements in x direction
Ny = 50;                     %# of elements in y direction
dx = X/(Nx);                %grid space in x
dy = Y/(Ny);                %grid space in y

%FEM parameters and matrices
Integration = [-0.57735026, 0.57735026];   %Integration points
Weights = [1.0, 1.0];        %Weights for integration
ngp = 2;                     %number of integration points in each direction
ne = Nx*Ny;                 %Total # of elements
np = (Nx+1)*(Ny+1);         %Total # of nodes
coords = zeros(np,2);       %Coordinates of the node point
norm_ele = zeros(4,2);      %element norm vectors for each of its nodes
SurfEle= zeros(ne,1);    %Holds elements that contain surface
Fe = zeros(4,1);
phivec = zeros(4,1);
phiSTvec = zeros(4,1);
M = zeros(np,1);           %Lump mass matrix
elenodes = zeros(ne,4);     %Connectivity matrix
phi = zeros(np,1);          %dof vector
Kxx = 10;                     %Thermal Conductivity
Tbc = 1.0;                  %BC for temperature
rho = 1.0;                  %material density
cp = 1.0;                   %specific heat
it = 100;                   %# of steps to reinitialize

%Simulation parameters;
dt = 1/2*min(dx*dx,dy*dy)/Kxx;
TotalTime = 5.0;
NumSteps = ceil(TotalTime/dt);
epsil = dx;

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

%Initialize Hyperbolic Level set
x0=0.5;
y0=0.5;
R=.25;
phic=((coords(:,1)-x0).*(coords(:,1)-x0)+(coords(:,2)-y0).*(coords(:,2)-y0)).^.5-R;
phiold = 0.5*(tanh(phic/(2*epsil))+1); %phih  = 0.5 is interface
%Matrix for contour plotting
phiMat = zeros(Nx+1, Ny+1);
phinew = phiold;

%***********Precalculate local stiffness***********************************
%Precalculate shape functions
Ke = zeros(4,4);
Ce = zeros(4,4);
u =  [0.0, 100.0];
for jj = 1:ngp
    I1 = Integration(jj);
    W1 = Weights(jj);
    for kk = 1:ngp
        I2 = Integration(kk);
        W2 = Weights(kk);
        [B,N, jacobian] = CalcShapeFunc(1,I1, I2, coords, elenodes);
        Ke = Ke + transpose(B)*B*W1*W2*det(jacobian);
        %For temporary speed up
        Ce = Ce + transpose(N)*u*B*W1*W2*det(jacobian);
    end%end for(kk)
end%end for(jj)

%Precalculate nodal B matrices
[B1,N,jacobian] = CalcShapeFunc(1, -1, -1, coords, elenodes);
[B2,N,jacobian] = CalcShapeFunc(1, 1, -1, coords, elenodes);
[B3,N,jacobian] = CalcShapeFunc(1, 1, 1, coords, elenodes);
[B4,N,jacobian] = CalcShapeFunc(1, -1, 1, coords, elenodes);
%       TEST
[B0,N0, jacobian] = CalcShapeFunc(1250,0, 0, coords, elenodes);

tau = 1*dx/(norm(u));

%***********Time simulation loop******************************************
for t = 1:NumSteps
    tic
    F = zeros(np,1);
    %-----------ADVECTION STEP----------------------------------------_%
    for ie = 1:ne
        Fe = zeros(4,1);
        for nn = 1:4
            phivec(nn,1) = phiold(elenodes(ie,nn),1);
        end%end for(nn)
        for jj = 1:ngp
            I1 = Integration(jj);
            W1 = Weights(jj);
            for kk = 1:ngp
                I2 = Integration(kk);
                W2 = Weights(kk);
                [B,N, jacobian] = CalcShapeFunc(ie,I1, I2, coords, elenodes);
                Fe = Fe + transpose(N)*u*B*phivec*W1*W2*det(jacobian) ...
                     + tau*transpose(B)*u'*u*B*W1*W2*det(jacobian)*phivec;  %upwind
                F = AssembleForce(ie,coords, elenodes, F, Fe);
            end%end for(kk)
        end%end for(jj)
    end%end for(ie)
    for ii = 1:np
        phinew(ii,1) = phiold(ii,1) - dt*F(ii,1)/M(ii);
    end%end for(ii)

    %--------------------------------------------------------------------%
    toc
    %
    max(phinew)
    %stores results in matrix
    phiold = phinew;
    
%******************* REINITIALIZATION STEP*******************************%

%     Calculate norms from time step !! AT ELEMENT CENTER
%     ict = 0;
%     %Find elements that contain interface
%     for ie = 1:ne
%         for nn = 1:4
%             phivec(nn,1) = phinew(elenodes(ie,nn),1);
%         end%end for(nn)
%         phicenter = N0*phivec;  %Calculate averaged level set
%         if phicenter < 0.7 && phicenter > 0.3 
%             ict = ict + 1;
%             SurfEle(ict,1) = ie;
%         end%end if
%     end%end for(ii)
%     phistar = phinew;
%     for tsteps = 1:it
%         F = zeros(np,1);
%         for ie = 1:ict
%             Fe = zeros(4,1);
%             ele = SurfEle(ie,1);
%             %Calculate norm vector for each node within element 'ie'
%             for nn = 1:4
%                 %Find corresponding nodal phistar
%                 phiSTvec(nn,1) = phistar(elenodes(ele,nn),1);
%             end%end for(nn)
%             
%             %Grab nodal phi
%             for mm = 1:4
%                 phivec(nn,1) = phiold(elenodes(ele,nn),1);
%             end%end for(nn)
%             
%             norm_ele(1,:) = (B1*phiSTvec)'./(norm(B1*phiSTvec));
%             norm_ele(2,:) = (B2*phiSTvec)'./(norm(B2*phiSTvec));
%             norm_ele(3,:) = (B3*phiSTvec)'./(norm(B3*phiSTvec));
%             norm_ele(4,:) = (B4*phiSTvec)'./(norm(B4*phiSTvec));
%             
%             %Loop over integration points
%             for jj = 1:ngp
%                 I1 = Integration(jj);
%                 W1 = Weights(jj);
%                 for kk = 1:ngp
%                     I2 = Integration(kk);
%                     W2 = Weights(kk);
%                     [B,N, jacobian] = CalcShapeFunc(ele,I1, I2, coords, elenodes);
%                     %Calculate integrated level set
%                     fbar = (1-phi_int)*norm_ele;
%                     phi_int = N*phivec;
%                     T1 = N'*N*fbar*B*phivec;
%                     T2 = N'*N*B'*fbar'*phivec;
%                     T3 = epsil*B'*B*(norm_ele*norm_ele')*phivec;
%                     Fe = Fe + (T1+T2+T3)*W1*W2*det(jacobian);
%                     F = AssembleForce(ele,coords, elenodes, F, Fe);
%                 end%end for(kk)
%             end%end for(jj)
%         end%end for(ie)
%         for ii = 1:np
%             phinew(ii,1) = phiold(ii,1) - dt*F(ii,1)/M(ii);
%         end%end for(ii)
%         phiold = phinew;
%     end%end for(it)
%************************ END REINITIALIZATION STEP ***********************
   
    %----------------Plot the interface----------------------------------%
    for jj = 1:Ny+1
        for ii = 1:Nx+1
           phiMat(ii,jj) = phinew(ii+(jj-1)*(Nx+1)); 
        end%end for(ii)
    end%end for(ii)
    contour(X, Y, phiMat, [0.5, 0.5], 'r', 'linewidth', 2)
    pause(0.1)
end%end for(t)
