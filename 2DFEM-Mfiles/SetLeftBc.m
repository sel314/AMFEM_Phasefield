function [phi] = SetLeftBc(phi, Nx, Ny, T)
%Specific BC's for Heat problem with boundaries on all edges
%Set initial conditions, 'T' on all the boundaries zero everywhere else
for ii = 1:Nx+1
    
    %Set temperature along the left edge
    phi((ii-1)*(Nx+1)+1) = T;
end%end for(ii)