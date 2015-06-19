function [phi] = SetBc(phi, Nx, Ny, T)
%Specific BC's for Heat problem with boundaries on all edges
%Set initial conditions, 'T' on all the boundaries zero everywhere else
for ii = 1:Nx+1
    
    %Set temperature along the bottom edge
    phi(ii) = T;
    
    %Set temperature along the right edge
    phi(ii*(Nx+1)) = T;
    
    %Set temperature along the left edge
    phi((ii-1)*(Nx+1)+1) = T;
    
    %Set temperature along the top edge
    phi((Nx+1)*Ny + ii) = T;
    
end%end for(ii)