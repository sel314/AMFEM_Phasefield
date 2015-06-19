function [norm_nodes] = GetNormals(np, ne, phiold, Integration, Weights, ...
                        ngp, norm_nodes, elenodes, coords, Mv)
addpath('~/2DFEM-Mfiles')
tolr = 1e-7;
phivec = zeros(4,1);
%Go through and calculate normals at the nodes
Fv = zeros(2*np,1);	    %Global force for vector solutions
for ie = 1:ne
    FeV = zeros(8,1);
    %Grab level set values at the nodes
    for nn = 1:4
        phivec(nn,1) = phiold(elenodes(ie,nn),1);
    end%end for(nn)
    %Loop through integration points
    for jj = 1:ngp
        I1 = Integration(jj);
        W1 = Weights(jj);
        for kk = 1:ngp
            I2 = Integration(kk);
            W2 = Weights(kk);
            [B,N, jacobian] = CalcShapeFunc(ie,I1, I2, coords, elenodes);
            %construct vector test function
            nlize = norm(B*phivec);
            %Calculate normalized gradient vector
            if nlize <= tolr
                nvec = [0;0];
            else
                nvec = (B*phivec)./nlize;
            end%end if
            Nv(1,:) = [N(1), 0 N(2), 0, N(3), 0, N(4), 0];
            Nv(2,:) = [0, N(1), 0, N(2), 0, N(3), 0, N(4)];
            FeV = FeV + Nv'*nvec*W1*W2*det(jacobian);
        end%end for(kk)
    end%end for(jj)
    Fv = AssembleForceVec(ie, coords, elenodes, Fv, FeV);
end%end for(ie)

for ii = 1:np*2
    norm_nodes(ii,1) = Fv(ii)/Mv(ii);
end%end for(ii)