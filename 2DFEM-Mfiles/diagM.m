%--------Subroutine to diagonalize element mass matrix--------
function [Me] = diagM(Me, ind)
for i = 1:ind
    MeD = 0; %summed up row of Me
    for j = 1:ind
        MeD = MeD + Me(i,j);
        Me(i,j) = 0;
    end%end for(j)
    Me(i,i) = MeD;
end%end for(i)