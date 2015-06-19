clc
close all
N = [5, 3, 1, 4];
B = [10 20 13 5; 20 9 13 8];
f = [10 13];
phi = [3 1 5 4];
tvec = zeros(4,1);
for I = 1:4
    for i = 1:2
        for J = 1:4
                tvec(I) = tvec(I) + N(I)*f(i)*phi(J)*B(i,J);
        end%end for(J)
    end%end for(i)
end%end for(I)
tvec
N'*f*B*phi'

tvec1 = zeros(4,1);

for I = 1:4
  for i = 1:2
      for k = 1:2
          for J = 1:4
              tvec1(I) = f(i)*B(i,I)*f(k)*B(k,J)*phi(J) + tvec1(I);
          end%end for(J)
      end%end for(k)
  end%end for(i)
end%end for(i)
tvec1
B'*f'*f*B*phi'

tvec2 = zeros(4,1);
for I = 1:4
    for i = 1:2
        for J = 1:4
            tvec2(I) = tvec2(I) + N(I)*B(i,J)*f(i)*phi(J);
        end%end for(J)
    end%end for(i)
end%end for(I)
tvec2
N'*f*B*phi'

tvec3 = zeros(4,1);
for I = 1:4
    for i = 1:2
        for j = 1:2
            for J = 1:4
                tvec3(I) = tvec3(I) + f(i)*B(i,I)*B(j,J)*f(j)*phi(J);
            end%end for(J)
        end%end for(j
    end%end for(i)
end%end for(I)
tvec3
B'*f'*f*B*phi'