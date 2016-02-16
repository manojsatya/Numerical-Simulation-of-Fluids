function [solU,solV] = postProcess(x,nx,ny)

uGrid = 1:(nx+1)*(ny+2);
vGrid = (nx+1)*(ny+2) + 1 : (nx+1) * (ny+2) + (nx+2)*(ny+1);
% pGrid = (nx+1)*(ny+2) + 1 + (nx+2)*(ny+1): (nx+1) * (ny+2) + (nx+2)*(ny+1) + (nx+2)*(ny+2);

U = x(uGrid);
V = x(vGrid);
% P = x(pGrid);

U = reshape(U, [nx+1, ny+2])'; 
V = reshape(V, [nx+2, ny+1])'; 
% P = reshape(P, [nx+2, ny+2])';

solU = zeros(ny,nx);
solV = zeros(ny,nx);
% solP = zeros(ny,nx);

for i = 1:nx
 for j = 1:ny
  solU(j,i) = 0.5* ( U(j+1,i) + U(j+1,i+1) );
 end
end

for i = 1:nx
 for j = 1:ny
  solV(j,i) = 0.5* ( V(j+1,i) + V(j+1,i+1) );
 end
end

end