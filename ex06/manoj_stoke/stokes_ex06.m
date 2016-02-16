clc
clear all


nx = 20; % Number of nodes in x direction
ny = 15; % Number of nodes in x direction

xlength = 1; % Domain length x-direction
ylength = 1; % Domain length y-direction

hx = xlength / nx; %grid size x-direction
hy = ylength / ny; %grid size x-direction

u_top = 1; % velocity top wall

tol = 1e-6; % Tolerance for gmres iteration

%%%%%%%%%%%%%%%%%%%%%%%%%%%%--------------Creating Matrix M---------%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 %%% ------Laplacian of u and v--------%%%%%%%%%
ex = (1/(hx*hx))* ones(nx+2,1); % add ghost layers
Uxx = spdiags([ex -2*ex ex], [-1 0 1], nx+2, nx+2);
ey = (1/(hy*hy))* ones(ny+2,1);% add ghost layers
Uyy = spdiags([ey, -2*ey ey], [-1 0 1], ny+2, ny+2);
Lu = -(kron(Uyy, speye(nx+2)) + kron(speye(ny+2), Uxx)); % Laplacian of u
Lv = -(kron(Uyy, speye(nx+2)) + kron(speye(ny+2), Uxx)); % Laplacian of v

 %%%%%%--------Gradient of p in x ---%%%%%%%
coeff_x = (1/(2*hx))* ones(nx+2,1);
Dx = spdiags([-1 * coeff_x 0*coeff_x 1*coeff_x], [-1 0 1], nx+2, nx+2);
grad_px = kron(speye(ny+2), Dx);

%%%%%%--------Gradient of p in y ---%%%%%%%
coeff_y = (1/(2*hy))* ones(ny+2,1);
Dy = spdiags([-1 * coeff_y 0*coeff_y 1*coeff_y], [-1 0 1], ny+2, ny+2);
grad_py = kron(Dy, speye(nx+2));

null_Matrix = zeros((nx+2)*(ny+2));

M_init = [Lu null_Matrix grad_px;null_Matrix Lv grad_py;grad_px' grad_py' null_Matrix];   


M = correction(M_init,nx,ny); % Editing matrix M_init to remove rows and columns not required in the grid
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


dlmwrite('Matrix_M.txt',full(M));


I = speye(size(M,1)); %%%%% For dirichlet conditions, Diagonal has to be one

[M,I] = setBoundaryU(M,I,nx,ny); 
M = setBoundaryV(M,I,nx,ny);
M = setBoundaryP(M,I,nx,ny);

rhs = zeros(size(M,1),1);

lidVelocity = (nx+1)*(ny+2) -(nx+1)+1 : (nx+1)*(ny+2);

rhs(lidVelocity) = u_top;

iterations = size(rhs,1);

[x,flag,relres,iter] = gmres(M,rhs,[],tol,iterations);


if(flag ==0)
    fprintf('gmres converged to the desired tolerance and iterations with residual %f\n',relres);
else
    fprintf('Solution did not converge');
end

[solU,solV] = postProcess(x,nx,ny);

quiver(solU,solV);
xlabel('xlength');
ylabel('ylength');

