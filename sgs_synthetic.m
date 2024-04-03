t1=clock;
% Parameters and Initialization
Lx = 2; Ly = 2; Lz = 2; % Domain size
sigma = 100; % Example value, adjust as necessary
Nx = 101; Ny = 101; Nz = 101; % Number of grid points in each dimension
dx = Lx/(Nx-1);
dy = Ly/(Ny-1);
dz = Lz/(Nz-1);
x = linspace(dx/2-1, 1-dx/2, Nx);
y = linspace(dy/2-1, 1-dy/2, Ny);
z = linspace(dz/2-1, 1-dz/2, Nz);
[X, Y, Z] = meshgrid(y, z, x);
% Initialize the variables F1 to F6 on the 3D grid
F1 = zeros(Nx, Ny, Nz);
F2 = zeros(Nx, Ny, Nz);
F3 = zeros(Nx, Ny, Nz);
F4 = zeros(Nx, Ny, Nz);
F5 = zeros(Nx, Ny, Nz);
F6 = zeros(Nx, Ny, Nz);
n_iter = 0;
thres = 1e-8;

% Boundary condition setup
% Boundary function F_b
F_b = @(p, q) (abs(p) <= 0.2 & abs(q) <= 0.2);

% Apply boundary conditions
for i = 1:Nx
    for j = 1:Ny
        for k = 1:Nz
            % F1 boundary condition at x = -1
            if F_b(y(j), z(k))
                F1(1, j, k) = 1;
            else
                F1(1, j, k) = 0;
            end

            % F3 boundary condition at y = -1
            if F_b(x(i), z(k))
                F3(i, 1, k) = 1;
            else
                F3(i, 1, k) = 0;
            end

            % F5 boundary condition at z = -1
            if F_b(x(i), y(j))
                F5(i, j, 1) = 1;
            else
                F5(i, j, 1) = 0;
            end
        end
    end
end

% F2, F4, and F6 boundary conditions are set to 0 at x=1, y=1, and z=1, respectively
F2(Nx, :, :) = 0; % F2 at x = 1
F4(:, Ny, :) = 0; % F4 at y = 1
F6(:, :, Nz) = 0; % F6 at z = 1


% Residual
res = residual(F1, F2, F3, F4, F5, F6, dx, dy, dz, sigma);
res_arr = [res];
avg_inner_iter = [];

while res > thres && n_iter < 5000
    n_nonlin_eqs = 0;
    n_inner_iter = 0;

    for i=2:Nx-1
        for j=2:Ny-1
            for k=2:Nz-1
                density = F1(i,j,k) + F2(i,j,k) + F3(i,j,k) + F4(i,j,k)+ F5(i,j,k)+F6(i,j,k);
                if density == 0
                    density = 1;
                end

                resF1 = (F1(i,j,k) - F1(i-1,j,k)) / (dx) - sigma * (density / 6 - F1(i,j,k));
                resF2 = (F2(i,j,k) - F2(i+1,j,k)) / (dx) - sigma * (density / 6 - F2(i,j,k));
                resF3 = (F3(i,j,k) - F3(i,j-1,k)) / (dy) - sigma * (density / 6 - F3(i,j,k));
                resF4 = (F4(i,j,k) - F4(i,j+1,k)) / (dy) - sigma * (density / 6 - F4(i,j,k));
                resF5 = (F5(i,j,k) - F5(i,j,k-1)) / (dz) - sigma * (density / 6 - F5(i,j,k));
                resF6 = (F6(i,j,k) - F6(i,j,k+1)) / (dz) - sigma * (density / 6 - F6(i,j,k));
                res = sqrt((resF1)^2 + (resF2)^2 + (resF3)^2 + (resF4)^2 + (resF5)^2 + (resF6)^2);

                n_nonlin_eqs = n_nonlin_eqs + 1;
                while res > 1e-3 * thres
                    F1(i,j,k) = (sigma * density / 6 +F1(i-1,j,k) / dx) / (1/dx + sigma);
                    F2(i,j,k) = (sigma * density / 6 +F2(i+1,j,k) / dx) / (1/dx + sigma);
                    F3(i,j,k) = (sigma * density / 6 +F3(i,j-1,k) / dy) / (1/dy + sigma);
                    F4(i,j,k) = (sigma * density / 6 +F4(i,j+1,k) / dy) / (1/dy + sigma);
                    F5(i,j,k) = (sigma * density / 6 +F5(i,j,k-1) / dz) / (1/dz + sigma);
                    F6(i,j,k) = (sigma * density / 6 +F6(i,j,k+1) / dz) / (1/dz + sigma);

                    density = F1(i,j,k) + F2(i,j,k) + F3(i,j,k) + F4(i,j,k)+ F5(i,j,k)+F6(i,j,k);
                    if density == 0
                        density = 1;
                    end

                    resF1 = (F1(i,j,k) - F1(i-1,j,k)) / (dx) - sigma * (density / 6 - F1(i,j,k));
                    resF2 = (F2(i,j,k) - F2(i+1,j,k)) / (dx) - sigma * (density / 6 - F2(i,j,k));
                    resF3 = (F3(i,j,k) - F3(i,j-1,k)) / (dy) - sigma * (density / 6 - F3(i,j,k));
                    resF4 = (F4(i,j,k) - F4(i,j+1,k)) / (dy) - sigma * (density / 6 - F4(i,j,k));
                    resF5 = (F5(i,j,k) - F5(i,j,k-1)) / (dz) - sigma * (density / 6 - F5(i,j,k));
                    resF6 = (F6(i,j,k) - F6(i,j,k+1)) / (dz) - sigma * (density / 6 - F6(i,j,k));
                    res = sqrt((resF1)^2 + (resF2)^2 + (resF3)^2 + (resF4)^2 + (resF5)^2 + (resF6)^2);

                    n_inner_iter = n_inner_iter + 1;

                end
            end
        end
    end

    for i=Nx-1:-1:2
        for j=Ny-1:-1:2
            for k=Nz-1:-1:2
                resF1 = (F1(i,j,k) - F1(i-1,j,k)) / (dx) - sigma * (density / 6 - F1(i,j,k));
                resF2 = (F2(i,j,k) - F2(i+1,j,k)) / (dx) - sigma * (density / 6 - F2(i,j,k));
                resF3 = (F3(i,j,k) - F3(i,j-1,k)) / (dy) - sigma * (density / 6 - F3(i,j,k));
                resF4 = (F4(i,j,k) - F4(i,j+1,k)) / (dy) - sigma * (density / 6 - F4(i,j,k));
                resF5 = (F5(i,j,k) - F5(i,j,k-1)) / (dz) - sigma * (density / 6 - F5(i,j,k));
                resF6 = (F6(i,j,k) - F6(i,j,k+1)) / (dz) - sigma * (density / 6 - F6(i,j,k));
                res = sqrt((resF1)^2 + (resF2)^2 + (resF3)^2 + (resF4)^2 + (resF5)^2 + (resF6)^2);

                n_nonlin_eqs = n_nonlin_eqs + 1;
                while res > 1e-3 * thres

                    F1(i,j,k) = (sigma * density / 6 +F1(i-1,j,k) / dx) / (1/dx + sigma);
                    F2(i,j,k) = (sigma * density / 6 +F2(i+1,j,k) / dx) / (1/dx + sigma);
                    F3(i,j,k) = (sigma * density / 6 +F3(i,j-1,k) / dy) / (1/dy + sigma);
                    F4(i,j,k) = (sigma * density / 6 +F4(i,j+1,k) / dy) / (1/dy + sigma);
                    F5(i,j,k) = (sigma * density / 6 +F5(i,j,k-1) / dz) / (1/dz + sigma);
                    F6(i,j,k) = (sigma * density / 6 +F6(i,j,k+1) / dz) / (1/dz + sigma);

                    density = F1(i,j,k) + F2(i,j,k) + F3(i,j,k) + F4(i,j,k)+ F5(i,j,k)+F6(i,j,k);
                    if density == 0
                        density = 1;
                    end

                    resF1 = (F1(i,j,k) - F1(i-1,j,k)) / (dx) - sigma * (density / 6 - F1(i,j,k));
                    resF2 = (F2(i,j,k) - F2(i+1,j,k)) / (dx) - sigma * (density / 6 - F2(i,j,k));
                    resF3 = (F3(i,j,k) - F3(i,j-1,k)) / (dy) - sigma * (density / 6 - F3(i,j,k));
                    resF4 = (F4(i,j,k) - F4(i,j+1,k)) / (dy) - sigma * (density / 6 - F4(i,j,k));
                    resF5 = (F5(i,j,k) - F5(i,j,k-1)) / (dz) - sigma * (density / 6 - F5(i,j,k));
                    resF6 = (F6(i,j,k) - F6(i,j,k+1)) / (dz) - sigma * (density / 6 - F6(i,j,k));
                    res = sqrt((resF1)^2 + (resF2)^2 + (resF3)^2 + (resF4)^2 + (resF5)^2 + (resF6)^2);

                    n_inner_iter = n_inner_iter + 1;

                end
            end
        end
    end

    res = residual(F1, F2, F3, F4, F5, F6, dx, dy, dz, sigma);
    res_arr = [res_arr res];

    n_iter = n_iter + 1;
    fprintf("Iter %d: Residual: %f, Averge number of inner iterations: %f\n", n_iter, res, n_inner_iter/n_nonlin_eqs);

    % if res < thres && n_iter > 1
    %     fprintf("Iter %d: Residual: %f, Averge number of inner iterations: %f\n ", n_iter, res, n_inner_iter/n_nonlin_eqs);
    %     break; % Exit the loop if the solution has converged
    % end

    residuals(n_iter) = res;
    avg_inner_iter = [avg_inner_iter n_inner_iter/n_nonlin_eqs];

end

% Visualization
F_density = F1 + F2 + F3 + F4 + F5 + F6;
slice(X, Y, Z, F_density, 0, 0, 0); % Adjust slice positions as needed
xlabel('X-axis');
ylabel('Y-axis');
zlabel('Z-axis');
title('Density Distribution');
colorbar;
drawnow; % Use drawnow for updating the figure window

residuals = residuals(1:n_iter-1); % Trim the residuals array
%Plotting the residuals on a logarithmic scale without markers
figure;
semilogy(1:length(residuals), residuals, 'LineWidth', 2);
xlabel('Iteration Number');
ylabel('Residual (L2 Norm)');
title('Convergence of the Iterative Method');
grid on; % Adds grid for better readability
drawnow;

t2=clock;
t=etime(t2,t1);
disp(['程序运行时间为：',num2str(t),'秒']);


function res = residual(F1, F2, F3, F4, F5, F6, dx, dy, dz, sigma)
density = F1 + F2 + F3 + F4 + F5 + F6;
density = density(2:end-1,2:end-1,2:end-1);
density(density == 0) = 1;

resF1 = -(F1(2:end-1,2:end-1,2:end-1) - F1(1:end-2,2:end-1,2:end-1)) / (dx) + sigma * (density / 6 - F1(2:end-1,2:end-1,2:end-1));
resF2 = -(F2(2:end-1,2:end-1,2:end-1) - F2(3:end,2:end-1,2:end-1)) / (dx) + sigma * (density / 6 - F2(2:end-1,2:end-1,2:end-1));
resF3 = -(F3(2:end-1,2:end-1,2:end-1) - F3(2:end-1,1:end-2,2:end-1)) / (dy) + sigma * (density / 6 - F3(2:end-1,2:end-1,2:end-1));
resF4 = -(F4(2:end-1,2:end-1,2:end-1) - F4(2:end-1,3:end,2:end-1)) / (dy) + sigma * (density / 6 - F4(2:end-1,2:end-1,2:end-1));
resF5 = -(F5(2:end-1,2:end-1,2:end-1) - F5(2:end-1,2:end-1,1:end-2)) / (dz) + sigma * (density / 6 - F5(2:end-1,2:end-1,2:end-1));
resF6 = -(F6(2:end-1,2:end-1,2:end-1) - F6(2:end-1,2:end-1,3:end)) / (dz) + sigma * (density / 6 - F6(2:end-1,2:end-1,2:end-1));

res = sqrt(norm(resF1, 'fro')^2 + norm(resF2, 'fro')^2 + norm(resF3, 'fro')^2 + norm(resF4, 'fro')^2 + norm(resF5, 'fro')^2 + norm(resF6, 'fro')^2);
end
