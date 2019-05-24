clear all; clc; clf; 
Nx = 101; Ny = 101;
x = linspace(0,10,Nx); y = linspace(0,10,Ny);
h = x(2)-x(1); dt = 0.2*h^2; T = 1; Nt = T/dt;
Du = 1; k2 = 11; 
Dv = 0.04; k1 = 7; 
%Dv = 0.02; k1 = 9; 
u(1:Nx,1:Ny) = 1 + 0.04*k2 + 0.1*rand(Nx,Ny);
v(1:Nx,1:Ny) = 0.2*k2 + 0.1*rand(Nx,Ny);
nu = u; nv = v; 
for n = 1:Nt-1
    for i = 2:Nx-1
        for j = 2:Ny-1
            nu(i,j) = u(i,j) + dt*(...
                  Du*(u(i+1,j) + u(i-1,j) - 4*u(i,j) + u(i,j+1) + u(i,j-1))/h^2 ...
                  + k1*(v(i,j) - u(i,j)*v(i,j)/(1 + (v(i,j))^2))...
                      );
            nv(i,j) = v(i,j) + dt*(...
                 Dv*(v(i+1,j) + v(i-1,j) - 4*v(i,j) + v(i,j+1) + v(i,j-1))/h^2 ...
                  + k2  - v(i,j) - 4*u(i,j)*v(i,j)/(1 + (v(i,j))^2));
        end
    end
    nu(1,:) = nu(2,:);  nu(end,:) = nu(end-1,:); 
    nu(:,1) = nu(:,2);     nu(:,end) = nu(:,end-1);  
    nv(1,:) = nv(2,:);  nv(end,:) = nv(end-1,:); 
    nv(:,1) = nv(:,2);     nv(:,end) = nv(:,end-1);  
    u = nu; v = nv;
     
    if mod(n,20) == 0
    figure(1); 
    subplot(2,1,1); mesh(x,y,nu'); title('u'); view(90,90);axis image
    subplot(2,1,2); mesh(x,y,nv'); title('v'); view(90,90);axis image
    end
end