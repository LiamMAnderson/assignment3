% ELEC 4700 Assignment 3 Part 2 & 3
% Liam Anderson 100941879
% Submission March 15 2020

global C

C.q_0 = 1.60217653e-19;             % electron charge
C.hb = 1.054571596e-34;             % Dirac constant
C.h = C.hb * 2 * pi;                % Planck constant
C.m_0 = 9.10938215e-31;             % electron mass
C.kb = 1.3806504e-23;               % Boltzmann constant
C.eps_0 = 8.854187817e-12;          % vacuum permittivity
C.mu_0 = 1.2566370614e-6;           % vacuum permeability

Mass = C.m_0*0.26;
k = C.kb;
T = 300;
vth = sqrt(2*((k*T)/(Mass))); % thermal velocity
dt = 1e-15;
t = 0;
step_size = 1e-9;
MTC = 0.2e-12;
step_size = 1e-9;
dv = step_size/vth;
ts = 1000;

boundaryX = 2e-7;
boundaryY = 1e-7;

% 3/2 ratio
L = 60;
W = 40;

G = sparse(L*W,L*W);
V = zeros(L*W,1);
con = zeros(L,W);

nx = L;
ny = W;

siginside = 1e-2;
sigoutside = 1;


for i = 1:nx
    for j = 1:ny
        if ((i > (nx/2 - nx/6) && i < (L/2 + nx/6)) && (j > (ny/2 + ny/6) || j < (ny/2 - ny/6)))  
            con(i,j) = siginside;
        else
            con(i,j) = sigoutside;
        end
    end
end


for i=1:nx
    for j=1:ny
        % MAPPING
        n = (j + (i-1)*ny);
        nxm = (j + (i-2)*ny);
        nxp = (j + (i)*ny);
        nyp = (j+1 + (i-1)*ny);
        nym = (j-1 + (i-1)*ny);
        G(n,:) = 0;
        % BOUNDARY CONDITIONS
        if (i == 1)
            V(n) = 1;
            G(n,n) = con(i,j);
        elseif (i == L)
            V(n) = 0;
            G(n,n) = con(i,j);
        elseif (j == 1)
            G(n,nxm) = (con(i-1,j) + con(i,j))/2;
            G(n,nxp) = (con(i+1,j) + con(i,j))/2;
            G(n,nyp) = (con(i,j+1) + con(i,j))/2;
            G(n,n) = -(G(n,nxm) + G(n,nxp) + G(n,nyp));
        elseif (j == W)
            G(n,nxm) = (con(i-1,j) + con(i,j))/2;
            G(n,nxp) = (con(i+1,j) + con(i,j))/2;
            G(n,nym) = (con(i,j-1) + con(i,j))/2;
            G(n,n) = -(G(n,nxm) + G(n,nxp) + G(n,nym));
        else
            G(n,nxm) = (con(i-1,j) + con(i,j))/2;
            G(n,nxp) = (con(i+1,j) + con(i,j))/2;
            G(n,nyp) = (con(i,j+1) + con(i,j))/2;
            G(n,nym) = (con(i,j-1) + con(i,j))/2;
            G(n,n) = -(G(n,nxm) + G(n,nxp) + G(n,nyp) + G(n,nym));
        end
        
    end
    
end

surfplot = zeros(L,W);

% Potential
F = G\V;
fprintf('The potential is %d V. ', F);

% e field
eField = mean(F)/boundaryX;

% accel
accel = force/(mass);

% new v
new_dv = accel*dv;



% BOTTLENECK SIMULATION

num_particles = 1000;

% random initial velocities
Vx = -sqrt(2*C.kb*T/mass + log(num_particles) + (1/2)*log(2*pi*C.kb*T/mass))...
    .*randn(1,num_particles);
Vy = -sqrt(2*C.kb*T/mass + log(num_particles) + (1/2)*log(2*pi*C.kb*T/mass))...
    .*randn(1,num_particles);   

% random initial positions
for i = 1:num_particles
    x(i) = rand()*boundaryX;
    y(i) = rand()*boundaryY;
    
    % INITIAL BOTTLENECK CHECK
    % Upper box (x-dimensions 0.8->1.2 and y-dimensions 0.6->1)
    if x(i) >= 0.6e-7 & x(i) <= 1.2e-7 & y(i) >= 0.6e-7
        x(i) = rand(1,1)*2e-7;
        y(i) = rand(1,1)*1e-7;
    end

    % Lower box (x-dimensions 0.8->1.2 and y-dimensions 0->0.4)
    if x(i) >= 0.6e-7 & x(i) <= 1.2e-7 & y(i) <= 0.4e-7
        x(i) = rand(1,1)*2e-7;
        y(i) = rand(1,1)*1e-7;
    end
    
    
end

figure(1);
clf

for t = 1:ts
    
    x(1:num_particles) = x(1:num_particles) + (dt .* Vx(1:num_particles));
    y(1:num_particles) = y(1:num_particles) + (dt .* Vy(1:num_particles));
    
    %t = t + dt;
    
    for i=1:1:num_particles
        
       Vx(i) =  Vx(i) + new_dv;
        
       % reflection for y
       if y(i) >= boundaryY || y(i) <= 0
           Vy(i) = - Vy(i);
       end
       % boundary condition for x
       if x(i) < 0
           x(i) = boundaryX;
       end
       if x(i) > boundaryX
           x(i) = 0;
       end
       
       % BOTTLENECK CHECK
       % Upper box
       if x(i) >= 0.8e-7 && x(i) <= 1.2e-7 && y(i) >= 0.6e-7
           if y(i) > 0.6e-7 % collision with side of upper box
               Vx(i) = -Vx(i);
           end
           if y(i) >= 0.6e-7 % collision with bottom of upper box
               Vy(i) = -Vy(i);
           end
       end

       % Lower box
       if x(i) >= 0.8e-7 && x(i) <= 1.2e-7 && y(i) <= 0.4e-7
           if y(i) < 0.4e-7 % collision with side of lower box
               Vx(i) = -Vx(i);
           end
           if y(i) <= 0.4e-7 % collision with top of lower box
               Vy(i) = -Vy(i);
           end
       end    
       
       
       % SCATTERING
       dv = step_size/vth;
       p_scat(i) = 1 - exp(-dv/MTC);
       random(i) = rand();
       if p_scat(i) > random(i); % if scatter assign new random velocity
           Vx(i) = rand()*vth*10*cos(2*pi*randn());
           Vy(i) = rand()*vth*10*sin(2*pi*randn());  
       end
       
    end
    
    
    % trace some particles for plotting
    particlestotrace = 10;
    for i=1:1:particlestotrace
        colors = hsv(particlestotrace);
        plot(x(i),y(i), 'o', 'markers', 1, 'color', colors(i,:), 'MarkerFaceColor', colors(i,:));
    end
   
    title(['Trajectory   (' ,num2str(particlestotrace), ' electrons)']);
    axis([0 boundaryX 0 boundaryY]);
    rectangle('Position', [0.8e-7 0 0.4e-7 0.4e-7]);
    rectangle('Position', [0.8e-7 0.6e-7 0.4e-7 0.4e-7]);
    hold on;
    pause(0.01);
end


% PLOTTING
for (i=1:L)
    for (j=1:W)
        n = j + (i-1)*W;
        surfplot(i,j) = F(n);
    end
end

figure(10)
surf(surfplot)
title('Potential (Surface) Plot')
xlabel('x')
ylabel('y')
zlabel('z')

[X, Y] = meshgrid(1:W,1:L);
[Ex, Ey] = gradient(surfplot);

figure(11)
quiver(X, Y)
title('EF (Quiver)')
xlabel('x')
ylabel('y')

density = [transpose(x), transpose(y)];
figure(12);
hist3(density, [25,25]);
hold on;
xlabel('X');
ylabel('Y');
zlabel('Electrons');
title('Density Map');

