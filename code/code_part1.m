% ELEC 4700 Assignment 3 Part 1
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

% ~~~~~~~~~~ PART 1 ~~~~~~~~~~

mass = C.m_0*0.26; % effective e mass
k = C.kb;
T = 300;
vth = sqrt(2*((k*T)/(mass))); % thermal velocity
step = 1e-15;
ts = 1000;
dt = 1e-15;
t = 0;
step_size = 1e-9;
dv = step_size/vth;

num_particles = 100;

boundaryX = 2e-7;
boundaryY = 1e-7;

MTC = 0.2e-12;

current = zeros(1,ts);
avgV = zeros(1,ts);
mu = zeros(1,ts);
temp = zeros(1,ts);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1(a) - If a voltage of 0.1V is applied across the x dimension of the semiconductor, what
% is the electric field on the electrons? You can assume it is constant over the
% semiconductor.

v = 0.5; % 0.5 V instead of 0.1 V to get a better trajectory plot

eField = v/boundaryX;
fprintf('The field is %d. ', eField);

% 1(b) - What is the force on each electron?

force = eField*(C.q_0);
fprintf('The force is %d. ', force);

% 1(c) - Calculate the acceleration on the electrons and use this in your model to update
% the velocity of each electron at each time step. Add the capability for the electrons
% to respond to a static electric field with both an x and a y component. Plot the
% trajectories of the electrons. They should be curved! Increase the electric field
% and see the curve! Be careful that your time step is appropriate!

accel = force/(mass);
fprintf('The acceleration is %d. ', accel);

% random initial velocities
Vx = -sqrt(2*C.kb*T/mass + log(num_particles) + (1/2)*log(2*pi*C.kb*T/mass))...
    .*randn(1,num_particles);
Vy = -sqrt(2*C.kb*T/mass + log(num_particles) + (1/2)*log(2*pi*C.kb*T/mass))...
    .*randn(1,num_particles);   

% random initial positions
for i = 1:num_particles
    x(i) = rand()*boundaryX;
    y(i) = rand()*boundaryY;
end

% update velocity using E field acceleration
new_dv = accel*(dv);

figure(1);
clf

for t=1:ts
    
    %t = t + dt;
    
    % update position
    x(1:num_particles) = x(1:num_particles) + (dt .* Vx(1:num_particles));
    y(1:num_particles) = y(1:num_particles) + (dt .* Vy(1:num_particles));
    
    
    for i=1:num_particles
        
       Vx(i) =  Vx(i) + new_dv;
        
       % reflection for y
       if y(i) >= boundaryY || y(i) <= 0
           Vy(i) = - Vy(i);
       end
       
       % boundary condition for x
       if x(i) <= 0
           x(i) = x(i) + boundaryX;
       end
       if x(i) >= boundaryX
           x(i) = x(i) - boundaryX;
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
    
    % Collect velocities for current plot
    v_collect(t) = sqrt(Vx(i)^2 + Vy(i)^2);
    
    % Current
    avgV = sum(v_collect)/num_particles;
    mu = (avgV/eField);
    current(t) = (C.q_0)*10e15*mu*eField/(boundaryX*boundaryY);
    
    % Temperature
    V2 = (Vx.^2 + Vy.^2);
    temp(t) = (mean(V2)*mass)/(2*k);

    % trace some particles for plotting
    particlestotrace = 10;
    for i=1:1:particlestotrace
        colors = hsv(particlestotrace);
        plot(x(i),y(i), 'o', 'markers', 1, 'color', colors(i,:), 'MarkerFaceColor', colors(i,:));
    end
    
    % trajectory plot
    title(['Trajectory (' ,num2str(particlestotrace), ' electrons)']);
    axis([0 boundaryX 0 boundaryY]);
    hold on;
    pause(0.01);
    
end

% PLOTS

% current plot
figure(2)
plot(linspace(1,ts,ts),current)
title('Current Plot')
xlabel('ts')
ylabel('Current/cm2')


% density plot
density = [transpose(x), transpose(y)];
figure(3);
hist3(density, [25,25]);
hold on;
xlabel('X');
ylabel('Y');
zlabel('Electrons');
title('Density Map');

% temperature plot
figure(4);
timevec = linspace(1,ts,ts);
plot(timevec,temp, 'r')
xlabel('Time');
ylabel('Temp (K)');
title('Avg. Temp');








    


