function [x, t, u] = calctemp(tmax, nt, xmax, nx, method, timeData, tempData)
% Function for modelling temperature in a space shuttle tile
% D N Johnston  14/02/24
%
% Input arguments:
% tmax   - maximum time (s)
% nt     - number of timesteps
% xmax   - total thickness (m)
% nx     - number of spatial steps
% method - solution method ('forward', 'backward' etc)
% timeData - time vector for surface temperatures (s)
% tempData - surface temperature vector (C or K)
%
% Return arguments:
% x      - distance vector (m)
% t      - time vector (s)
% u      - temperature matrix (C or K)
%
% For example, to perform a  simulation with 501 time steps
%   [x, t, u] = calctemp(4000, 501, 0.05, 21, 'forward', timeData, tempData);
%

% Set material properties and derived values (LI-900)
% Obtained from NASA document: Structures and Materials: Space Shuttle Tiles, Grades 5-8 - NASA
% Note that we're assuming constant properties.
thermCon = 0.0484; % W/m K; 0.028 BTU/ft/hr/F at 70F and 1 atm
density  = 144;    % kg/m^3; 9 lb/ft^3
specHeat = 628;    % J/kg/K; 0.15 Btu/lb/F
thermDiff = thermCon/(density*specHeat);



% Initialise everything.
dt = tmax / (nt-1);
t = (0:nt-1) * dt;
dx = xmax / (nx-1);
x = (0:nx-1) * dx;
u = zeros(nt, nx);
p = (thermDiff*dt)/(dx)^2;
% Use interpolation to get outside temperature at time vector t 
% and store it as left-hand boundary vector L.

timeData(1, length(timeData) + 1) = 4000;
tempData(1, length(tempData)+1) = 0;

L = interp1(timeData, tempData, t, "linear", "extrap");

% set initial conditions equal to boundary temperature at t=0.
u(1, :) = L(1);
R = 0;


% Select method and run simulation.
switch method
    case 'forward'
        u(:, 1) = L; % Outside boundary condition
        u(:, nx) = 0; % Inside boundary condition; set to zero as a starting point.
        % You need to put your solution code here.
        i = 2:nx; % set up index vectors
        im = 1:nx-1;
        ip = [3:nx nx-1];

               
        for n=1:nt-1
            u(n+1, i) = (1 - 2 * p) * u(n, i) + p * (u(n, im) + u(n, ip));
        end
     
    case 'dufort-frankel'
        u(:, 1) = L;
        i = 2:nx; % set up index vectors
        im = 1:nx-1;
        ip = [3:nx nx-1];

        u(:, nx) = 0;

        for n=1:nt-1
            if n == 1
                nm = 1;
            else
                nm = n-1;
            end
            u(n+1,i) = ((1 - 2*p)*u(nm,i) + 2*p*(u(n,im) + u(n,ip))) / (1 + 2*p);
        end
    case 'backward'
        
        for n=1:nt-1
            
            u(:, 1) = L; 

            indexVec = 2:nx-1;
            
            B(1) = 1;
            C(1) = 0;
            D(1) = L(n+1);
           
            A(indexVec) = -p;
            B(indexVec) = 1 + 2*p;
            C(indexVec) = -p;
            D(indexVec) = u(n, indexVec);
           
            A(nx) = -2*p;
            B(nx) = 1 + 2*p;
            D(nx) = u(n,nx);

            u(n+1,:) = tdm(A,B,C,D);
            
        end
    case 'Crank-Nicolson'
        for n=1:nt-1
            
            u(:, 1) = L; 

            indexVec = 2:nx-1;

            B(1) = 1;
            C(1) = 0;
            D(1) = L(n+1);

            A(indexVec) = -p / 2;
            B(indexVec) = 1 + p;
            C(indexVec) = -p/2;
            D(indexVec) = (p/2)*u(n, indexVec) + (1-p)*u(n, indexVec) + (p/2)*u(n, indexVec);

            A(nx) = -p;
            B(nx) = 1 + p;
            D(nx) =p*u(n,nx-1) + (1 - p)*u(n,nx);


            u(n+1,:) = tdm(A,B,C,D);
        end
    otherwise
        error (['Undefined method: ' method])
end

figure (2);
surf(x, t, u);

% this gives a better appearance
shading interp

% rotate plot. 
view(140,30)

xlabel('\itx\rm - m');
ylabel('\itt\rm - s');
zlabel('\itu\rm - deg C');
% title(['DuFort-Frankel, \itp\rm = ' num2str(p)])


% Tri-diagonal matrix
function x = tdm(a,b,c,d)
n = length(b);


for i = 2:n
    factor = a(i) / b(i-1);
    b(i) = b(i) - factor * c(i-1);
    d(i) = d(i) - factor * d(i-1);
end

x(n) = d(n) / b(n);


for i = n-1:-1:1
    x(i) = (d(i) - c(i) * x(i+1)) / b(i);
end

    