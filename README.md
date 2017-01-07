MATLAB CODE WHICH MODELS A WAVE BY SOLVING THE WAVE EQUATION USING FINITE DIFFERENCE SCHEME  
---------------------------------------------------------------------------------------------------------------------------
Nx=100;     % number of grid points in the x-axis
dx=1e-3;    % grid spacing in the x-axis (m)
Ny=100;     % number of grid points in the y-axis
dy=1e-3;    % grid spacing in the y-axis (m)
c=1500;     % speed of sound (m/s)
Nt=100;     % number of time steps
CFL=0.5;    % Courant Friedrichs Lewy number
dt=CFL*dx/c; % set the size of the time step
t=(1:100)*dt*Nt;    % time duration over Nt amount of time steps  


x=(1:Nx)*dx;    % create the grid axis in the x-axis
y=(1:Ny)*dy;    % create the grid axis in the y-axis
x_pos=(Nx/2)*dx;    % set the position of the source in the x-axis
y_pos=(Ny/2)*dy;    %set the position of the source in the y-axis

% set the initial pressure to be a Gaussian Function 
variance=(2*dx)^2;  
gaussian_x=exp(-(x-x_pos).^2/(2*variance));
gaussian_y=exp(-(y-y_pos).^2/(2*variance));
p_n=gaussian_x'*gaussian_y;

% set pressure at (n - 1) to be equal to n (this implicitly sets the
% initial particle velocity to be zero)
p_nm1=p_n;

% preallocate the pressure at (n + 1) (this is updated during the time
% loop)
p_np1=zeros(size(p_n));


% a nested for loop which iterates over time
for n=1:Nt
    
    % two nested for loop which iterates over a 2D space
    for i=2:Nx-1
    for j=2:Ny-1

    % calculates the next value of pressure using 
    % finite difference solution using past and present values of pressure    
    p_np1(i,j)=2*p_n(i,j)-p_nm1(i,j)+...
    ((c*dt/dx)^2)*(p_n(i+1,j)-2*p_n(i,j)+p_n(i-1,j))+...
    ((c*dt/dy)^2)*(p_n(i,j+1)-2*p_n(i,j)+p_n(i,j-1));

    % two end commands underneath causing the end of the 
    % iteration over the 2D grid spacing         
    end
    end
    
    
    
p_nm1=p_n;                  % copy the value at n to (n - 1)
p_n=p_np1;                  % copy the pressure at (n + 1) to n
    

surf(p_np1)                 % plot the pressure field over 3D, where the density of
                            % pressure is indicated by the colour of the field
                            
axis([0 100 0 100 -1 1]);   % set the limit of all three axis in order for it to 
                            % stop auto-ranging
                            
title(['n = ' num2str(n)]); % title which indicates the time frame 
xlabel('x-axis')            % command to label the x-axis
ylabel('y-axis')            % command to label the y-axis
zlabel('Pressure Magnitude')% command to label the z-axis
drawnow;                    % force the plot to update
pause(0.1);                 % briefly pause before updating the loop



% end commands underneath causing the end of the 
% iteration after Nt amount of time steps
end 
