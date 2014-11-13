ball_start_state = [1; 3; 0; 0];

x_hist = [];
t_hist = [];
measurement_hist = [];
error_hist = [];

x = ball_start_state;
dt = 0.01;
radius = 1;
damping = 0.8;
gravity = [0; -9.8];
figure
for t=0:dt:1
  x_hist = [x_hist x];
  measurement_hist = [measurement_hist; x(1:2)];
  t_hist = [t_hist t];
  
  
  %% SLAM
  % Big optimization
  A = eye(size(x_hist, 2)*4+size(x_hist, 2)*2, size(x_hist, 2)*4);
  b = zeros(size(A, 1), 1);
  for i=1:size(x_hist, 2)
    if i==1
      % First state is set as prior
      b(1:4) = ball_start_state;
    else
      starti = 1+4*(i-1); endi = starti+3;
      % x(i, 1) = x(i-1) + dt*x(i-1, 3)
      % sim for x2
      % sim for x3, x4 except it's dependence on gravity
      A(starti:endi, starti-4:endi-4) = -[1 0 dt 0;
                                     0 1 0 dt;
                                     0 0 1 0;
                                     0 0 0 1];
      b(starti:endi) = [0;0;0;gravity(2)*dt];
    end
    % And measurements
    last_section_start = size(x_hist, 2)*4+2*i-1;
    A(last_section_start:last_section_start+1, 1+(i-1)*4:1+(i-1)*4+1) = eye(2);
  end
  % And then measurements
  b(1+size(x_hist, 2)*4:end) = measurement_hist;
  
  % Optimal
  % A.'A*θ = A.'*b
  % so
  optimal_state = (A.'*A)^(-1)*A.'*b;
  
  pos_guess_this_cycle = optimal_state(end-3:end);
  error = x - pos_guess_this_cycle
  error_hist = [error_hist error];
  
  
  %% DYNAMICS AND DRAWING
  % Reverse if we clip the ground
  if (x(2) < (dt*abs(x(4))) + radius)
    x(4) = -x(4)*damping;
  end
  % Update pos
  x(1:2) = x(1:2) + dt*x(3:4);
  % Update vel
  x(3:4) = x(3:4) + dt*gravity;
  
  % DRAW
  clf;
  axisAnnotation('ellipse',...               % draw circle
    [x(1); x(2); 0; 0] + radius*[-1;-1;2;2],...  % [x y w h]
    'FaceColor','r');         gurobi               % make it red
  line([-5,5]*radius,[0,0],'Color','k','LineWidth',1.5);
  axis equal;
  axis(radius*[-5 5 -.5 9.5]);
  drawnow
  t
end

snopt
