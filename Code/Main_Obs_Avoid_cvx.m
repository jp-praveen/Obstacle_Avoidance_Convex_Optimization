%% Obstacle Avoidance

function [mse_x, mse_y, rmse, eflag_new, err_pos, end_time, dataTransfer] = Main_Obs_Avoid_cvx(probParam, probFlag, probFit)

%% Useful plot settings
% set(0, 'DefaultAxesFontSize', 26, 'DefaultAxesFontWeight','bold')
% set(0, 'DefaultTextFontSize', 26, 'DefaultTextFontWeight','bold')
% % set(gcf, 'PaperPositionMode', 'auto');
% set(gca, 'LooseInset', get(gca,'TightInset'));
% % % set(gcf, 'Color', 'none');         % Set the figure background to transparent
% % set(gca, 'Color', 'none');         % Set the axes background to transparent
% % print('-clipboard', '-dmeta');  % Copy as vector-based EMF with transparency
% % set(gcf, 'InvertHardcopy', 'off');  % Prevent MATLAB from inverting colors
% % % set(gcf, 'Renderer', 'opengl');  % Use OpenGL renderer for better transparency support

%% Problem Parameters
N = probParam.N ;
i_max = probParam.i_max;
alphaProb = probParam.alpha;

H = probParam.H;
w= probParam.w;
r_tol = probParam.r_tol;

% Final time
tf= probParam.tf;

tVec = linspace(0,tf,N);

% Obstacle-1 information
R1= probParam.R1;
H1= probParam.H1;
obsCenter1 = probParam.obsCenter1;

% Obstacle-2 information
R2= probParam.R2;
H2= probParam.H2;
obsCenter2 = probParam.obsCenter2;

% max tilt angle
theta_max= probParam.theta_max;

% control/acceleration
g = [0 0 -9.81];
ai = -g;
af = -g;

% position boundary condition
ri= probParam.ri;
rf= probParam.rf;

% velocity boundary condition
vi= probParam.vi;
vf= probParam.vf;

% acceleration/control min max
umin= probParam.umin;
umax = probParam.umax;

% Initiating error for successive convexification
err_pos = 10;
sol_old = [];
i = 1;

%% Use if Initial Guess is needed
% % Initial guess for trajectory (straight line from ri to rf)
% r_guess = [linspace(ri(1), rf(1), N)', linspace(ri(2), rf(2), N)', linspace(ri(3), rf(3), N)'];
% 
% % Initial guess for velocity (zero velocity)
% v_guess = repmat(vi, N, 1);
% 
% % Initial guess for control input (constant acceleration)
% u_guess = repmat((ai + af) / 2, N, 1);
% 
% % Initial guess for slack variables
% s_guess = repmat(umax, N, 1);  % set to the upper control bound
% s1_guess = 0;  % Initial guess for s1 
% nu1_guess = 0.001;%repmat(0.001, N, 1);
% nu2_guess = 0.001;%repmat(0.001, N, 1);
% 
% initialGuessStruct.r_guess = r_guess;  % Initial guess for r
% initialGuessStruct.v_guess = v_guess;  % Initial guess for v
% initialGuessStruct.u_guess = u_guess;  % Initial guess for u
% initialGuessStruct.s_guess = s_guess;  % Initial guess for s
% initialGuessStruct.s1_guess = s1_guess;  % Initial guess for s1

% Plotting settings
colors = lines(10); 
legendEntries = {};

% Plot Flags
plotIterFlag = probFlag.plotIterFlag;
plotFlag = probFlag.plotFlag;
plotFlagSpline = probFlag.plotFlagSpline;
animationFlag = probFlag.animationFlag;
animation_time = probFlag.animation_time;
videoWriteFlag = probFlag.videoWriteFlag;

% B-Spline fit
numPieces = probFit.numPieces;
order = probFit.order;

for j=1:1:10
tic
    while true
    
        % Solve the problem without obstacle
        [sol_old,fval,eflag] = obstacle_relaxed_cvx_foh(N, tf, g, ri, rf, vi, vf, ai, af, umin, umax, theta_max);
    
        if strcmp(eflag, 'Solved')
            
            rPrev = sol_old.r;
            if plotFlag == 1
                plotIterativeTrajectory(sol_old.r, sol_old.u, 1, 0,colors(1, :),tVec)
                legendEntries{end+1} = ['Relaxed Cons. '];
            end

    
            while i < i_max && err_pos > r_tol
    
                % Successive convexification for avoidaing obstacles
                [sol_new,fval_new,eflag_new] = obstacle_cvx_foh(N, tf, w, obsCenter1,obsCenter2, R1, R2, sol_old, g, ri, rf, vi, vf, ai, af, umin, umax, theta_max);
                err_pos  = max(vecnorm((sol_new.r - rPrev)', "inf"));
                
                rPrev = sol_new.r;
                sol_old = sol_new;
                if plotFlag == 1
                    plotIterativeTrajectory(sol_old.r, sol_old.u, err_pos, i,colors(i+1, :), tVec)
                    legendEntries{end+1} = ['Iteration ', num2str(i), ' Max pos error ', num2str(err_pos)];
                end
                i = i+1;
            end
    
        else
            tf = tf*alpha;
            i = 1;
        end
    
        if err_pos < r_tol
            break;
        end
        i = 1;
    end
if err_pos < r_tol
    end_time=toc;
    break;
end
end

if i < 4
    noObsTraj = 1;
else
    noObsTraj=0;
end

% Parametric representation of the optimized trajectories
[mse_x, mse_y, rmse, coeffs_x, coeffs_y]  = splineFit(numPieces, tVec, sol_new.r(:,1), sol_new.r(:,2), plotFlagSpline);

% Solution data structure
vx_convex = sol_new.v(:,1);
vy_convex = sol_new.v(:,2);

v_convex = sol_new.v;
u_convex = sol_new.u;
r_convex = sol_new.r;
x_convex = r_convex(:,1);
y_convex = r_convex(:,2);

a_convex = u_convex + g;
a_x = a_convex(:,1);
a_y = a_convex(:,2);

% Calculate distances from each trajectory point to the obstacle1 center
distance_to_obstacle1 = sqrt((x_convex - obsCenter1(1)).^2 + (y_convex - obsCenter1(2)).^2);
% Find the index of the minimum distance
[minDistance1, minIndex1] = min(distance_to_obstacle1);

% Extract the waypoint coordinates
xWpClosest1 = x_convex(minIndex1);
yWpClosest1 = y_convex(minIndex1);
vxWpClosest1 = vx_convex(minIndex1);
vyWpClosest1 = vy_convex(minIndex1);
tWpClosest1 = tVec(minIndex1);

% Calculate distances from each trajectory point to the obstacle1 center
distance_to_obstacle2 = sqrt((x_convex - obsCenter2(1)).^2 + (y_convex - obsCenter2(2)).^2);
% Find the index of the minimum distance
[minDistance2, minIndex2] = min(distance_to_obstacle2);

% Extract the waypoint coordinates
xWpClosest2 = x_convex(minIndex2);
yWpClosest2 = y_convex(minIndex2);
vxWpClosest2 = vx_convex(minIndex2);
vyWpClosest2 = vy_convex(minIndex2);
tWpClosest2 = tVec(minIndex2);

XIniObs1Traj = x_convex(1:minIndex1);
YIniObs1Traj = y_convex(1:minIndex1);
tVec1 = tVec(1:minIndex1);
[~, ~, ~, XIniObs1Trajcoeffs_x, XIniObs1Trajcoeffs_y]  = splineFit(numPieces, tVec1, XIniObs1Traj, YIniObs1Traj, plotFlagSpline);
dataTransfer.XIniObs1Trajcoeffs_x = XIniObs1Trajcoeffs_x;
dataTransfer.XIniObs1Trajcoeffs_y = XIniObs1Trajcoeffs_y;

XObs1Obs2Traj = x_convex(minIndex1:minIndex2);
YObs1Obs2Traj = y_convex(minIndex1:minIndex2);
if length(XObs1Obs2Traj) < 8
    noSuffPt = 1;
else
    noSuffPt=0;
end

if noSuffPt == 0
    tVec2 = tVec(minIndex1:minIndex2);%linspace(tVec(minIndex1),tVec(minIndex2),15);
    [~, ~, ~, XObs1Obs2Trajcoeffs_x, XObs1Obs2Trajcoeffs_y]  = splineFit(numPieces, tVec2, XObs1Obs2Traj, YObs1Obs2Traj, plotFlagSpline);
    dataTransfer.XObs1Obs2Trajcoeffs_x = XObs1Obs2Trajcoeffs_x;
    dataTransfer.XObs1Obs2Trajcoeffs_y = XObs1Obs2Trajcoeffs_y;
end

XObs2EndTraj = x_convex(minIndex2:end);
YObs2EndTraj = y_convex(minIndex2:end);
tVec3 = tVec(minIndex2:end);%linspace(tVec(minIndex2),tf,15);
[~, ~, ~, XObs2EndTrajcoeffs_x, XObs2EndTrajcoeffs_y]  = splineFit(numPieces, tVec3, XObs2EndTraj, YObs2EndTraj, plotFlagSpline);
dataTransfer.XObs2EndTrajcoeffs_x = XObs2EndTrajcoeffs_x;
dataTransfer.XObs2EndTrajcoeffs_y = XObs2EndTrajcoeffs_y;

dataTransfer.xWpClosest1 = xWpClosest1;
dataTransfer.yWpClosest1 = yWpClosest1;
dataTransfer.vxWpClosest1 = vxWpClosest1;
dataTransfer.vyWpClosest1 = vyWpClosest1;
dataTransfer.tWpClosest1 = tWpClosest1;

dataTransfer.xWpClosest2 = xWpClosest2;
dataTransfer.yWpClosest2 = yWpClosest2;
dataTransfer.vxWpClosest2 = vxWpClosest2;
dataTransfer.vyWpClosest2 = vyWpClosest2;
dataTransfer.tWpClosest2 = tWpClosest2;

dataTransfer.coeffs_x = coeffs_x;
dataTransfer.coeffs_y = coeffs_y;

dataTransfer.x_convex = x_convex;
dataTransfer.y_convex = y_convex;

dataTransfer.noSuffPt = noSuffPt;
dataTransfer.noObsTraj = noObsTraj;

dataTransfer.uConvex = u_convex;

if plotFlag == 1
    set(0, 'DefaultAxesFontSize', 26, 'DefaultAxesFontWeight','bold')
    set(0, 'DefaultTextFontSize', 26, 'DefaultTextFontWeight','bold')
    figure(1)
    box on
    axis tight
    axis padded
    % axis equal
    hold on
    draw_circle(obsCenter1(1,1),obsCenter1(1,2),R1,'k');
    draw_circle(obsCenter2(1,1),obsCenter2(1,2),R2,'k');
    lgd1 = legend(legendEntries, 'Location', 'Best');
    set(lgd1, 'FontSize', 16, 'FontWeight', 'bold');  % Adjust the legend font size and weight
    hold off;
    
    figure(2)
    box on
    axis tight
    axis padded
    hold on
    lgd2 = legend(legendEntries, 'Location', 'Best');
    set(lgd2, 'FontSize', 16, 'FontWeight', 'bold');  % Adjust the legend font size and weight
    hold off;
end

%% Plots
% Position
if plotFlag == 1
    figure(3)
    plot3(sol_new.r(:,1), sol_new.r(:,2), sol_new.r(:,3), '--', 'LineWidth',3,'MarkerSize',3, 'Color', 'b', 'DisplayName', 'Convex Trajectory' )
    hold on
    % plot3(sol_new.r(1,1), sol_new.r(1,2), sol_new.r(1,3),'.','MarkerEdgeColor','r','MarkerFaceColor',[.49 1 .63],'MarkerSize',2)
    % hold on
    % plot3(sol_new.r(end,1), sol_new.r(end,2), sol_new.r(end,3),'.','MarkerEdgeColor','r','MarkerFaceColor',[.49 1 .63],'MarkerSize',2)
        xlabel('X [m]')
    ylabel('Y [m]')
    zlabel('Z [m]')
    hold on
    grid on
    % quiver3(sol_new.r(1:end,1), sol_new.r(1:end,2), sol_new.r(1:end,3),sol_new.u(:,1), sol_new.u(:,2), sol_new.u(:,3),0.5, 'Color','k');
    draw_circle(obsCenter1(1,1),obsCenter1(1,2),R1,'g');
    draw_circle(obsCenter2(1,1),obsCenter2(1,2),R2,'k');
    axis equal
    plot(xWpClosest1, yWpClosest1,'.','MarkerEdgeColor',[0 1 1],'MarkerFaceColor',[0 1 1],'MarkerSize',20)
    plot(xWpClosest2, yWpClosest2,'.','MarkerEdgeColor',[.7 .75 .93],'MarkerFaceColor',[.7 .75 .93],'MarkerSize',20)
    axis padded
end

if animationFlag == 1
    x = sol_new.r(:, 1);  % X positions from the trajectory
    y = sol_new.r(:, 2);  % Y positions from the trajectory
    z = sol_new.r(:, 3);  % Z positions from the trajectory
    
    % Roll, pitch, yaw are assumed to be zeros or can be derived
    roll = zeros(size(x));  % Assuming no roll for simplicity
    pitch = zeros(size(x));  % Assuming no pitch for simplicity
    yaw = zeros(size(x));    % Assuming no yaw for simplicity
    
    obsCenters = [obsCenter1; obsCenter2];  % 2 obstacles at different positions
    obsRadii = [R1, R2];  % Radii of the obstacle
    
    % Call the animation function
    drone_Animation(x, y, z, roll, pitch, yaw, obsCenters, obsRadii,animation_time, videoWriteFlag);
    
end

%% Plot Functions
function plotIterativeTrajectory(rPlot, uPlot, tolerance, iteration, color, tVec)
    figure(1)
    % Hold the figure to plot multiple iterations on the same figure
    hold on;

    % Plot the trajectory
    plot3(rPlot(:,1), rPlot(:,2), rPlot(:,3), '-s','Color', color, 'LineWidth', 3);
    
    % Plot quiver for control (control vectors)
    % quiver3(r(:,1), r(:,2), r(:,3), u(:,1), u(:,2), u(:,3), 0.5, 'Color', color);
    
    % Set labels and grid
    xlabel('X [m]');
    ylabel('Y [m]');
    zlabel('Z [m]');
    % title('3D Trajectories with Control Quivers');
    grid on;
    hold off;  % Release the figure hold for the next update

    % Acceleration
    figure(2)
    hold on
    plot(tVec,vecnorm(uPlot'), '-s', 'LineWidth',3, 'MarkerSize',6, 'Color', color)
    xlabel('Time [s]')
    ylabel('||u|| [m/s^2]')
    grid on
    hold off
    end

function h = draw_circle(x, y, r, color)
    % Inputs:
    %   x - x-coordinate of the circle center
    %   y - y-coordinate of the circle center
    %   r - radius of the circle
    %   color - color of the circle, default is 'b' (blue)
    %
    % Output:
    %   h - handle to the plot
    %
    % Example usage:
    %   draw_circle(0, 0, 5, 'r');  % Draw a red circle at (0,0) with radius 5

    if nargin < 4
        color = 'b';  % Default color is blue
    end

    theta = linspace(0, 2*pi, 100);  % Create 100 points for smoothness
    x_circle = x + r * cos(theta);   % X points of the circle
    y_circle = y + r * sin(theta);   % Y points of the circle

    h = plot(x_circle, y_circle, 'Color', color, 'LineWidth', 3);  % Draw the circle
    axis equal;  % Make sure the aspect ratio is equal
    hold on;     % Keep the plot
end

function animation = drone_Animation(x, y, z, roll, pitch, yaw, obsCenters, obsRadii,animation_time, videoFlag)
    % Animate the quadcopter movement along a 3D trajectory
    % with plotted obstacles.
    % 
    % Inputs:
    %   x, y, z - Quadcopter trajectory positions
    %   roll, pitch, yaw - Quadcopter orientation angles
    %   obsCenters - Obstacles' center positions (Nx3 matrix for N obstacles)
    %   obsRadii - Radii of obstacles (1xN vector)

    %% Video recording setup
    if videoFlag
        videoFileName = 'drone_animation.mp4';  % File name for the video
        videoObj = VideoWriter(videoFileName, 'MPEG-4');  % Create a VideoWriter object
        videoObj.FrameRate = 2;  % the frame rate
        videoObj.Quality = 100;  % maximum quality (0 to 100)
        open(videoObj);  % Open the video file for writing
    end

    %% Define design parameters
    D2R = pi/180;
    b   = 0.3;   % the length of total square cover by whole body of quadcopter in meter
    a   = b/3;   % the legth of small square base of quadcopter(b/4)
    H   = 0.06;  % hight of drone in Z direction (4cm)
    H_m = H+H/2; % hight of motor in z direction (5 cm)
    r_p = b/4;   % radius of propeller
    %% Conversions
    ro = 45*D2R;                   % angle by which rotate the base of quadcopter
    Ri = [cos(ro) -sin(ro) 0;
          sin(ro) cos(ro)  0;
           0       0       1];     % rotation matrix to rotate the coordinates of base 
    base_co = [-a/2  a/2 a/2 -a/2; % Coordinates of Base 
               -a/2 -a/2 a/2 a/2;
                 0    0   0   0];
    base = Ri*base_co;             % rotate base Coordinates by 45 degree 
    to = linspace(0, 2*pi);
    xp = r_p*cos(to);
    yp = r_p*sin(to);
    zp = zeros(1,length(to));

    %% Define Figure plot
    fig1 = figure('pos', [0 50 800 600]);
    % fig1 = figure('Position', [100, 100, 1920, 1080]);  % Set figure size to 1920x1080
    hg   = gca;
    % view(68,53);
    view(40, 45)
    % view(0, 90);  %  top-down view
    box on
    grid on;
    % axis equal;
    % axis tight
    xlim([-1.5 1.5] + [min(x) max(x)]); % Dynamic x-limits based on trajectory
    ylim([-1.5 1.5] + [min(y) max(y)]); % Dynamic y-limits based on trajectory
    zlim([0 3.5] + [min(z) max(z)]);    % Dynamic z-limits based on trajectory
    % title('Quadcopter Drone Animation with Trajectory and Obstacles')
    xlabel('X [m]');
    ylabel('Y [m]');
    zlabel('Z [m]');
    hold(gca, 'on');

    %% Plot obstacles (assumes spherical obstacles)
    numObstacles = size(obsCenters, 1);  % Number of obstacles
    for i = 1:numObstacles
        % [X, Y, Z] = sphere(20);  % Create sphere for obstacle
        % surf(obsCenters(i, 1) + obsRadii(i)*X, ...
        %      obsCenters(i, 2) + obsRadii(i)*Y, ...
        %      obsCenters(i, 3) + obsRadii(i)*Z, ...
        %      'FaceColor', 'b');  % Plot obstacle
        draw_tree_model_stl(obsCenters(i, 1), obsCenters(i, 2), 0.009, 'tree_model.stl')
    end

    %% Initialize the trajectory line to plot incrementally
    h_traj = plot3(nan, nan, nan, 'r.-', 'LineWidth', 2.5);  % Initialize trajectory line
    
    %% Design Different parts of the quadcopter
    % design the base square
    drone(1) = patch([base(1,:)],[base(2,:)],[base(3,:)],'r');
    drone(2) = patch([base(1,:)],[base(2,:)],[base(3,:)+H],'r');
    alpha(drone(1:2),0.7);
    % design 2 perpendicular legs of quadcopter 
    [xcylinder ycylinder zcylinder] = cylinder([H/2 H/2]);
    drone(3) = surface(b*zcylinder-b/2,ycylinder,xcylinder+H/2,'facecolor','b');
    drone(4) = surface(ycylinder,b*zcylinder-b/2,xcylinder+H/2,'facecolor','b'); 
    alpha(drone(3:4),0.6);
    % design 4 cylindrical motors 
    drone(5) = surface(xcylinder+b/2,ycylinder,H_m*zcylinder+H/2,'facecolor','r');
    drone(6) = surface(xcylinder-b/2,ycylinder,H_m*zcylinder+H/2,'facecolor','r');
    drone(7) = surface(xcylinder,ycylinder+b/2,H_m*zcylinder+H/2,'facecolor','r');
    drone(8) = surface(xcylinder,ycylinder-b/2,H_m*zcylinder+H/2,'facecolor','r');
    alpha(drone(5:8),0.7);
    % design 4 propellers
    drone(9)  = patch(xp+b/2,yp,zp+(H_m+H/2),'c','LineWidth',0.5);
    drone(10) = patch(xp-b/2,yp,zp+(H_m+H/2),'c','LineWidth',0.5);
    drone(11) = patch(xp,yp+b/2,zp+(H_m+H/2),'p','LineWidth',0.5);
    drone(12) = patch(xp,yp-b/2,zp+(H_m+H/2),'p','LineWidth',0.5);
    alpha(drone(9:12),0.3);

    %% Create a group object and parent surface for transformation
    combinedobject = hgtransform('parent',hg );
    set(drone,'parent',combinedobject)

    %% Animation loop
    for i = 1:length(x)
        % Update the quadcopter's transformation: translation and rotation
        translation = makehgtform('translate', [x(i) y(i) z(i)]);
        rotation1 = makehgtform('xrotate', (pi/180)*(roll(i)));
        rotation2 = makehgtform('yrotate', (pi/180)*(pitch(i)));
        rotation3 = makehgtform('zrotate', (pi/180)*(yaw(i)));

        % Apply the transformations
        set(combinedobject,'matrix', translation * rotation3 * rotation2 * rotation1);

        % Incrementally update the trajectory
        set(h_traj, 'XData', x(1:i), 'YData', y(1:i), 'ZData', z(1:i));

        if videoFlag
            % Capture frame using exportgraphics and read it back
            exportgraphics(gcf, 'current_frame.png', 'Resolution', 300);  % High-quality capture
            frame = imread('current_frame.png');  % Read the exported image
            frame = imresize(frame, [2294, 2294]);
            % Ensure the dimensions are multiples of 2 (required by the video codec)
            [height, width, ~] = size(frame);
            newHeight = floor(height / 2) * 2;
            newWidth = floor(width / 2) * 2;
        
            % Resize the frame only if necessary
            if height ~= newHeight || width ~= newWidth
                frame = imresize(frame, [newHeight, newWidth]);  % Resize to ensure multiples of 2
            end
        
            % Write the frame to the video file
            writeVideo(videoObj, frame);
        end

        % Pause for animation effect
        drawnow;
        pause(animation_time);
    end
    %% Close the video file
    if videoFlag
        close(videoObj);  % Close the video file after recording
    end
    function h_fig = draw_tree_model_stl(x, y, scale, model_file)
    % Import and display a 3D tree model from an STL file.
    %
    % Inputs:
    %   x, y       - position to place the tree
    %   scale      - scaling factor for the model
    %   model_file - path to the 3D model file (STL format)
    %
    % Output:
    %   h_fig - handle to the figure
    %
    % Example usage:
    %   draw_tree_model(0, 0, 1, 'tree_model.stl');

    % Read the 3D model file
    fv = stlread(model_file);

    % Extract the faces and vertices
    faces = fv.ConnectivityList;
    vertices = fv.Points;

    % Scale and position the model
    vertices = vertices * scale;
    vertices(:,1) = vertices(:,1) + x;
    vertices(:,2) = vertices(:,2) + y;

    % Plot the model
    % figure;
    patch('Faces', faces, 'Vertices', vertices, ...
          'FaceColor',       [0.5 0.8 0.5], ...
          'EdgeColor',       'none',        ...
          'FaceLighting',    'gouraud',     ...
          'AmbientStrength', 0.15);
    camlight('headlight');
    material('dull');
    grid on;

    h_fig = gcf;
end
end
end 
