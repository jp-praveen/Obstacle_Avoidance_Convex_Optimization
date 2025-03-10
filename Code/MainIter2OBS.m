%% Author: Praveen Jawaharlal Ayyanathan
% This code generates the training and testing data used in the work: 
% Praveen J. Ayyanathan}, Tyler Flegel, Ehsan Taheri, Neural Network Based 
% Collision-Free Trajectory Generator for Quadcopters, 2025 Scitech Forum, 
% Orlando, FL, January 6-10, 2025.

clear
close all
clc

% B-Spline Parameters
order = 4;
numPieces = 4;

% Successive Convex Optimization Parameters
N=101; % Number of nodes
i_max = 20;  
alpha = 1.5;
H = diag([1 1 0]);  % Circular Obstacle Avoidance Region
w = 1e-5;
r_tol = 0.05;  % Stopping criteria for Suceesive Convexification
H1 = 1.5; % Height of obstacle1 for plotting
H2 = 1.5; % Height of obstacle1 for plotting
umin = 8.5; % acceleration/control min 
umax = (9.218251*4)/2.28; % acceleration/control max

% Hover Conditions
vi = [0 0 0];
vf = [0 0 0];

% Max Tilt Angle for UAV
theta_max = 40*pi/180;

% Visualization flags
plotIterFlag = 0;
plotFlag = 1;
plotFlagSpline = 0;
animationFlag = 0;
animation_time = 0;
videoWriteFlag = true;

% Flag data structure
probFlag.plotIterFlag = plotIterFlag;
probFlag.plotFlag = plotFlag;
probFlag.plotFlagSpline = plotFlagSpline;
probFlag.animationFlag = animationFlag;
probFlag.animation_time = animation_time;
probFlag.videoWriteFlag = videoWriteFlag;

% Parameter data structure
probParam.N = N;
probParam.alpha = alpha;
probParam.i_max = i_max;
probParam.H = H;
probParam.w = w;
probParam.r_tol = r_tol;
probParam.H1 = H1;
probParam.H2 = H2;
probParam.umin = umin;
probParam.umax = umax;
probParam.vi = vi;
probParam.vf = vf;
probParam.theta_max = theta_max;

% B-Spline parameter data structure
probFit.order = order;
probFit.numPieces = numPieces;

% Total training data
numSamples = 500000;  % Adjust based on your requirements

% Pre-Allocating array sizes
mseXArr = zeros(numSamples,1);
mseYArr = zeros(numSamples,1);
rmseArr = zeros(numSamples,1);
endTimeArr = zeros(numSamples,1);
eflagArr = cell(numSamples,1);
errPosArr = zeros(numSamples,1);

% Generate andom initial positions
x_i_samples = -2 + (4 * rand(numSamples, 1));  % [-2, 2]
y_i_samples = -4 + (1 * rand(numSamples, 1));  % [-4, -3]

% Generate random final positions
x_f_samples = -2 + (4 * rand(numSamples, 1));  % [-2, 2]
y_f_samples = 3.5 + (0.5 * rand(numSamples, 1));  % [3.5, 4]

% Generate random obstacle positions
x_o_1_samples = -1 + (2 * rand(numSamples, 1));  % [-1, 1]
y_o_1_samples = -1.25 + (0.5 * rand(numSamples, 1));  % [-1.25, -0.75]

% Generate random obstacle positions
x_o_2_samples = -1 + (2 * rand(numSamples, 1));  % [-1, 1]
y_o_2_samples = 1.5 + (0.25 * rand(numSamples, 1));  % [1.5, 1.75]

% Generate random final times
t_f_samples = 7 + (2 * rand(numSamples, 1));  % [7, 9]

% Generate random radius
R_samples = 0.4 + (0.15 * rand(numSamples, 1));  % [0.4, 0.55]

% Neural Network input dimension
input_dim = 11;

% B-Spline Coefficient dimension 
coeffs_dim = 7;

% Neural Network parameter dimension
waypoint_dim = 5;

% Pre-Allocating array sizes
inputs = zeros(numSamples, input_dim);
errFlag = zeros(numSamples, 1);

% Neural Network Trajectory output
outputs_coeffs_x = zeros(numSamples, coeffs_dim);
outputs_coeffs_y = zeros(numSamples, coeffs_dim);

% Neural Network Parameter output
outputs_waypoint_closest_1 = zeros(numSamples, waypoint_dim);
outputs_waypoint_closest_2 = zeros(numSamples, waypoint_dim);

% Neural Network Segment input
inputsParTraj1 = zeros(numSamples, 7); 
inputsParTraj2 = zeros(numSamples, 9);
inputsParTraj3 = zeros(numSamples, 7);

% Neural Network Segment output
outputsParTraj1 = zeros(numSamples, 14); 
outputsParTraj2 = zeros(numSamples, 14); 
outputsParTraj3 = zeros(numSamples, 14); 


ii = 1;
iii = 1;
for i = 1:numSamples

    % Boundary conditions for Successive Convex Optimization
    ri = [x_i_samples(i), y_i_samples(i), 0];
    rf = [x_f_samples(i), y_f_samples(i), 0];
    obsCenter1 = [x_o_1_samples(i), y_o_1_samples(i), 0];
    obsCenter2 = [x_o_2_samples(i), y_o_2_samples(i), 0];
    tf = t_f_samples(i);
    R1 = R_samples(i);
    R2 = R_samples(i);
    
    % Add boundary conditions to the parameter data structure
    probParam.R1 = R1;
    probParam.R2 = R2;
    probParam.ri = ri;
    probParam.rf = rf;
    probParam.obsCenter1 = obsCenter1;
    probParam.obsCenter2 = obsCenter2;
    probParam.tf = tf;

    % Formulate the Convex problem
    [mse_x, mse_y, rmse, eflag_new, err_pos, end_time, dataTransfer] = Main_Obs_Avoid_cvx(probParam,probFlag, probFit);
    
    % Check for convergence before saving data
    % For ii<50001 look for "curved" trajectories
    if strcmp(eflag_new, 'Solved') && dataTransfer.noObsTraj ==0 && ii<50001 && dataTransfer.noSuffPt ==0

        % Mean square and root mean square errors between the convex and B-Spline fit trajectories
        % in x- and y- direction
        mseXArr(ii) = mse_x;
        mseYArr(ii) = mse_y;
        rmseArr(ii) = rmse;
        eflagArr{ii} = eflag_new;
        errPosArr(ii) = err_pos;
        endTimeArr(ii) = end_time;
         
        % Input to the Trajectory and Parameter Neural Network models
        inputs(ii, :) = [ri(1:2), rf(1:2), obsCenter1(1:2), R1,obsCenter2(1:2), R2, tf]; 

        % Output to the Trajectory Neural Network model
        outputs_coeffs_x(ii, :) = dataTransfer.coeffs_x;
        outputs_coeffs_y(ii, :) = dataTransfer.coeffs_y;

        % Output to the Parameter Neural Network model
        outputs_waypoint_closest_1(ii, :) = [dataTransfer.xWpClosest1, dataTransfer.yWpClosest1, dataTransfer.vxWpClosest1, dataTransfer.vyWpClosest1, dataTransfer.tWpClosest1];
        outputs_waypoint_closest_2(ii, :) = [dataTransfer.xWpClosest2, dataTransfer.yWpClosest2, dataTransfer.vxWpClosest2, dataTransfer.vyWpClosest2, dataTransfer.tWpClosest2];
        
        % Input and output to the Segment Neural Network model
        inputsParTraj1(ii, :) = [ri(1:2), dataTransfer.xWpClosest1, dataTransfer.yWpClosest1, dataTransfer.vxWpClosest1, dataTransfer.vyWpClosest1, dataTransfer.tWpClosest1];
        outputsParTraj1(ii,:) = [dataTransfer.XIniObs1Trajcoeffs_x; dataTransfer.XIniObs1Trajcoeffs_y];
        inputsParTraj2(ii, :) = [dataTransfer.xWpClosest1, dataTransfer.yWpClosest1, dataTransfer.vxWpClosest1, dataTransfer.vyWpClosest1,dataTransfer.xWpClosest2, dataTransfer.yWpClosest2, dataTransfer.vxWpClosest2, dataTransfer.vyWpClosest2, dataTransfer.tWpClosest2-dataTransfer.tWpClosest1];
        outputsParTraj2(ii,:) = [dataTransfer.XObs1Obs2Trajcoeffs_x; dataTransfer.XObs1Obs2Trajcoeffs_y];
        inputsParTraj3(ii, :) = [dataTransfer.xWpClosest2, dataTransfer.yWpClosest2, dataTransfer.vxWpClosest2, dataTransfer.vyWpClosest2,rf(1:2), tf];
        outputsParTraj3(ii,:) = [dataTransfer.XObs2EndTrajcoeffs_x; dataTransfer.XObs2EndTrajcoeffs_y];
        ii = ii+1;

    % For ii>50000 look for any trajectory
    elseif strcmp(eflag_new, 'Solved') && ii>50000 && dataTransfer.noSuffPt ==0
        % Mean square and root mean square errors between the convex and B-Spline fit trajectories
        % in x- and y- direction
        mseXArr(ii) = mse_x;
        mseYArr(ii) = mse_y;
        rmseArr(ii) = rmse;
        eflagArr{ii} = eflag_new;
        errPosArr(ii) = err_pos;
        endTimeArr(ii) = end_time;

        % Input to the Trajectory and Parameter Neural Network models
        inputs(ii, :) = [ri(1:2), rf(1:2), obsCenter1(1:2), R1,obsCenter2(1:2), R2, tf]; 
        outputs_coeffs_x(ii, :) = dataTransfer.coeffs_x;
        outputs_coeffs_y(ii, :) = dataTransfer.coeffs_y;
        
        % Output to the Trajectory Neural Network model
        outputs_waypoint_closest_1(ii, :) = [dataTransfer.xWpClosest1, dataTransfer.yWpClosest1, dataTransfer.vxWpClosest1, dataTransfer.vyWpClosest1, dataTransfer.tWpClosest1];
        outputs_waypoint_closest_2(ii, :) = [dataTransfer.xWpClosest2, dataTransfer.yWpClosest2, dataTransfer.vxWpClosest2, dataTransfer.vyWpClosest2, dataTransfer.tWpClosest2];
        
        % Input and output to the Segment Neural Network model
        inputsParTraj1(ii, :) = [ri(1:2), dataTransfer.xWpClosest1, dataTransfer.yWpClosest1, dataTransfer.vxWpClosest1, dataTransfer.vyWpClosest1, dataTransfer.tWpClosest1];
        outputsParTraj1(ii,:) = [dataTransfer.XIniObs1Trajcoeffs_x; dataTransfer.XIniObs1Trajcoeffs_y];
        inputsParTraj2(ii, :) = [dataTransfer.xWpClosest1, dataTransfer.yWpClosest1, dataTransfer.vxWpClosest1, dataTransfer.vyWpClosest1,dataTransfer.xWpClosest2, dataTransfer.yWpClosest2, dataTransfer.vxWpClosest2, dataTransfer.vyWpClosest2, dataTransfer.tWpClosest2-dataTransfer.tWpClosest1];
        outputsParTraj2(ii,:) = [dataTransfer.XObs1Obs2Trajcoeffs_x; dataTransfer.XObs1Obs2Trajcoeffs_y];
        inputsParTraj3(ii, :) = [dataTransfer.xWpClosest2, dataTransfer.yWpClosest2, dataTransfer.vxWpClosest2, dataTransfer.vyWpClosest2,rf(1:2), tf];
        outputsParTraj3(ii,:) = [dataTransfer.XObs2EndTrajcoeffs_x; dataTransfer.XObs2EndTrajcoeffs_y];
        ii = ii+1;
    else
        % Record non-convergence
        errFlag(iii,1)= i;
        iii=iii+1;
    end
    if ii > 55000
        break
    end
end

numValidSamples = ii - 1;

% Trim arrays to include only valid samples
inputs = inputs(1:numValidSamples, :);
outputs_coeffs_x = outputs_coeffs_x(1:numValidSamples, :);
outputs_coeffs_y = outputs_coeffs_y(1:numValidSamples, :);
outputs_waypoint_closest_1 = outputs_waypoint_closest_1(1:numValidSamples, :);
outputs_waypoint_closest_2 = outputs_waypoint_closest_2(1:numValidSamples, :);

% Set the random seed for reproducibility 
rng(42);

% Generate a random permutation of indices
indices = randperm(numValidSamples);

% Define the split ratio (e.g., 80% training, 20% testing)
trainRatio = 0.8;
numTrain = floor(trainRatio * numValidSamples);

% Split the indices
trainIndices = indices(1:numTrain);
testIndices = indices(numTrain+1:end);

% Split inputs for the Trajectory and Parameter Neural Network models
inputs_train = inputs(trainIndices, :);
inputs_test = inputs(testIndices, :);

% Split outputs for the Trajectory and Parameter Neural Network models
outputs_coeffs_x_train = outputs_coeffs_x(trainIndices, :);
outputs_coeffs_x_test = outputs_coeffs_x(testIndices, :);

outputs_coeffs_y_train = outputs_coeffs_y(trainIndices, :);
outputs_coeffs_y_test = outputs_coeffs_y(testIndices, :);

outputs_waypoint_closest_train_1 = outputs_waypoint_closest_1(trainIndices, :);
outputs_waypoint_closest_test_1 = outputs_waypoint_closest_1(testIndices, :);

outputs_waypoint_closest_train_2 = outputs_waypoint_closest_2(trainIndices, :);
outputs_waypoint_closest_test_2 = outputs_waypoint_closest_2(testIndices, :);

% Split the inputs and outputs for the Segment Neural Network
inputsParTraj1_train = inputsParTraj1(trainIndices, :);
inputsParTraj2_train = inputsParTraj2(trainIndices, :);
inputsParTraj3_train = inputsParTraj3(trainIndices, :);

outputsParTraj1_train = outputsParTraj1(trainIndices, :);
outputsParTraj2_train = outputsParTraj2(trainIndices, :);
outputsParTraj3_train = outputsParTraj3(trainIndices, :);

inputsParTraj1_test = inputsParTraj1(testIndices, :);
inputsParTraj2_test = inputsParTraj2(testIndices, :);
inputsParTraj3_test = inputsParTraj3(testIndices, :);

outputsParTraj1_test = outputsParTraj1(testIndices, :);
outputsParTraj2_test = outputsParTraj2(testIndices, :);
outputsParTraj3_test = outputsParTraj3(testIndices, :);

% Save training data
save('training_data_2obs_2_edit.mat', 'inputs_train', 'outputs_coeffs_x_train', 'outputs_coeffs_y_train', ...
    'outputs_waypoint_closest_train_1', 'outputs_waypoint_closest_train_2');

save('training_data_parTraj_2_edit.mat','inputsParTraj1_train','inputsParTraj2_train','inputsParTraj3_train','outputsParTraj1_train',...
    'outputsParTraj2_train', 'outputsParTraj3_train');

% Save testing data
save('testing_data_2obs_2_edit.mat', 'inputs_test', 'outputs_coeffs_x_test', 'outputs_coeffs_y_test', ...
    'outputs_waypoint_closest_test_1', 'outputs_waypoint_closest_test_2');

save('testing_data_parTraj_2_edit.mat','inputsParTraj1_test','inputsParTraj2_test','inputsParTraj3_test','outputsParTraj1_test',...
    'outputsParTraj2_test', 'outputsParTraj3_test');

