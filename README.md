# Neural Network Based Collision-Free Trajectory Generator for Quadcopters
Link to Paper: https://arc.aiaa.org/doi/abs/10.2514/6.2025-0565

In this work, I developed a neural network based framework for collision-free trajectory generation for quadcopters navigating in environments with static obstacles. Using trajectories generated through a successive convexification method, three neural network (NN) models: 1) NN-Trajectory, 2) NN-Parameter, and 3) NN-Segment were developed to approximate and reconstruct optimal trajectories for one- and two-obstacle scenarios. Hyperparameter tuning and architecture optimization for the NNs were performed using Optuna ensuring efficient performance. A cascaded structure using NN-Parameter and NN-Segment is proposed to enhance the interpretability of the framework. Additionally, a safety margin correction step was introduced to dynamically adjust position-level information of parameter points (points closest to the obstacles), ensuring collision-free trajectories with negligible loss of accuracy or optimality.
The framework was validated across multiple scenarios demonstrating its capability.

This repository contains to MATLAB code that was used to obtain the training and testing data based on the successive convexification technique outlined in "
Szmuk, M., Pascucci, C. A., Dueri, D., and Açikmeşe, B., Convexification and real-time on-board optimization for agile quad-rotor maneuvering and obstacle avoidance, 2017 IEEE/RSJ International Conference on Intelligent Robots and Systems (IROS), IEEE, 2017, pp. 4862–4868"

The code used for Neural Network hyperparameter optimization and training can be found here: https://github.com/jp-praveen/Neural_Network_Trajectory_Generator