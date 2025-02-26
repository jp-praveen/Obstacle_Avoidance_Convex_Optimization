function [sol,fval,eflag] = obstacle_relaxed_cvx_foh(N, tf, g, ri, rf, vi, vf, ai, af, umin, umax, thetaMax)

% Arena
xmin = -50;
xmax = 50;
ymin = xmin;
ymax = xmax;
zmin = -5;
zmax = 5;

% discretised time
t = tf/N;

cvx_begin
    cvx_quiet(true)
    cvx_solver ecos

    % Define variables
    variable r(N,3)
    variable v(N,3)
    variable u(N,3)
    variable s(N)
    variable s1

    % Objective
    minimize(s1)

    subject to
        % S1 constraints
        norm([2*s; 1-s1]) <= (s1+1);
        s1 >= -1;

        % Slack constraints
        
        norms(u, 2, 2) <= s;        % Vectorized norm for each row of u
        s <= umax;                 % Element-wise comparison for s
        s >= umin;                 % Element-wise comparison for s
        s .* cos(thetaMax) <= u(:,3); % Element-wise comparison for each row of u(:,3)

        % Hover constraints on control vec
        u(1,:) == ai;
        u(N,:) == af;

        % Velocity constraints
        v(1,:) == vi;
        v(N,:) == vf;
        v(2:N,:) == v(1:N-1,:) + t * (u(1:N-1,:) + 0.5 * (u(2:N,:) - u(1:N-1,:)) + repmat(g, N-1, 1));

        % Position constraints
        r(1,:) == ri;
        r(2:N,:) == r(1:N-1,:) + t*v(1:N-1,:) + (t^2/2) * (u(1:N-1,:) + 0.5 * (u(2:N,:) - u(1:N-1,:)) + + repmat(g, N-1, 1));
        r(N,:) == rf;

        % Arena Constraints
        
        r(:,1) <= xmax; r(:,1) >= xmin;
        r(:,2) <= ymax; r(:,2) >= ymin;
        r(:,3) <= zmax; r(:,3) >= zmin;

        % if useInitialGuessFlag == 1
        %     r == initialGuessStruct.r_guess;
        %     v == initialGuessStruct.v_guess;
        %     u == initialGuessStruct.u_guess;
        %     s == initialGuessStruct.s_guess;
        %     s1 == initialGuessStruct.s1_guess;
        %  end

cvx_end

% Prepare output
sol.r = r;
sol.v = v;
sol.u = u;
sol.s = s;
sol.s1 = s1;

fval = cvx_optval;
eflag = cvx_status;
output = [];