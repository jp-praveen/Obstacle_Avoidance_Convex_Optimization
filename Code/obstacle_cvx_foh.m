function [sol,fval,eflag] = obstacle_cvx_foh(N, tf, w, obsCenter1,obsCenter2, R1, R2, sol_old, g, ri, rf, vi, vf, ai, af, umin, umax, thetaMax)
rPrev = sol_old.r;

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
    variable nu1
    variable nu2

    % Objective
    minimize(s1*w + nu1 + nu2)

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


        % Obstacle-1 Constraints
        H = diag([1 1 0]);
        Del_r1 = (rPrev - obsCenter1);   % Difference for all rows
        xi_1 = vecnorm(H * Del_r1', 2, 1)'; % Compute norm along each row
        zi_1 = (H' * H * Del_r1') ./ xi_1'; % Compute zi_1 for all rows
        ( xi_1 + sum((zi_1' .* (r - rPrev))')' ) >= (R1 - nu1); % Vectorized constraint for all rows
        
        % Obstacle-1 Constraints
        Del_r2 = (rPrev - obsCenter2);   % Difference for all rows
        xi_2 = vecnorm(H * Del_r2', 2, 1)'; % Compute norm along each row
        zi_2 = (H' * H * Del_r2') ./ xi_2'; % Compute zi_2 for all rows
        ( xi_2 + sum((zi_2' .* (r - rPrev))')' ) >= (R2 - nu2); % Vectorized constraint for all rows

        % Relaxation constraints
        nu1 >= 0;
        nu2 >= 0;

        % if useInitialGuessFlag == 1
        %     r == sol_old.r;
        %     v == sol_old.v;
        %     u == sol_old.u;
        %     s == sol_old.s;
        %     s1 == sol_old.s1;
        %     nu1 == sol_old.nu1;
        %     nu2 == sol_old.nu2;
        %     % cvx_solver_settings('init_val', x0)
        %  end

cvx_end

% Prepare output
sol.r = r;
sol.v = v;
sol.u = u;
sol.s = s;
sol.s1 = s1;
sol.nu1 = nu1;
sol.nu2 = nu2;

fval = cvx_optval;
eflag = cvx_status;
output = [];