function [c, ceq] = runnerPathCst2(t, x, u, params)
% x = [px; py; vx; vy; z]
% u = [ax; ay]

ax = u(1,:);
ay = u(2,:);

% --- Inequality: accel magnitude <= aMax ---
a2      = ax.^2 + ay.^2;
%c_accel = a2 - params.aMax^2;   % <= 0

% No additional path equalities now (slacks removed)
c   = [];
ceq = [];
end
