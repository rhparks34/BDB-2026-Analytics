function xDot = runnerDynamics_simple4(t, x, u, params)
% x = [px; py; vx; vy; z]
% u = [F; lambda]
%
% F      : power-like target (instantaneous effort level)
% lambda : tradeoff in [0,1] between
%          - lambda = 0 => aDist: accel along velocity (max-speed option)
%          - lambda = 1 => aAlign: pure lateral accel to turn toward goal
%
% Acceleration direction aDir is a blend of these two; magnitude r is chosen
% along aDir so that P(r*aDir) ≈ F using the power model P = P_tan + P_lat.

px = x(1,:); 
py = x(2,:);
vx = x(3,:);
vy = x(4,:);
z  = x(5,:); %#ok<NASGU>

F_ctrl      = u(1,:);   % target power
lambda_ctrl = u(2,:);   % tradeoff

N = size(x,2);

% Time derivatives
pDotX = vx;
pDotY = vy;
vDotX = zeros(1,N);
vDotY = zeros(1,N);
zDot  = zeros(1,N);

% Parameters
m       = params.m;
k       = params.k;
gamma   = params.gamma;
vEps    = params.vEps;
eps_abs = 1e-2;                   % smoothing in P_tan

pTarget = params.pTarget(:);      % [px_T; py_T]

fzeroOpts = optimset('Display','off');   % keep inner solver quiet

for j = 1:N

    Fj      = F_ctrl(j);
    lambdaj = lambda_ctrl(j);

    % --- Velocity and unit directions ---
    vxj = vx(j);
    vyj = vy(j);
    v2  = vxj^2 + vyj^2;
    vNorm = sqrt(v2 + vEps^2);

    % Direction to goal
    p    = [px(j); py(j)];
    g    = pTarget - p;
    gNorm = norm(g);

    if gNorm < 1e-8
        gHat = [1; 0];            % arbitrary if at target
    else
        gHat = g / gNorm;
    end

    % Unit velocity: if very small speed, just align with goal
    if vNorm < 1e-8
        vHat = gHat;
    else
        vHat = [vxj; vyj] / vNorm;
    end

    % Normal vector (rotate vHat by +90 degrees)
    n = [-vHat(2); vHat(1)];

    % Decompose goal direction in (vHat, n) basis
    c_v = dot(gHat, vHat); %#ok<NASGU>  % not used directly here, but kept for clarity
    c_n = dot(gHat, n);

    % -------------------------------
    % 1) Define the two canonical directions
    % -------------------------------


    % (a) Distance / "momentum" direction along the tangent
    % If velocity is helping (c_v >= 0), accelerate along vHat.
    % If velocity is hurting (c_v < 0), brake along -vHat to kill momentum.
    if c_v >= 0
        dDist = vHat;        % preserve/use momentum
    else
        dDist = -vHat;       % maximum negative tangential accel (brake)
    end

    % (a) Distance / max-speed direction: purely along velocity
    %dDist = vHat;   % unit, positive tangential, zero lateral

    % (b) Alignment direction: purely lateral, rotating toward goal
    % (b) Alignment direction: max turning for given budget Fj,
    % with tangential accel chosen to cancel drag in the forceTan term:
    %   m*aTan + k*v^2 = 0  => aTan_cancel = -(k/m)*v^2
    if abs(c_n) < 1e-8
        dAlign = dDist;  % already aligned in this frame
    else
        nHat = sign(c_n) * n;        % unit normal toward goal-side
    
        aTan_cancel = -(k/m) * v2;   % drag-cancel tangential accel (scalar)
    
        % All remaining budget goes to lateral turning:
        % P_lat = gamma*m*||aNorm||^2  => aNorm_max = sqrt(Fj/(gamma*m))
        aNorm_max = sqrt(max(Fj,0) / (gamma*m));
    
        aAlign_vec = aTan_cancel * vHat + aNorm_max * nHat;
    
        if norm(aAlign_vec) < 1e-12
            dAlign = nHat;
        else
            dAlign = aAlign_vec / norm(aAlign_vec);  % direction for blending
        end
    end


    % -------------------------------
    % 2) Blend directions via lambda
    % -------------------------------
    %lambda_abs = abs(lambdaj);
    lambdaj = max(0, min(1, lambdaj));  % clamp to [0,1]
    %lambda_abs = abs()
    dRaw = (1 - lambdaj)*dDist + lambdaj*dAlign;
    if norm(dRaw) < 1e-8
        dRaw = dDist;
    end
    aDir = dRaw / norm(dRaw);    % unit direction of acceleration

    % -------------------------------
    % 3) Solve for magnitude r along aDir so P(r*aDir) ≈ F
    % -------------------------------
    if Fj <= 0
        % No effort: zero acceleration
        rBest = 0;
        P_val = power_cost(0, 0, vxj, vyj, m, k, gamma, vEps, eps_abs);
    else
        % Define scalar function f(r) = P(r*aDir) - F
        fun = @(r) power_cost(r*aDir(1), r*aDir(2), vxj, vyj, ...
                              m, k, gamma, vEps, eps_abs) - Fj;

        % Evaluate at r = 0
        f0 = fun(0);

        if f0 >= 0
            % Even zero accel is already too "expensive" for this F:
            % best we can do is r = 0.
            rBest = 0;
            P_val = f0 + Fj;  % fun(0) + Fj = P(0)
        else
            % Find an upper bracket where P(r_hi) > F
            r_hi = 1.0;
            f_hi = fun(r_hi);

            iter = 0;
            while f_hi <= 0 && iter < 20
                r_hi = 2*r_hi;
                f_hi = fun(r_hi);
                iter = iter + 1;
            end

            if f_hi <= 0
                % Couldn't bracket root (F extremely large vs P), just take r_hi
                rBest = r_hi;
                P_val = power_cost(rBest*aDir(1), rBest*aDir(2), ...
                                   vxj, vyj, m, k, gamma, vEps, eps_abs);
            else
                % Use fzero on [0, r_hi]
                rBest = fzero(fun, [0, r_hi], fzeroOpts);
                rBest = max(rBest, 0);
                P_val = power_cost(rBest*aDir(1), rBest*aDir(2), ...
                                   vxj, vyj, m, k, gamma, vEps, eps_abs);
            end
        end
    end

    aFinal = rBest * aDir;
    ax = aFinal(1);
    ay = aFinal(2);

    % Dynamics
    vDotX(j) = ax;
    vDotY(j) = ay;
    zDot(j)  = P_val;   % power-like cost

end

xDot = [pDotX;
        pDotY;
        vDotX;
        vDotY;
        zDot];

end

% =====================================================================
% Helper: your original P = P_tan + P_lat for a single (ax, ay)
% =====================================================================
function P = power_cost(ax, ay, vx, vy, m, k, gamma, vEps, eps_abs)

    v2    = vx^2 + vy^2;
    denom = sqrt(v2 + vEps^2);

    if denom < 1e-8
        vHatX = 1; vHatY = 0;
    else
        vHatX = vx / denom;
        vHatY = vy / denom;
    end

    % Tangential component
    aTan   = ax * vHatX + ay * vHatY;   % scalar tangential
    aTanX  = aTan * vHatX;
    aTanY  = aTan * vHatY;

    % Normal components
    aNormX = ax - aTanX;
    aNormY = ay - aTanY;

    % Tangential "force" term
    forceTan = m * aTan + k * v2;
    P_tan    = sqrt(forceTan^2 + eps_abs^2);

    % Lateral penalty
    P_lat = gamma * m * (aNormX^2 + aNormY^2);

    P = P_tan + P_lat;
end