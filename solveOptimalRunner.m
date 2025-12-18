function out = solveOptimalRunner(xi,vi, xf,vf, opts)
% solve a single minimum-time segment with |a(t)|<=1
% between fixed endpoints, with known start/end velocities.
%
% inputs:
%   xi, xf : 1x2 positions in world coords
%   vi, vf : 1x2 velocities in world coords (same time units as output)
%   opts : struct with optional fields:
%     N (default 200)
%     FL_init initial guess [a b c S]
%     iters (default 20) Newton steps in MatchCoords
%
% output struct fields:
%   X, Y : world coordinates
%   T : segment time S
%   FLC : [a b c S]
%   vframe_i : start velocity in frame
%   vframe_f : end velocity in frame

  if nargin<5, opts = struct; end
  if ~isfield(opts,'N'), opts.N = 200; end
  if ~isfield(opts,'iters'), opts.iters = 20; end
  if ~isfield(opts,'FL_init'), opts.FL_init = [0.5 0.0 0.5 1.2]; end

  % Globals for BVP
  clear global N L
  global N L
  N = opts.N;
  M = -2*diag(ones(N-1,1)) + diag(ones(N-2,1),1) + diag(ones(N-2,1),-1);
  L = inv(M);

  % Build segment frame: translate xi to origin; rotate so xf-xi aligns +x; scale by length Lseg
  e = (xf - xi);
  Lseg = norm(e);
  if Lseg==0, error('xi and xf must differ.'); end
  ex = e / Lseg;
  R  = [ex; [-ex(2) ex(1)]];      % rows are ex and ey (CCW)
  % World->frame: p' = [(p - xi)*R'] / Lseg
  % Velocities rotate but do NOT scale with Lseg (space is scaled, not time).
  vi_f = (vi*R');   vi_f = vi_f / Lseg;
  vf_f = (vf*R');   vf_f = vf_f / Lseg;

  % Target velocity components (alph1,bet1 at start; alph2,bet2 at end)
  Target = [vi_f(1) vi_f(2) vf_f(1) vf_f(2)];

  % Solve for [a,b,c,S] in the canonical frame (0,0)->(1,0)
  FLC = MatchCoords(opts.FL_init, Target, opts.iters);
  S   = FLC(4);

  % Recover canonical path, then map back to world
  [~,~,~,~, Xc, Yc] = abcSfun(FLC);         % 0..1 path, N+1 points
  % Frame->world: p = xi + Lseg * [Xc Yc] * R
  P = [Xc(:) Yc(:)] * R;
  Xw = xi(1) + Lseg * P(:,1);
  Yw = xi(2) + Lseg * P(:,2);

  out = struct('X',Xw(:).', 'Y',Yw(:).', 'T',S, 'FLC',FLC, ...
               'vframe_i',vi_f, 'vframe_f',vf_f);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% === Core routines from the paper (lightly cleaned but unchanged)      ===
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [alph1, bet1, alph2, bet2, X, Y ] = abcSfun(FLC)
% Given linear accel field (At+B) normalized to |a(t)|=1 and time S,
% solve path from (0,0) to (1,0) with FD BVP.
  global N L
  a = FLC(1); b = FLC(2); c = FLC(3); d = 1.0; S = FLC(4);

  T = linspace(0,S,N+1);
  D = sqrt((a*T+b).^2 + (c*T+d).^2);
  F = (S/N)^2*(c*T+d)./D;
  F(N) = F(N)-1;
  G = (S/N)^2*(a*T+b)./D;

  X = zeros(1,N+1); Y = zeros(1,N+1);
  X(1)=0; X(N+1)=1; X(2:N)= (L * F(2:N)')';
  Y(1)=0; Y(N+1)=0; Y(2:N)= (L * G(2:N)')';

  alph1 = (X(2)-X(1))*N/S;
  alph2 = (X(N+1)-X(N))*N/S;
  bet1  = (Y(2)-Y(1))*N/S;
  bet2  = (Y(N+1)-Y(N))*N/S;
end

function J = abcSjako(FLC)
  del = .1;
  [a1 b1 a2 b2] = abcSfun(FLC-[del 0 0 0]);
  [c1 d1 c2 d2] = abcSfun(FLC+[del 0 0 0]);
  J(1:4,1) = ([c1 d1 c2 d2]-[a1 b1 a2 b2])/(2*del);

  [a1 b1 a2 b2] = abcSfun(FLC-[0 del 0 0]);
  [c1 d1 c2 d2] = abcSfun(FLC+[0 del 0 0]);
  J(1:4,2) = ([c1 d1 c2 d2]-[a1 b1 a2 b2])/(2*del);

  [a1 b1 a2 b2] = abcSfun(FLC-[0 0 del 0]);
  [c1 d1 c2 d2] = abcSfun(FLC+[0 0 del 0]);
  J(1:4,3) = ([c1 d1 c2 d2]-[a1 b1 a2 b2])/(2*del);

  [a1 b1 a2 b2] = abcSfun(FLC-[0 0 0 del]);
  [c1 d1 c2 d2] = abcSfun(FLC+[0 0 0 del]);
  J(1:4,4) = ([c1 d1 c2 d2]-[a1 b1 a2 b2])/(2*del);
end

function FLCoordsMatch = MatchCoords(FLCoords,Target,iter)
  % Damped Gauss-Newton (Levenberg–Marquardt) with backtracking
  lambda = 1e-2;    % initial damping
  for i = 1:iter
    [a1,b1,a2,b2] = abcSfun(FLCoords);
    r = [Target(1)-a1; Target(2)-b1; Target(3)-a2; Target(4)-b2];

    J = abcSjako(FLCoords);

    % Solve (J'J + lambda I) * d = J' r
    A = J.'*J + lambda*eye(4);
    g = J.'*r;

    % If A is ill-conditioned, increase damping
    if rcond(A) < 1e-10
      lambda = max(lambda*10, 1e-6);
      A = J.'*J + lambda*eye(4);
    end
    d = A \ g;

    % Line-search on step size
    t = 1.0;
    base_err = norm(r);
    accepted = false;
    while t >= 1/64
      trial = FLCoords + t*d.';
      % keep time positive
      if trial(4) <= 1e-6, trial(4) = 1e-6; end

      [ta1,tb1,ta2,tb2] = abcSfun(trial);
      r_new = [Target(1)-ta1; Target(2)-tb1; Target(3)-ta2; Target(4)-tb2];
      if norm(r_new) < base_err
        FLCoords = trial;
        lambda   = max(lambda/2, 1e-8);   % decrease damping if progress
        accepted = true;
        break
      else
        t = t/2;                          % backtrack
      end
    end
    if ~accepted
      lambda = min(lambda*4, 1e3);        % increase damping if stuck
    end

    % early exit if residual small
    if norm(r) < 1e-8, break; end
  end
  FLCoordsMatch = FLCoords;
end