function J = runnerTerminalCost2(t0, x0, tF, xF, targetPos)
% Minimize squared distance from final position to target

pxF = xF(1,:);
pyF = xF(2,:);

J = 2*(pxF - targetPos(1)).^2 + (pyF - targetPos(2)).^2;
end