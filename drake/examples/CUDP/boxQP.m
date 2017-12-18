function [x,Pfree,Lfree,error] = boxQP(Q,f,lb,ub,x0)
% Minimize 0.5*x'*Q*x + f'*x  s.t. lb<=x<=ub using a projected Newton's
% method. In addition to the solution, it returns the projection matrix
% that projects onto the free subspace, as well as the Cholesky factor of
% the free subspace Hessian.
%
% Adapted from:
% Y. Tassa, N. Mansard, and E. Todorov, "Control-Limited Differential
% Dynamic Programming", ICRA 2014


% Algorithm parameters
maxIter        = 100;       % maximum number of iterations
minGrad        = 1e-8;      % minimum norm of non-fixed gradient
minRelImprove  = 1e-8;      % minimum relative improvement
stepDec        = 0.6;       % factor for decreasing stepsize
minStep        = 1e-16;     % minimal stepsize for linesearch
Armijo         = 0.1;   	% Armijo parameter (fraction of linear improvement required)

% Function that clamps a vector within the bounds
function xc = clamp(x)
    xc = max(lb, min(ub, x));
end

% ------------------------ Problem Setup ------------------------ %

n       = size(Q,1);
In      = eye(n);
Pfree   = In;
Lfree   = zeros(n);
free    = zeros(n,1);
clamped = zeros(n,1);
error   = 0;

% Make sure initial guess is feasible
if nargin > 5 && numel(x0) == n
    x = clamp(x0);
else
    x = clamp(zeros(n,1));
end

% initial objective value
val = 0.5*x'*Q*x + f'*x;

% ------------------------ Solver Loop ------------------------ %
for iter = 1:maxIter
    
    % Unconstrained gradient
    g = Q*x + f;
    
    % Find free and clamped dimensions
    old_free = free;
    clamped((x == lb)&(g>0)) = 1;
    clamped((x == ub)&(g<0)) = 1;
    free = ~clamped;
    
    % Projection onto free subspace
    if(any(free))
        Pfree = In(free,:);
    else
        Pfree = 0;
        break; % Stop if everything is clamped
    end
    
    % Free subspace Hessian
    if free ~= old_free
        Qff = Pfree*Q*Pfree';
        [Lfree, error] = chol(Qff, 'lower');
        if error
            break; % Stop if Hessian isn't positive-definite
        end
    end
    
    % Stop if the free subspace gradient norm is below the tolerance
    gnorm  = norm(g(free));
    if gnorm < minGrad
        break;
    end
    
    % Calculate descent direction (minimize over free subspace)
    d = -Pfree'*Lfree'\(Lfree\(g(free)));
    
    % Armijo line search
    step  = 1;
    nstep = 0;
    xc    = clamp(x+step*d);
    valc    = 0.5*xc'*Q*xc + f'*xc;
    expected = d'*g;
    while (valc - val)/(step*expected) < Armijo
        step  = step*stepDec;
        nstep = nstep+1;
        xc    = clamp(x+step*d);
        valc  = 0.5*xc'*Q*xc + f'*xc;
        if step<minStep
            break;
        end
    end
    
    % Check improvement
    if(val - valc < minRelImprove*abs(val))
        x = xc;
        break;
    else
        val = valc;
        x = xc;
    end
end

end