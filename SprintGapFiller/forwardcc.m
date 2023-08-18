function [flux,z] = forwardcc(model,core,tol)
% model must have fields 'S' and 'rev'
% core is a binary vector indicating whether core or non-core reaction
[m,n] = size(model.S);
n_core = sum(core);% number of core reactions

% objective
f = -1*[zeros(n,1);unifrnd(1,1.1,n_core,1)];

% equalities
Aeq = [model.S, sparse(m,n_core)];
beq = zeros(m,1);
csenseeq = repmat('E',m,1); % equality

% inequalities
temp1 = speye(n);
temp2 = speye(n_core);
Aineq = [temp1(core,:),-1*temp2 ];
bineq = zeros(n_core,1);
csenseineq = repmat('G',n_core,1); % greater than

% bounds
lb = model.lb;
ub = model.ub;
lb = [lb;-Inf(n_core,1)];
ub = [ub; tol*ones(n_core,1)];

% Set up LP problem
LPproblem.A=[Aeq;Aineq];
LPproblem.b=[beq;bineq];
LPproblem.lb=lb;
LPproblem.ub=ub;
LPproblem.c=f;
LPproblem.osense=1;%minimise
LPproblem.csense = [csenseeq; csenseineq];

solution = solveCobraLP(LPproblem);
if solution.stat~=1
    fprintf('%s%s\n',num2str(solution.stat),' = solution.stat')
    fprintf('%s%s\n',num2str(solution.origStat),' = solution.origStat')
    warning('LP solution may not be optimal')
end
x=solution.full;

if ~isempty(x)
    flux = x(1:n);
    z = x(n+1:end);
else
    flux=nan(n,1);
    z=nan(n_core,1);
end

