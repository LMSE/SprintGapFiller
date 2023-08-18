function ConsReacIDS = speedcc(model,tol)

if nargin < 2 || isempty(tol)
    tol=1e-4; 
end
[m,n] = size(model.S);
IrR = model.ub<=0; % reactions irreversible in reverse direction
temp = model.ub(IrR);

model.S(:,IrR) = -model.S(:,IrR);
model.ub(IrR) = -1*model.lb(IrR);
model.lb(IrR) = -1*temp;
rev = ones(n,1);
rev(model.lb>=0) = 0;
model.rev=rev;
LPS=0;

prev_rxns = false(numel(model.rxns),1);
temp_core = true(numel(model.rxns),1);
while sum(temp_core)~=sum(prev_rxns)
    prev_rxns = temp_core;
    LPS = LPS+2;
    [flux1,~] = forwardcc(model,temp_core,tol);
    temp_core(abs(flux1)>=0.99*tol)=false;
    [flux2,~] = reverse(model,temp_core,tol);
    temp_core(abs(flux2)>=0.99*tol)=false;
end

ConsReacIDS=setdiff([1:numel(model.rxns)],find(temp_core))';