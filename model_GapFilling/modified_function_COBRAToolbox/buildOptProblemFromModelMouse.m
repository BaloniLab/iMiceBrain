function optProblem = buildOptProblemFromModelMouse(model)
% modified from buildOptProblemFromModel.m
% Builds a COBRA Toolbox LP, QP, RLP, or RQP problem structure from a COBRA Toolbox model structure.
%
% USAGE:
%    optProblem = buildOptProblemFromModel(model)
%
% INPUT:
%    model:     A COBRA model structure with at least the following fields:
%                - `.S` - The stoichiometric matrix
%                - `.c` - Objective coefficient vector
%                - `.lb` - Lower bound vector
%                - `.ub` - Upper bound vector              
%                - `.b` - Right-hand side vector (accumulation/depletion)
%                - `.csense` - Constraint sense ('E' for equality, 'L' for ?, 'G' for ?)
%
% OUTPUT:
%    optProblem: A COBRA optProblem structure with:
%                - `.A`: LHS matrix (stoichiometric constraints)
%                - `.b`: RHS vector (accumulation/depletion)
%                - `.c`: Objective coefficients
%                - `.lb`: Lower bound vector
%                - `.ub`: Upper bound vector
%                - `.osense`: Objective sense (`-1` for maximize, `1` for minimize)
%                - `.csense`: Constraint sense
%
%    The function has been modified to remove optional inputs for verification and parameters.

[nMet, nRxn] = size(model.S);

% Ensure required fields exist
if ~isfield(model, 'c')
    model.c = zeros(nRxn, 1);
end

if ~isfield(model, 'b')
    model.b = zeros(nMet, 1);
end

if ~isfield(model, 'csense')
    model.csense = repmat('E', nMet, 1); % Default to equality constraints
end

% Construct optProblem structure
optProblem.A = model.S;
optProblem.b = model.b;
optProblem.ub = model.ub;
optProblem.lb = model.lb;
optProblem.c = model.c;
optProblem.csense = model.csense;

% Set objective sense
[~, optProblem.osense] = getObjectiveSense(model);

end
