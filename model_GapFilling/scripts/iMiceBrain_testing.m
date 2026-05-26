% testing iMiceBrain functionality against metabolic tasks 
initCobraToolbox('false');
changeCobraSolver('gurobi', 'all');
setRavenSolver('gurobi');
model = readCbModel('/Users/eso1993/Library/CloudStorage/Box-Box/PFOS_Project/drafts_statistics_AND_consensus_draft/models/final_models/iMiceBrain.mat');

%% adopting from ftINIT algorithm in testing model
% prepINITModel needs a rev domain in model structure 
if ~isfield(model, 'rev')
    model.rev = model.lb < 0 & model.ub > 0;
end

% perform testing
taskStruct = parseTaskList('/Users/eso1993/Library/CloudStorage/Box-Box/Papers preparations/CP2 paper/supplmentary files/metabolicTasks_Essential_mouse_model.xls');
closedModel = closeModel(model);
[taskReport, essentialRxns, taskStructure, essentialFluxes]=checkTasks(closedModel,{},{},{},{},taskStruct);
%% fastFVA 
[minFlux, maxFlux, optsol, ret, fbasol, fvamin, fvamax, statussolmin, statussolmax] = fastFVA (model, 100, 'max', 'ibm_cplex', model.rxns, model.S, '', '', '');
