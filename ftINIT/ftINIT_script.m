%% Constructing the mouse brain-draft, metabolic reconstruction using ftINIT in matlab following this tutorial (https://github.com/SysBioChalmers/RAVEN/blob/main/INIT/ftINIT.m)
% tutorial and documentation (https://sysbiochalmers.github.io/Human-GEM-guide/getting_started/)
% working with ftINIT needs MATLAB vresion of R2019b and later
% Initialize the COBRA Toolbox
initCobraToolbox('false');
changeCobraSolver('gurobi', 'all');
setRavenSolver('gurobi');

%% read the parental model first and test it against metabolic tasks 
model = readCbModel('/Users/eso1993/z/CloudStorage/Box-Box/PFOS_Project/mouse_model_iMM1865/iMM1865_updated.mat');

% ftINIT needs a rev domain in model structure 
if ~isfield(model, 'rev')
    model.rev = model.lb < 0 & model.ub > 0;
end

% test the model against metabolic tasks list 
taskStruct = parseTaskList('/Users/eso1993/Library/CloudStorage/Box-Box/PFOS_Project/ftINIT_mouse_models/input_data_ftINIT/metabolicTasks_Essential_mouse_model.xls');
prepData = prepINITModel(model, taskStruct, {}, false, {}, 'e', true); % 'e' the symbol for extracellular compartment in iMM1865 
prepData.essentialRxns = unique([prepData.essentialRxns; {'BIOMASS_reaction'}]); % define the objctive function as an essential reaction

%% read the brain transcriptomics data and preprocess them
% we will use TPM normalized transcriptomics data 
% read the brain transcriptomics data 
% we will determine the threshold depending on mean value of expression data
brain = readtable('/Users/eso1993/Library/CloudStorage/Box-Box/PFOS_Project/ftINIT_mouse_models/input_data_ftINIT/TPM_normalized_brain_mapped_against_model.csv');
brain_data = table2cell(brain(:, 2:end)); % Exclude the first column (gene IDs)
brain_data = cell2mat(brain_data); % Convert cell array to numeric matrix
brain_data_mean = mean(brain_data, 2); % Calculate row-wise mean

% set the quantiles:
brain_local_thresholds_30 = quantile(brain_data, 0.3, 2); % 30th quantile per gene (local lower threshold)
brain_local_thresholds_90 = quantile(brain_data, 0.9, 2); % 90th quantile per gene (local upper threshold)

% Compute the mean of the local thresholds:
brain_gmin = mean(brain_local_thresholds_30); % Global lower threshold
brain_gmax = mean(brain_local_thresholds_90); % Global upper threshold

% Classify each gene into one of three categories based on gene_expression:
TAS_brain = zeros(size(brain_data_mean)); % Initialize Transcript Activity Scores (TAS)
TAS_brain(brain_data_mean < brain_gmin) = -1; % Inactive
TAS_brain(brain_data_mean > brain_gmax) = 1;  % Active
TAS_brain(brain_data_mean >= brain_gmin & brain_data_mean <= brain_gmax) = 0; % Intermediate



%% set the input for ftINIT
genes = table2cell(brain(:,1)); % either brain or liver or kidney, all have the same first column for gene IDs
genes = cellstr(string(genes));
arrayData.genes = genes; % get the column with gene id
arrayData.tissues = {'Brain', 'Liver', 'Kidney'}; % set the tissues for the models
arrayData.levels = [TAS_brain, TAS_liver, TAS_kidney]; % map the scores
arrayData.threshold = 0; % Set the threshold for active/inactive
paramsFT = 5000;

%% run ftINIT
% run the two steps
% The first step excludes most of the reactions without gene rules (GPRs) from the problem, and the second step determines which of those reactions should be removed
modelBrain = ftINIT(prepData, arrayData.tissues{1}, [], [], arrayData, {}, getHumanGEMINITSteps('1+1'), false, true);
modelBrain = updateGenes(modelBrain); % remove inactive genes
modelBrain.id = 'ftINIT_iMiceBrain'; % change the model ID
modelBrain.name = 'ftINIT_iMiceBrain'; % change the model name
[metbrain] = findMetsFromRxns(modelBrain, modelBrain.rxns);
[metbrain] = regexprep(metbrain, '_(.)$', '[$1]'); % convert the format of mets and compartments
modelBrain.mets = metbrain; % update the mets field within the model

%% export the models as matlab 
outputbrain = '/Users/eso1993/Desktop/ftINIT_iMiceBrain.mat';
save(outputbrain, 'modelBrain');

%% export models in Xml file format
writeCbModel(modelBrain, 'format', 'sbml', 'fileName', '/Users/eso1993/Library/CloudStorage/Box-Box/PFOS_Project/ftINIT_mouse_models/ftINIT_iMiceBrain.xml');
