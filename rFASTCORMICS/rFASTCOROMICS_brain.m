%% Building Brain tissue-specific metabolic reconstructions using rFASTCORMICS
% tutorial avaialble at github (https://github.com/sysbiolux/rFASTCORMICS?utm_source=chatgpt.com)
% USES the RNAseq data which are preferably FPKM transformed but can work
% with TPM normalized data
initCobraToolbox('false');
changeCobraSolver('ibm_cplex', 'all'); % it works with cplex solver 

% load the model and create consistent model by running FASTCC (Vlassis et al., 2014)
load('/Users/egabal/Library/CloudStorage/Box-Box/PFOS_Project/mouse_model_iMM1865/iMM1865_updated.mat');
model.description = 'iMM1865';
mouse_model = fastcc_4_rfastcormics(ratGEM, 1e-4,0); % create the consistent model
Cmodel = removeRxns(modez, model.rxns(setdiff(1:numel(model.rxns),mouse_model)));

% read the TPM normalized genes count of brain tissue 
brain = readtable('/Users/eso1993/Library/CloudStorage/Box-Box/PFOS_Project/ftINIT_mouse_models/input_data_ftINIT/TPM_normalized_brain_mapped_against_model.csv');
brain_data = table2cell(brain(:, 2:end)); % Exclude the first column (gene IDs)
brain_data = cell2mat(brain_data); % Convert cell array to numeric matrix
colnames = brain.Properties.VariableNames(:, 2:end); % get the column names which are the sample IDs 

% data discretization
discretized = discretize_FPKM_Mouse(brain_data, colnames);

% identify the genes names from the data file 
rownames = table2cell(brain(:,1)); 
rownames = cellstr(string(rownames));

% get the dictionary for mapping the genes 
dico = readtable('/Users/eso1993/Library/CloudStorage/Box-Box/PFOS_Project/rFASTCOROMICS_mouse_models/input/dico_data_biomRT.csv');
dico.entrezgene_id = cellstr(string(dico.entrezgene_id)); % convert to string of cell array

% to construct one brain-specific generic metabolic network from 93 samples, we will assign the genes 
% a score of 1 has to be found in at > 70% of the samples to be considered
% expressed otherwise non-expressed or uncertain
consensus_proportion = 0.7;

% identify the objective function reaction to be kept 
biomass_rxn = {'BIOMASS_reaction'};
already_mapped_tag = 0;
epsilon = 1e-5; %avoid small number errors

% set reactions that won't be penalized such as that of genes encoding for
% transporter proteins 
unpenalizedSystems = {'Transport, endoplasmic reticular';
    'Transport, extracellular';
    'Transport, golgi apparatus';
    'Transport, mitochondrial';
    'Transport, peroxisomal';
    'Transport, lysosomal';
    'Transport, nuclear'};
Cmodel.subSystems = cellfun(@char, Cmodel.subSystems, 'UniformOutput', false); % ensure that the elemnts are character vector
unpenalized = Cmodel.rxns(ismember(Cmodel.subSystems,unpenalizedSystems));

% initialize the optional settings 
optional_settings = struct();
optional_settings.unpenalized = unpenalized; % assign the unpenalized reactions

% enforce the reaction of the objective function for reconstructed model
optional_settings.func = {'BIOMASS_reaction'}; 

% run rFASTCOROMICS
[modelBrain] = fastcormics_RNAseq_Mouse(Cmodel, discretized, rownames, dico, ...
    biomass_rxn, already_mapped_tag, consensus_proportion, epsilon, optional_settings);

% update the genes of model to remove inactive genes 
modelBrain = updateGenes(modelBrain);
modelBrain.modelID = 'rFASTCORMICS_Mouse_brain';
modelBrain.modelName = 'rFASTCORMICS_Mouse_brain';
modelBrain.description = 'mouse brain-tissue specific metabolic network';

% export the model in matlab format
outputbrain = '/Users/eso1993/Desktop/FASTCORMICS_iMiceBrain.mat';
save(outputbrain, 'modelBrain');

outputbrain = '/Users/eso1993/Desktop/mCADRE_iMiceBrain.xml';
save(model, 'modelBrain');

% export model in xml format
modelBrain = rmfield(modelBrain, 'modelAnnotation');
writeCbModel(modelBrain, 'format', 'sbml', 'filename', '/Users/egabal/Library/CloudStorage/Box-Box/PFOS_Project/rFASTCOROMICS_mouse_models/FASTCORMICS_iMiceBrain.xml');