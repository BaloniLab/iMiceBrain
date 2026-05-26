% calculating the reactions scores for reconstruction of tissue-drafts
% using transcriptomics data for CORDA workflow
% will use the TPM normalized data 
tissue_data = readtable('/Users/eso1993/Library/CloudStorage/Box-Box/PFOS_Project/ftINIT_mouse_models/input_data_ftINIT/TPM_normalized_brain_mapped_against_model.csv');
data = table2cell(tissue_data(:, 2:end)); % Exclude the first column (gene IDs)
data = cell2mat(data); % Convert cell array to numeric matrix
data = data +1; % add a pseudo value of 1 to avoind zeros values 
data_mean = mean(data, 2); % Calculate row-wise mean

% read the model
initCobraToolbox('false');
changeCobraSolver('gurobi', 'all');
model = readCbModel('/Users/eso1993/Library/CloudStorage/Box-Box/PFOS_Project/mouse_model_iMM1865/iMM1865_updated.mat');

% map the expression data 
genes = table2cell(tissue_data(:,1));
genes = cellstr(string(genes));
expressionData.gene = genes; % construct the expressionData structure with first gene ids
expressionData.value = data_mean; % assign the mean expression values 
[expressionRxns, parsedGPR, gene_used] = mapExpressionToReactions(model, expressionData);

% set the threshold:
threshold_25 = prctile(data_mean, 25); % 25th percentile
threshold_30 = prctile(data_mean, 30); % 30th percentile % in case of ONLY liver data, it will be set as the lower threshold
threshold_50 = prctile(data_mean, 50); % 30th percentile % in case of ONLY liver data, it will be set as the upper threshold 

% set the CORDA scores as the following 
% Create an empty array of the same size as expressionRxns
corda_scores = nan(size(expressionRxns)); % Initialize all scores as NaN

% Assign scores based on conditions:
corda_scores(isnan(expressionRxns)) = -1; % NaN reactions get a score of -1
corda_scores(expressionRxns == 0) = 0; % Score 0 for reactions with expression data = 0
corda_scores(expressionRxns > 0 & expressionRxns <= threshold_25) = 1; % Score 1 for rxns below 25th percentile
corda_scores(expressionRxns > threshold_25 & expressionRxns <= threshold_30) = 2; % Score 2 for rxns between 25th and 30th percentile
corda_scores(expressionRxns > threshold_30) = 3; % Score 3 for rxns above 30th percentile

%% Ensure BIOMASS_reaction has the highest score
biomass_index = find(strcmp(model.rxns, 'BIOMASS_reaction'));
if ~isempty(biomass_index)
    corda_scores(biomass_index) = 3;
end % Set the BIOMASS reaction to high confidence

% Export CORDA scores along with reaction IDs to a CSV file
reaction_ids = model.rxns; 
corda_scores_table = table(reaction_ids, corda_scores, 'VariableNames', {'Reaction_ID', 'CORDA_Score'});

% Specify the file path
path = '/Users/eso1993/Desktop/brain_rxns_scores.csv';

% Write to CSV
writetable(corda_scores_table, path);