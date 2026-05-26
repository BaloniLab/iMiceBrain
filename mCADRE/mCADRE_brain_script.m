% initiate cobra toolbox first and get the gurobi license activated 
initCobraToolbox('false');
changeCobraSolver('ibm_cplex', 'LP');
changeCobraSolver('ibm_cplex', 'MILP');

% now getting the model from its path
% working first with the brain data 
% Load the model file
iMM1865model = readCbModel('/Users/egabal/Library/CloudStorage/Box-Box/Esraa_Gabal-PFAS/iMM1865_updated.mat');

% before working with mCADRE get its repo on terminal through 
% type down this command in terminal: 
% ( git clone https://github.com/jaeddy/mcadre )
brain_ubiquity_scores = readtable('/Users/egabal/Library/CloudStorage/Box-Box/PFOS_Project/tidy_data_brain/brain_ubiquity_model_all.csv');

% get the gene list from the expression data [similar ID and order in model]
G = table2cell(brain_ubiquity_scores(:,1));
% the values of gene IDs are encountered as numeric
% convert them to a character array
G = cellfun(@num2str, G, 'UniformOutput', false);

% get the ubiquity score for the chosen quantile (30th)
U_input = table2cell(brain_ubiquity_scores(:,3));
U = cellfun(@double, U_input); 

% Get GPRxn and GPRmat from the model
[GPRrxns, GPRmat] = parse_gprs(iMM1865model);

% create U_GPR to map gene expression scores to the model reactions 
U_GPR = map_gene_scores_to_rxns(iMM1865model, G, U, GPRmat);

% map the high confience of the reactions 
is_C_H = map_high_conf_to_rxns(iMM1865model, GPRmat, GPRrxns, '');

% Calculate expression-based evidence
E_X = calc_expr_evidence(iMM1865model, GPRrxns, U_GPR, is_C_H);


%% calculating the confidence score from the expression data
% 1: Create expressionData structure
expressionData = struct();
% 2: Assign gene cell array to expressionData
expressionData.gene = G;
% 3: assign the ubiquit score values vector to expressionData
% using the ubiquit scores of the 30th quantile
expressionData.value = U;
% Determines the expression data associated to each reaction present in the model
confidenceScores = mapExpressionToReactions(iMM1865model, expressionData);

% option to perform functional check for the nucleotide salvage pathway (1)
% or not (0)
salvageCheck = 1;

% will not detemine the predefined high confidence reactions NOW
C_H_genes = [];

% choosing the fastcc
method = 2;

% Threshold for consistency
epsilon = 1e-6;

%% run mCADRE as the following:
% ensure mcadre to use fastcc
[consistentRxns, modelConsistent] = fastcc(iMM1865model, epsilon);

% we have Cplex solver so use fastcc, it is the value 2
[PM, GM, C, NC, Z, model_C, pruneTime, cRes] = mcadre(iMM1865model, G, U, confidenceScores, salvageCheck, C_H_genes, method);


%% writing the pruned brain draft reconstruction in a matlab file in a specific path
outputfilepath = '/Users/egabal/Desktop/iMiceBrain.mat'; % write the filename in the path
iMiceBrain = writeCbModel(PM, 'mat', outputfilepath);

% get the reactions with zero expression (measured zero in expression and also missing in the expression data)
% convert first from cell array to a table
Z_array = array2table(Z);

% define the path and write the excel file in [xls]
Zpath_brain = '/Users/eso1993/Desktop/missing_reactions_brain_draft.xls';
writetable(Z_array ,Zpath_brain);

% save the noncore reactions of the generic iMM1865 model
NC_array = array2table(NC);
NC_path = '/Users/eso1993/Desktop/mCADRE_output_brain_draft_reconstruction/Noncore_reactions_in_iMM1865.xls';
writetable(NC_array, NC_path);

% save the core reactions of the generic iMM1865 model
C_array = array2table(C);
C_path = '/Users/eso1993/Desktop/mCADRE_output_brain_draft_reconstruction/Core_reactions_in_iMM1865.xls';
writetable(C_array, C_path);

% update the gene list to remove those unused genes from the prunned draft
iMiceBrain = readCbModel('/Users/eso1993/Desktop/iMiceBrain.mat');
[modelNew] = updateGenes(iMiceBrain);

% save the update draft 
outputfilepath = '/Users/eso1993/Desktop/iMiceBrain.mat'; % write the filename in the path
iMiceBrain = writeCbModel(modelNew, 'mat', outputfilepath);

% get the genes and corresponding metabolic subsystems
% add this to the brain-specific draft 
iMiceBrain_updated_genes = readCbModel('/Users/eso1993/Desktop/mCADRE_iMiceBrain.mat');
[GenSubSystem] = findSubSysGen(iMiceBrain_updated_genes);
iMiceBrain_updated_genes.GenSubSystem = [GenSubSystem];

% get the genes and subsystems in a csv file by copying the values from
% iMiceBrain.GenSubSystem and paste to an excel sheet then
% save it as a csv

%% export the model in xml file 
model = readCbModel('/Users/egabal/Library/CloudStorage/Box-Box/PFOS_Project/mCADRE_output_brain_draft_reconstruction/mCADRE_iMiceBrain.mat');
model = rmfield(model, 'modelAnnotation');
writeCbModel(model, 'format', 'sbml', 'filename', '/Users/egabal/Library/CloudStorage/Box-Box/PFOS_Project/mCADRE_output_brain_draft_reconstruction/mCADRE_iMiceBrain.xml');
