% modifiying some fields in iMiceBrain-APOE 
initCobraToolbox('false');
changeCobraSolver('gurobi', 'all');

% reading the model
model = readCbModel('/Users/eso1993/Library/CloudStorage/Box-Box/PFOS_Project/drafts_statistics_AND_consensus_draft/models/final_models/iMiceBrain.mat');

% now assign the sbo terms 
%model = annotateSBOTerms(model);
%model.metSBOTerms = model.metSBOTerms'; % transpose to vertical order 
%model.geneSBOTerms = model.geneSBOTerms';

% change the metabolites format in model so as the demand/sink reactions
model.mets = regexprep(model.mets, '\[([a-zA-Z0-9_]+)\]', '($1)');
sink_rxns = startsWith(model.rxns, 'SK_');

model.rxns(sink_rxns) = strrep(model.rxns(sink_rxns), 'SK_', 'sink_');
model.rxns = regexprep(model.rxns, '_([a-z]{1,2})$', '($1)');
model.rxns = regexprep(model.rxns, '__([LDR])', '_$1');

for i = 1:length(model.mets)
    met = model.mets{i};

    % Skip if metabolite already has -L or -D enantiomer naming
    if contains(met, '-L') || contains(met, '-D') || contains(met, '-B')
        continue;
    end

    % Otherwise, replace all underscores with hyphens
    model.mets{i} = strrep(met, '_', '-');
end
model.mets = regexprep(model.mets, '--', '-');
model.mets = regexprep(model.mets, '-hs', '_hs');
rxn_idx = strcmp(model.rxns, 'EX_4hpro_LT(e)');
model.rxns(rxn_idx) = {'EX_4HPRO'};

% to avoid the mistaken dupliacted reaction in sanity check 
%rxn_idx = strcmp(model.rxns, 'DM_Rtotal3[e]');
%model.rxns(rxn_idx) = {'DM_Rtotal3_ALT[e]'};

% remove it temprarily 
model = removeRxns(model, 'DM_Rtotal3_ALT[e]');

% check for the NA in rxnGeneMatrix 
%model = buildRxnGeneMat (model);
%emptyCols = find(all(model.rxnGeneMat == 0, 1));  % Find columns with all zeros

% theer are genes with no GPR rules 
%if isempty(emptyCols)
    %disp('? No empty columns in rxnGeneMat.');
%else
    %fprintf('? %d empty columns found in rxnGeneMat.\n', numel(emptyCols));
    %disp('Corresponding genes:');
    %disp(model.genes(emptyCols(:)));  % Show the gene IDs with no mapped reactions
%end

% remove them and update the model
%model2 = updateGenes(model2);

% adjustments for the cocnistsency checkup 
epsilon = 1e-6;
[~, fluxConsistentRxns] = fastcc(model, epsilon);

% Get reaction IDs that were skipped in consistency testing (-1)
inconsistent_rxnIDs = model.rxns(~fluxConsistentRxns);

% check the reactions which were skipped in tetsing (-1 value)
rxns_minus1 = model.rxns(fluxConsistentRxns == -1);
rxn_idx = find(fluxConsistentRxns == -1);
bounds = [model.lb(rxn_idx), model.ub(rxn_idx)];
T = table(rxns_minus1, bounds(:,1), bounds(:,2), 'VariableNames', {'Reaction', 'LB', 'UB'});

% just temprarily open these reactions for checking flux consistency again
model.lb(rxn_idx) = -1000;
model.ub(rxn_idx) = 1000;
[~, fluxConsistentRxns] = fastcc(model, epsilon); % now all metabolic reactions are flux consistent

% perform sanity check 
if ~isfield(model, 'rev')
    model.rev = model.lb < 0 & model.ub > 0;
end

resultsFileName = 'iMiceLiver_SanityCheck_Results.xlsx';  
resultsPath = '/Users/eso1993/Desktop/';               

ExtraCellCompIn  = 'e';
ExtraCellCompOut = 'e';
runSingleGeneDeletion = 1;

param = struct();
param.testFluxConsistency = true;
param.testRev = 1;
param.checkDuplicates = 1;

[TableChecks, Table_csources, CSourcesTestedRxns, TestSolutionNameOpenSinks, TestSolutionNameClosedSinks] = ...
    performSanityChecksonReconMouse(model, resultsFileName, ExtraCellCompIn, ExtraCellCompOut, ...
    runSingleGeneDeletion, resultsPath, param);

% check for the duplicated reaction
method = 'FR';           % Method: flux reaction comparison
removeFlag = 0;          % Just detect, don?t remove

% Run duplicate check
[~, removedRxnInd, keptRxnInd] = checkDuplicateRxn(model, method, removeFlag);
duplicatedRxnID = model.rxns{6033};
disp(['Duplicated Reaction ID: ', duplicatedRxnID]);

% investigate the ATP production from different resources 
[Table_csources,TestedRxns,PercTestedRxns] = testATPYieldFromCsources(model, 'iMiceBrain');

%% detect deadend mets 
dead = detectDeadEnds(model); % no deadend mets

%% working on other fields within the model
model = readCbModel('/Users/eso1993/Library/CloudStorage/Box-Box/PFOS_Project/drafts_statistics_AND_consensus_draft/models/final_models/iMiceBrain.mat');

% inspect that all model genes are involved in GPR rules 
all_model_genes = string(model.genes);

% Combine all grRules into one big string
all_grRules = join(string(model.grRules), ';');

% Find which genes do NOT appear in any grRule
genes_not_in_grRules = all_model_genes(~contains(all_grRules, all_model_genes));

% Display count and preview
disp("Genes in model.genes but NOT in any grRule:");
disp(genes_not_in_grRules);

%% detect the non-brain genes in model and associated gpr 
no_brain_genes = readtable('/Users/eso1993/Desktop/genes_not_in_brain_mouse_rnaseq_genes.csv');
no_brain_genes = string(no_brain_genes{:,1});

% Initialize categories
genes_in_or = [];
genes_in_and = [];
genes_in_single = [];

% Check each gene
for g = 1:length(no_brain_genes)
    gene = no_brain_genes(g);
    if ~any(model.genes == gene)
        continue
    end

    for i = 1:length(model.grRules)
        gpr = model.grRules{i};
        if isempty(gpr)
            continue
        end

        if contains(gpr, gene)
            if contains(gpr, ' or ')
                genes_in_or(end+1,1) = gene;
            elseif contains(gpr, ' and ')
                genes_in_and(end+1,1) = gene;
            else
                genes_in_single(end+1,1) = gene;
            end
            break
        end
    end
end

% Remove duplicates
genes_in_or = unique(genes_in_or);
genes_in_and = unique(genes_in_and);
genes_in_single = unique(genes_in_single);

% Summary
fprintf('\nSummary:\n');
fprintf('Genes in OR rules: %d\n', numel(genes_in_or));
fprintf('Genes in AND rules: %d\n', numel(genes_in_and));
fprintf('Genes in single-gene rules: %d\n', numel(genes_in_single));

%% now remove only the genes under OR rule and update model accordingly 
% Combine the two sets
genes_to_remove = unique(genes_in_or);
genes_to_remove = cellstr(string(genes_to_remove));

% Remove from model
model_new = removeGenesFromModel(model, genes_to_remove);

%% inspect if these genes were also removed from GPR
% Initialize tracker
still_present_genes = {};

for i = 1:length(genes_to_remove)
    gene = genes_to_remove{i};
    
    % Check if the gene still appears in any grRule
    if any(contains(model_new.grRules, gene))
        still_present_genes{end+1,1} = gene;
    end
end

% Display results
if isempty(still_present_genes)
    disp('? All specified genes have been removed from GPRs.');
else
    fprintf('?? %d gene(s) still appear in GPRs:\n', numel(still_present_genes));
    disp(still_present_genes);
end

% All specified genes have been removed from GPRs.

%% export the model
% export model and check fva 
file = '/Users/eso1993/Desktop/iMiceBrain.mat';
save(file, 'model_new');

writeCbModel(model_new, 'format', 'sbml', 'fileName', '/Users/eso1993/Desktop/iMiceBrain.xml');

%% performing single gene knocout simulation with both FBA and iMOMA 
test_gene = readtable('/Users/egabal/Library/CloudStorage/Box-Box/Papers preparations/CP2 paper/supplmentary files/gene_knockout_candidates_iMiceBrain.xlsx');
lethal_genes = test_gene{:,1};            
lethal_genes = lethal_genes(~isnan(lethal_genes));
lethal_genes = string(lethal_genes);

% FBA method
[grRatioDble, grRateKO, grRateWT]  = doubleGeneDeletion(model, 'FBA', model.genes, model.genes, 'False');
imagesc(grRatio)

% IMoma 
[grRatio, grRateKO, grRateWT, hasEffect, delRxns, fluxSolution] = singleGeneDeletion(model, 'lMOMA', lethal_genes, 'False');
imagesc(grRatio)

% viable genes
viable_genes = test_gene{:,2};            
viable_genes = viable_genes(~isnan(viable_genes));
viable_genes = string(viable_genes);

% FBA method
[grRatio, grRateKO, grRateWT, hasEffect, delRxns, fluxSolution] = singleGeneDeletion(model, 'FBA', viable_genes, 'False');
imagesc(grRatio)

% lMOMA method
[grRatio, grRateKO, grRateWT, hasEffect, delRxns, fluxSolution] = singleGeneDeletion(model, 'lMOMA', viable_genes, 'False');
imagesc(grRatio)



