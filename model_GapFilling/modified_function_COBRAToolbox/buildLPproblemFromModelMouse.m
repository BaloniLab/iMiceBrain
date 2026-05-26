function optProblem = buildLPproblemFromModel(model, verify)
% Wrapper function to replace deprecated buildLPproblemFromModel with buildOptProblemFromModel
% as followed from the github (https://github.com/opencobra/cobratoolbox/blob/master/deprecated/buildLPproblemFromModel.m)

warning('Function buildLPproblemFromModel is deprecated. Using buildOptProblemFromModel instead.');

optProblem = buildOptProblemFromModel(model, verify);

end
