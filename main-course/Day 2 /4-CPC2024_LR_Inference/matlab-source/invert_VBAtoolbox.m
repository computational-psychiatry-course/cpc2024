function [post, out] = invert_VBAtoolbox(data)
% Bayesian logistic regression using the VBA toolbox (variational laplace)
% -------------------------------------------------------------------------
% This script is a minimal demo showing how to run the inference only given
% the specification of the model prediction, letting the toolbox do the all
% the work.
% -------------------------------------------------------------------------

%% =========================================================================
% Model definition
% =========================================================================
% Note that in the toolbox, parameters of static models are called "phi"

% mapping between input and response
function [gx] = g_logistic(~,param,input,~)
    gx = VBA_sigmoid(param * input);
end

% number of parameters
dim.n_phi = 1;

% indicate we are fitting binary data
options.sources.type = 1;
  
%% =========================================================================
% Inference
% =========================================================================

% specify the prior 
options.priors.muPhi = data.prior.mu;
options.priors.SigmaPhi = data.prior.sigma;

% call the inversion routine
[post, out] = VBA_NLStateSpaceModel (data.y, data.x, [], @g_logistic, dim, options);



end
