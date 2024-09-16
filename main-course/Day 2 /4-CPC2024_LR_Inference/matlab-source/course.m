%% Dice

z = randi(6, 1e7, 1);
z(1:12)'

%% Sum of dice rolls

n_squares = 63;

parfor iGame = 1:1e6
    z = randi (6, n_squares, 1); 
    n_turns_to_win(iGame) = find (cumsum (z) >= n_squares, 1);
end

mean (n_turns_to_win) % e_6(63) = 18.47619

%% simulate data with logistic model
rng(5); 
theta = 1;
data = experiment(theta,true); 

data.prior.mu = 0;
data.prior.sigma = 4;

%% Sampling from a Gaussian distribution
rng(10); 
[samples,lh]=evidence(true,data);

pause()
mu0_hat = mean(samples) 
  
pause()
sigma_0_hat = var(samples) 

pause() 
log_evidence = log(mean(lh))

%% Markov Chain simulation
monte_carlo('MC',true);  

%% Metropolis Hastings
rng(561);   
samples = metropolis_hastings(data);

pause() 
mu_hat = mean(samples)

pause()
sigma_hat = var(samples)

pause()
cred_interval = quantile(samples,[.025 .975])

%% Variational Laplace using the VBA Toolbox
[post, out] = invert_VBAtoolbox(data);

logEvidence = out.F
posterior_mu = post.muPhi
posterior_sigma = post.SigmaPhi
