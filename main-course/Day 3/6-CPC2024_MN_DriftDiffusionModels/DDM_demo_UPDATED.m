%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                      Drift diffusion model
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Goal: Show how drift diffusion model can make predictions about choice
% and RT 

% Additinoal Resources: 
% In the folder that contains this demo you will find several talks given
% by MRN on drift diffusion modeling. There are also recorded lectures
% available on the NEUR 1660 canvas site on DDM modeling. 
% Long Ding has a great tutorial on DDM on her website along with links to
% other papers/resources:
% https://longdecision.github.io/DDM_tutorial/

% Things to learn:
% How DDM selects choice & RT
% How basic parameters of DDM (drift rate/threshold) affect error
% distributions.
% How more complex parameters (across trial variability in start
% point/accumulation rate) affect choice/RT distributions. 

% How to use this script:
% Set trial-to-trial variability parameters to zero. 
% Set nSims to a small number (ie. 5 or 10).
% Then set "doPlots" to true...
% Run script and see what diffusion process looks like. 
% change drift rate/threshold and see how the picture of the diffusion
% process changes...
% Then increase nSims to a big number (>= 1000) and set "doPlots" to false 
% Run simulations and look at Accuracy and RT distributions
% Make predictions for how things will change if you increase drift rate/
% threshold
% Then make changes, test predictions. 
% 

% Try to answer the questions below: 

% 0) What are parameters? How will they affect model?
% 1) Look at while loop -- How many times will it iterate?
% 2) Where is momentary evidence? 
% 3) How is it accumulated?
% 4) When will it stop accumulating?
% 5) Run model, look at output. 
% 6) change drift rate/starting point and see how it affects RT/accuracy. 

% And if you want to get really into things, answer these... but for these
% ones you may need to simulate a *lot* of trials to get stable "answers":

% 7) how do RTs differ for error versus correct trials?
% 8) what happens when we add variability to start point? 
% 9) what happens when we add variability to drift rate? 



% NOTE: Drift rate in this demo is always oriented toward the "correct"
% option rather than to a concrete choice... so if you copy this code to
% make a diffusion model of your particular task, keep in mind that you'll
% probably want to have a drift rate pointed toward a particular option
% (ie. choose left versus choose right) and thus need to compute accuracy a
% bit differently. 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                      Drift Diffusion model!!!!
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Demo parameters:
nSims=1000
doPlots=false;

% Set basic parameters of drift diffusion model:
A     = .005 % Drift rate
y0    =  0   % Starting point
c     = .2   % Noise (standard deviation of momentary evidence)
z     =  10;  % Decision threshold
ndt   = 200;  % Non decision time (in ms)

% Set trial-to-trial variability parameters:
SP_noise=0; % 2
A_var   =0; % 0.01

% Set parameters of simulation 
dt    = 0.02; % how large steps should we simulate in time? 

clear rt
for j =1:nSims % Loop through "trials"
    
    % Initialize things that change within trial (VARIABLES):
    y = y0+normrnd(0,SP_noise); % initialize evidence:
    t=1;   % start on first timestep.

    A_trial=A+normrnd(0, A_var); % set drift rate for this trial
    while abs(y(t)) < z    % Loop until accumulated evidence exceeds threshould (z)
        r=randn(1); % sample a random normal variable for the random walk...
        dW=sqrt(dt).*r; % scale the variance of the variable by the size of the timestep
        dy= A_trial.*dt +c.*dW; % compute change in accumulated evidence (accumulated signal = A.*dt; accumulated noise = c.*dW)
        y(t+1)=y(t)+dy;
        t=t+1;
    end
    

% If we are showing plots, show the diffusion process. 

if doPlots==true
plotLength=800;

true_t=(1:t).*dt;
hold on
plot([0, plotLength], [0, 0], '-k')
plot(true_t, y, 'b')
ylim([-z-1, z+1])

ylabel('Evidence for correct option')
xlabel('Time')
xlim([0, plotLength])

plot([0, plotLength], [z, z], '--k')
plot([0, plotLength], [-z, -z], '--k')
end

rt(j)=t.*dt+ndt;
isAccurate(j)=y(end)>0;

end

%% Show data summary plots:


% Pllot RT separately for correct/incorrect trials:

figure(1)
subplot(2, 1, 1)

hist(rt(isAccurate==1), 50)
xlim([0, 4000])
ylabel('Correct trials')
set(gca, 'fontSize', 24)

subplot(2, 1, 2)

hist(rt(isAccurate==0), 50)
xlim([0, 4000])
ylabel('Error trials')

set(gca, 'fontSize', 24)
xlabel('Reaction time')


% Plot accuracy

figure(2)
hold on
plot([0, 2], [.5, .5], '--k')
bar(mean(isAccurate), 'b')
ylim([0,1])
ylabel('Accuracy')
set(gca, 'fontSize', 24)





