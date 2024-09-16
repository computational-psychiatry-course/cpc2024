function theta = metropolis_hastings(data)


N = 1e3;

% proposal distribution
proposal_sigma = .3;

% number examples
Nex = 4;

% init plot
display = init_display();
display.lims = [-1.5 3.5];

%% sampling loop
%  ---------------------------------------------------------------------

% starting value
theta = 0 ;
ljoint_current = log_posterior(data,theta);
% iterate
for t = 1 : N  
   
    % propose new sample
    proposal = theta(t) + proposal_sigma * randn();

    % compute un-normalized posterior
    ljoint_proposal = log_posterior(data,proposal);

    % do we get warmer?
    accept_prob = min(1, exp(ljoint_proposal-ljoint_current));
    accept =  accept_prob > rand(1) ;
    
    
    if t == 1
        display = plot_chain(display);
        display = plot_kernel(display,proposal_sigma);
        plot_joint (display, data);
        display = plot_histogram(display);
    end
    if t <= Nex
       pause
    end
    show_proposal(display, proposal);
    if t <= Nex
       pause
    end
    check_proposal(display, accept);
    
    
    % if yes, accept, otherwise, stay in place
    if accept
        theta(t+1) = proposal;
        ljoint_current = ljoint_proposal;
    else
        theta(t+1) = theta(t); 
    end        

   
  drawnow
  if t <= Nex
       pause
  end
  show_jump(display, theta(end));
  if t <= Nex
       pause
  end
       show_trace(display, theta);
       %update_histogram(theta);
       
       drawnow
   
end

end


function display = init_display()
    % figure
    % ---------------------------------------------------------------------
    display.fig = new_figure();

    % subplots
    % ---------------------------------------------------------------------
    % joint
    display.plot.joint = subplot(7,1,1); 
    % kernel
    display.plot.kernel = subplot(7,1,2);
    display.plot.kernel.Box = 'off';
    display.plot.kernel.NextPlot = 'add';
    % chain
    display.plot.chain = subplot(7,1,3:6);
    display.plot.chain.NextPlot = 'add';
    display.plot.chain.Box = 'off';
    display.plot.chain.YColor = 'w';
    % histogram
    display.plot.hist = subplot(7,1,7);
    display.plot.hist.Box = 'off';
    
end

function fig = new_figure ()
    fig = figure('Color','w','ToolBar','none','WindowStyle','docked'); %,'MenuBar','none'
end

function plot_joint (display, data)
    % compute log_joint
    npost = 1e3;      
    thetas = linspace(-20,20,npost);
    for i=1:npost
        lup(i)=log_posterior(data,thetas(i));
    end
    lup = exp(lup);
    % plot
    plot(display.plot.joint,thetas,lup,'Color',[.77 .16 .63], 'LineWidth',2);
    display.plot.joint.YLabel.String = 'p(\theta) p(y|\theta)';
    display.plot.joint.YTick = [];
    display.plot.joint.YLim = [min(lup) max(lup)];
    display.plot.joint.XLim = display.lims;
    display.plot.joint.Box = 'off';
    
end

function display = plot_kernel(display,proposal_sigma)
        theta_prop = linspace(-20, 20 , 1e3);        
        proppdf = normpdf(theta_prop,0,proposal_sigma);
        display.kernel.density = plot(display.plot.kernel,theta_prop,proppdf,'k', 'LineWidth',2);
        display.kernel.mean = plot(display.plot.kernel,[0 0],[0 normpdf(0,0,proposal_sigma)],':k', 'LineWidth',1);
        display.plot.kernel.Visible = 'off';
        display.plot.kernel.XLim = display.lims;
end

function display = plot_chain(display)

	display.current = scatter(display.plot.chain,0,0,100,'k','filled');   
    
	display.proposal = scatter(display.plot.chain,0,0,100,'filled');
    display.proposal.Visible = 'off';
    
    display.chain = plot(display.plot.chain,NaN,NaN,'.','Color',.5*[1 1 1], 'LineWidth',2,'MarkerSize',15);
    
    display.plot.chain.Visible = 'off';
    display.plot.chain.XLim = display.lims;
    display.plot.chain.YLim = [-50 1];
end

function display = plot_histogram(display)
    edges = display.lims(1) : .1 : display.lims(2);
    display.histogram = histogram(display.plot.hist,[],edges,'Normalization','pdf','EdgeColor','none','FaceColor',[.3 .3 .3]);
    display.plot.hist.Visible = 'off'; 
    display.plot.hist.YLim = [0 1];
    display.plot.hist.XLim = display.lims;
end

function show_proposal(display, theta)
    display.proposal.XData = theta;
    display.proposal.MarkerFaceColor = [.93 .70 .12];
    display.proposal.Visible = 'on';
end

function check_proposal(display, accept)
    if accept
        display.proposal.MarkerFaceColor = 'g';
    else
        display.proposal.MarkerFaceColor = 'r';
    end
end

function show_jump(display, theta)
    display.proposal.Visible = 'off';
    display.current.XData = theta;
    
    display.kernel.mean.XData = [theta theta];
    display.kernel.density.XData = display.kernel.density.XData + theta - median(display.kernel.density.XData); 
end

function show_trace(display, theta)
    set(display.chain, ...
            'XData', theta, ...
            'YData', 1:numel(theta) ...
            );
        
    display.current.YData = numel(theta);
    display.proposal.YData = numel(theta);
    display.plot.chain.YLim = [min(numel(theta)-50,0) numel(theta)+1];
    
    display.histogram.Data = theta;
    display.plot.hist.YLim = [0 max(display.histogram.Values)];
        
end

function s = sigmoid (z)
    s = 1 ./ (1 + exp(-z));
end

function lup = log_posterior (data,theta)
    lup = -0.5 * (theta-data.prior.mu).^2/data.prior.sigma ... % Gaussian prior
        + sum( ... % aggregate over obsevations
            data.y.*log(sigmoid(data.x*theta)) + (1-data.y).*log(1-sigmoid(data.x*theta)) ... %binomial log-likelihood
        );
end