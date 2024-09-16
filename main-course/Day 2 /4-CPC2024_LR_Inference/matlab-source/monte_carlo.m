function theta=monte_carlo(mode,display,m)

%% initialization
%  ---------------------------------------------------------------------
% number of samples
if display
    N = 3e2; 
else
    N = 3e6;
end

% proposal distribution
proposal_sigma = .3;

% number examples
Nex = 5;

% init plot
if display
    [handles,Slim]=init_display();
end

%% compute un-normalized log-posterior

mysig = @(z) 1 ./ (1 + exp(-z));

function lup=log_posterior(m,th)
    lup = -.5*(th-m.prior.mu).^2/m.prior.sigma ... % Gaussian prior
        + sum( ... % aggregate over obsevations
            m.y.*log(mysig(m.x*th)) + (1-m.y).*log(1-mysig(m.x*th)) ... %binomial log-likelihood
        );
end

%% sampling loop
%  ---------------------------------------------------------------------

% starting value
theta = 0 ;

% iterate
for t=1:N-1  
   
   if display, display1(); end
   
   % draw new theta
   switch mode
       % classic sampling
       case 'SP'
            proposal =  proposal_sigma*randn(1);
            accept=true;
       % markov chain
       case 'MC'
            proposal = theta(t) + proposal_sigma*randn(1);
            accept=true;
       % Metropolis Hasting
       case 'MH'
            % propose new sample
            proposal = theta(t) + proposal_sigma*randn(1)
            
            % compute un-normalized posterior
            old = log_posterior(m,theta(t));
            new = log_posterior(m,proposal);
            
            % do we get warmer?
            accept_prob = min(1, exp(new-old));
            accept =  accept_prob > rand(1) ;
            % if yes, accept, otherwise, stay in place
            if ~accept
                proposal = theta(t); 
            end        
   end
   
   % store new sample
   theta(t+1) = proposal;
   
   if display, display2(); end
   
end

% clean up (burnin, decimation)
if strcmp(mode,'MH')
    pause();
    theta=monte_carlo_diagnostic(theta,display);
end

%% subfunctions
%  ================================================================================
    function [handles,Slim]=init_display()
            handles.plot=newFigure();
            handles.proposal = [NaN NaN NaN];
            handles.hX = NaN;
            if mode=='MH'
                Slim = 3;
            else
                Slim=1;
            end
    end
    
    function display1()
        % plot histroy of values
        plot_trace();
        if t<=Nex
            pause
        end
        
        % plot proposal ditribution
        plot_proposal() ;
        if t<=Nex
            pause
        end
    end

    function display2()
        % plot proposition
        plot_proposed_x() ;
        if t<=Nex
            pause
        end
        % show acceptance
        if strcmp(mode,'MH')
                plot_log_posterior()
            if t==1
                pause
            end
            %fprintf(' move by %f; a = %04.3f: %s\n', proposal-theta(t),accept_prob,elvis(accept,'accept','reject'));
            if accept
                set(handles.proposal(3),'MarkerFaceColor','g');
            else
                set(handles.proposal(3),'MarkerFaceColor','r');
            end
            if t<=Nex
                pause
            end
        end 
        % plot visited density
        plot_histogram() ;   
        drawnow
    end

    function h=newFigure()
        figure('Color','w','ToolBar','none','WindowStyle','docked'); %,'MenuBar','none'
        h(1)=subplot(6,1,6)  ; box off ; 
        h(1).Visible = 'off'; 
        h(2)=subplot(6,1,1); box off ; hold on
        h(2).Visible = 'off';
        h(3)=subplot(6,1,2:5); box off ;
        h(3).YColor = 'w';
     end

    function plot_proposal()
        try, delete(handles.proposal); end

        switch mode
            case 'SP'
                theta0=0;
            otherwise
                theta0 = theta(end);
        end
        theta_prop = (-Slim:.01:Slim);
        
        proppdf = normpdf(theta_prop,theta0,proposal_sigma);
        handles.proposal(1)=plot(handles.plot(2),theta_prop,proppdf,'k', 'LineWidth',2);
        hold on
        handles.proposal(2)=plot(handles.plot(2),[theta0 theta0],[0 normpdf(theta0,theta0,proposal_sigma)],':k', 'LineWidth',1);
        hold off
        handles.plot(2).Visible = 'off';
        set(handles.plot(2),'XLim',[-Slim Slim]);

    end

    function plot_log_posterior()
        subplot(handles.plot(4))
        npost = 100;      
        theta_post=linspace(-Slim,Slim,npost);
        for i=1:npost
            lup(i)=log_posterior(m,theta_post(i));
        end
        lup = exp(lup);
        plot(theta_post,lup,'m', 'LineWidth',1);
        handles.plot(4).XColor = 'k';
        ylabel('p(\theta) p(y|\theta)')
        %set(handles.plot(4),'YTick',[])
        ylim([min(lup) max(lup)])
        box off
        xlim([-Slim Slim]);
    end

    function plot_proposed_x()
        try, delete(handles.proposal(3)); end
        handles.proposal(3)=scatter(handles.plot(2),theta(end),0,100,[.93 .70 .12],'filled');
        while abs(theta(end)) > Slim
            Slim = Slim+1;
        end
        set(handles.plot(2),'XLim',[-Slim Slim]);
    end

    function plot_trace()
         try, delete(handles.hX); end 
         subplot(handles.plot(3));
         plot(theta(2:end),(N-numel(theta)+2):N,'.','Color',.3*[1 1 1], 'LineWidth',1,'MarkerSize',15);
         hold on
         handles.hX=scatter(theta(end),N,100,'k','filled');
         hold off
         ylim([N-max(numel(theta),50) N]);
         xlim([-Slim Slim]);
         box off ;
         set(gca,'YTick',[]);
         handles.plot(3).YColor = 'w';
    end

    function plot_histogram()
        f=histogram(handles.plot(1),theta(2:end),-Slim:.05:Slim,'Normalization','pdf','EdgeColor','none','FaceColor',.3*[1 1 1]);
        handles.plot(1).Visible = 'off'; 
        set(handles.plot(1),'YLim',[0 max(f.Values)]);
        set(handles.plot(1),'XLim',[-Slim Slim]);
    end
end
    