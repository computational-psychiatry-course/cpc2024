function [theta, ll]=evidence(display,m)

%% initialization
%  ---------------------------------------------------------------------
% number of samples

N = 700; 

% number examples
Nex = 3;
Slim = 5;


%% compute un-normalized log-posterior

mysig = @(z) 1 ./ (1 + exp(-z));

function lup=log_likelihood(m,th)
    sg = mysig(m.x*th);
    lup = sum( ... % aggregate over obsevations
            m.y.*log(sg) + (1-m.y).*log(1-sg) ... %binomial log-likelihood
        );
end

function lup=log_prior(m,th)
    lup = -.5*(th-m.prior.mu).^2/m.prior.sigma;
end


% init plot

[MLE, max_LL] = fminsearch(@(x) - log_likelihood(m, x), 0);
maxL = exp(-max_LL);

%%
dt = .01;
theta_prop = (-Slim:dt:Slim); % XXX

for i=1 : numel(theta_prop)
    lpd(i)=log_prior(m,theta_prop(i));
end
lpd = exp(lpd);
lpd = 5*lpd/nansum(lpd*dt);
    
for i=1 : numel(theta_prop)
 lld(i)=log_likelihood(m,theta_prop(i));
end
lld = exp(lld);
lld = lld/nansum(lld*dt);

maxL = max(lld);
lpd = lpd / maxL;
lld = lld / maxL;

max_proba = max(max(lpd),max(lld));

if display
    [handles]=init_display();
end


plot(handles.plot(2),theta_prop,lpd,'Color',[.184 .333 .592], 'LineWidth',2);
plot(handles.plot(2),theta_prop,lld,'Color',[.70 .085 0], 'LineWidth',2);



handles.plot(2).YLim = [-eps max_proba];

handles.likelihood = plot(handles.plot(2),[NaN NaN],[0 0],'--','Color',[.70 .085 0], 'LineWidth',1);


%% sampling loop
%  ---------------------------------------------------------------------


% starting value
theta = [] ;
ll = [] ;

% iterate
for t=1:N 
   
   
  
   proposal =  sqrt(m.prior.sigma)*randn(1);
            
   % store new sample
   theta(t) = proposal;
   ll(t) = log_likelihood(m, proposal);
   
   display1();
   display2(); 
   
end

ll = exp(ll);

%% subfunctions
%  ================================================================================
    function [handles]=init_display()
            handles.plot=newFigure();
            handles.proposal = [NaN NaN NaN NaN];
            handles.hX = NaN;
            
            
            %
            handles.histB = histogram(handles.plot(1),[],-Slim:.1:Slim,'Normalization','pdf','EdgeColor','none','FaceColor',[.1 .1 .5]);
            handles.plot(1).Visible = 'off'; 
            handles.plot(1).YLim = [0 1];
            handles.plot(1).XLim = [-Slim Slim];
            
            %
            handles.histR=histogram(handles.plot(4), [] ,linspace(0, 1,25),'Normalization','pdf','EdgeColor','none','FaceColor',[.5 .1 .1]);
            handles.histR.Orientation = 'horizontal';
            handles.plot(4).Visible = 'off'; 
            handles.plot(4).XLim = [0 1];
            handles.plot(4).YLim = [0 1];
                        
            % trace
            handles.trace = plot(handles.plot(3),NaN,NaN,'.','Color',.3*[1 1 1]);
            handles.trace.MarkerSize = 15;
            hold on
            handles.hX = plot(handles.plot(3),NaN,NaN,'ko');
            hold off
            handles.plot(3).YLim = [-50 0];
            handles.plot(3).XLim = [-Slim Slim];
            %box off ;
            handles.plot(3).YTick = [];
            handles.plot(3).YColor = 'w';

            % data
            plot(handles.plot(5),m.x, m.y, 'k.','MarkerSize',13);
            hold on
            xx = linspace(min(m.x), max(m.x), 100);
            gx = sigmoid(xx,m.theta);
            plot(handles.plot(5), xx, gx, 'Color',[.5 .5 .5], 'LineWidth',2);
            handles.plot(5).YTick = [];
            handles.plot(5).YColor = 'w';
            handles.plot(5).Box = 'off';
            
            handles.fit = plot(handles.plot(5), xx, nan*gx,'Color',[.70 .085 0], 'LineWidth',2);
            handles.plot(5).YLim = [-.03 1];

    end
    
    function display1()
        % plot histroy of values
        if t<=Nex
            pause
        end
        plot_trace();
        if t<=Nex
            pause
        end
        plot_fit();
       
    end



    function display2()
        
        % plot visited density
        if t<=Nex
            pause
        end
        plot_histogram() ;
        
        
        drawnow
    end

    function h=newFigure()
        figure('Color','w','ToolBar','none','WindowStyle','docked'); %,'MenuBar','none'
        % bottom histogram
        h(1)=subplot('Position',[.05 .05 .45 .13]) ; box off ; 
        h(1).Visible = 'off'; 
        % distribution
        h(2)=subplot('Position',[.05 .8 .45 .18]); box off ; hold on
        h(2).Visible = 'off';
        % trace
        h(3)=subplot('Position',[.05 .2 .45 .6],'NextPlot','add'); box off ;
        h(3).YColor = 'w';
        % right histogram
        h(4)=subplot('Position',[.51 .8 .09 .18]); box off ;
        h(4).Visible = 'off';
        
        % fit 
        h(5)=subplot('Position',[.70 .8 .26 .18]); box off ;
        h(5).Visible = 'off';
    end


    function plot_fit()
        xx = linspace(min(m.x), max(m.x), 100);
        handles.fit.YData = sigmoid(xx, theta(end));
        
        set(handles.likelihood, 'XData', theta(end)*[1 1], 'YData', [0 interp1(theta_prop,lld, theta(end))]);
        
    end

    function plot_trace()


        set(handles.trace, ...
            'XData', theta, ...
            'YData', 1:numel(theta) ...
            )
        handles.plot(3).YLim = [min(numel(theta)-50,0) numel(theta)+1];
        
        set(handles.hX, ...
            'XData', theta(end), ...
            'YData', numel(theta) ...
            )
    end

    function plot_histogram()
      handles.histB.Data = theta;
      handles.plot(1).YLim = [0 max(handles.histB.Values)];

      handles.histR.Data = exp(ll+max_LL);
      handles.plot(4).XLim = [0 max(handles.histR.Values)]; 
      handles.plot(4).YLim = handles.plot(2).YLim;
    end
end
    