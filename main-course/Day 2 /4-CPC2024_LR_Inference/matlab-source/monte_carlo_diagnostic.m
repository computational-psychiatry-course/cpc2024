function theta = monte_carlo_diagnostic(theta,display)

%
N = numel(theta);
BURNIN   = 50;
DECIMATE = 10 ;

%%
if display
    figure('Color','w','ToolBar','none','WindowStyle','docked'); %,'MenuBar','none'

    % plot timeseries
    subplot(3,1,1)
    h=scatter(1:numel(theta),theta,'.k');
    h.CData = zeros(N,3); 
    xlabel('sample number');
    ylabel('sample value');
    xlim([0 numel(theta)]);

    % discard burnin
    pause
    h.CData(1:BURNIN,1) = 1;
    pause

    % compute and plot autocorrelation
    max_lag = 20;
    for lag = 1:max_lag
        c(lag) = corr(theta(1:end-lag+1)',theta(lag:end)');
    end
    subplot(3,1,2)
    plot(1:max_lag,c);
    xlabel('lag')
    ylabel('autocorrelation')
    xlim([1 max_lag])
    box off
    pause

    % show decimation
    h.CData(:,1) = 1;
    h.CData(BURNIN:DECIMATE:end,1) = 0;

    pause

    % plot histogram
    subplot(3,1,3)
    histogram(theta(BURNIN:DECIMATE:end),'Normalization','pdf','EdgeColor','none','FaceColor',.3*[1 1 1])
    box off
    mu = mean(theta(BURNIN:DECIMATE:end));
    sigma = std(theta(BURNIN:DECIMATE:end));
    s=sprintf('mean = %3.2f\nstd = %3.2f',mu, sigma);
    xlabel('theta')
    ylabel('approximate density')

    % plot sufficient statistics
    pause
    text(mu,mean(get(gca,'Ylim')),s,'FontWeight','bold','FontSize',20)

end
%%
theta = theta(BURNIN:DECIMATE:end);