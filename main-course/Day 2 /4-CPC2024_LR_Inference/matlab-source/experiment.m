function m = experiment(theta,display)

m.theta = theta;

m.x = -5:.2:5;
mysig = @(z) 1 ./ (1 + exp(-z));
m.g = mysig(m.x*m.theta);

m.y = +(rand(1,numel(m.x)) < m.g) ;

if display
    figure('Color','w');
    plot(m.x, VBA_sigmoid(m.x*m.theta),'Color',[.5 .5 .5],'LineWidth',2);
    hold on
    plot(m.x,m.y,'k.','MarkerSize',14);
    hold off
    ylim([-.1 1.1])
    set(gca,'XTick',[])
    set(gca,'YTick',[0 1])
    xlabel('number of slides')
    ylabel('sleepiness')
    box off
end







