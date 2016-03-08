function sim = displayPoints(aPts, bPts, sim)
    if(nargin < 3)
        sim = 1000;
    end

    set(gca, 'FontSize', 12)

    hold off
    plot(bPts(:,1),bPts(:,2),'go');
    hold on
    plot(aPts(:,1),aPts(:,2),'b.');
    axis equal
    
    if(sim < 1000)
        title(['similarity = ', num2str(sim)]);
    end
    
    bxmin = min(bPts(:,1));
    bxmax = max(bPts(:,1));
    bymin = min(bPts(:,2));
    bymax = max(bPts(:,2));
    
    xlim([bxmin - 1.5, bxmax + 1.5])
    ylim([bymin - 1.5, bymax + 1.5])

    pbaspect([1, 1, 1])
    drawnow

    unused = 1;
end