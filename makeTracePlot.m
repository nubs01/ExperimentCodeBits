function fh = makeTracePlot(X,Y,SD,timeBreaks,timeTicks,timeTickLabels,titleStr,xStr,yStr,plotOpt)
    traceAxSize = [0.04 0.15 0.92 0.65];
    traceFigSize = [67 434 1735 505];

    fh = figure('Position',traceFigSize,'Visible','off');
    shadedErrorPlot(X,Y,SD,plotOpt)
    for tb = 1:size(timeBreaks,1)
        wid = timeBreaks(tb,2)-timeBreaks(tb,1);
        ylim = get(gca,'ylim');
        hei = diff(ylim);
        rectangle('Position',[timeBreaks(tb,1) ylim(1) wid hei],'FaceColor',[1 1 1],'EdgeColor','none')
    end
    set(gca,'XTick',timeTicks,'XTickLabels',timeTickLabels)
    xlim([X(1) X(end)])
    title(titleStr)
    ylabel(yStr)
    xlabel(xStr)
    set(gca,'Position',traceAxSize)
    set(gca,'TickLength',[.005 .1])
