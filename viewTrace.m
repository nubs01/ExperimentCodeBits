function h = viewTrace(x,t,winSize)
h = figure;
fPos = get(gcf,'Position');
nTrace = size(x,2);
for l=1:nTrace
    subplot(nTrace,1,l)
    plot(t,x(:,l))
    xlim([t(1) t(1)+winSize])
end
setappdata(h,'winNum',1)
backButton = uicontrol(h,'Style','pushbutton','String','<<','Position',[0.01 0.01 50 25],'Callback',['for l=1:' num2str(nTrace) ';subplot(' num2str(nTrace) ',1,l);set(gca,''XLim'',get(gca,''XLim'')-' num2str(winSize) ');end;'],'Units','normalized');
nextButton = uicontrol(h,'Style','pushbutton','String','>>','Position',[fPos(3)-50 0.01 50 25],'Callback',['for l=1:' num2str(nTrace) ';subplot(' num2str(nTrace) ',1,l);set(gca,''XLim'',get(gca,''XLim'')+' num2str(winSize) ');end;']','Units','normalized');
h = {h,backButton,nextButton};

