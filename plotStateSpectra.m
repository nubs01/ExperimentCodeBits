function fh = plotStateSpectra(spectra,freq,states,fh,plotRow,nTets,tetNum,stateOrder)

    figure(fh);
    n1 = numel(stateOrder);
    n2 = (plotRow-1)*n1+1;
    cc = 'bmgcryk';
    for k=1:numel(stateOrder)
        idx = strcmp(states,stateOrder{k});
        S = spectra(idx,:);
        subplot(nTets,numel(stateOrder),n2+k-1)
        hold on
        if isempty(S)
            text(0.01,0.5,sprintf('Error 404: No %s episodes found.',stateOrder{k}),'FontSize',11)
        else
            plot(freq,10*log10(S),'k-','LineWidth',1,'Color',[0 0 0]+0.75);
            mh = plot(freq,10*log10(mean(S,1)),[cc(k) '-'],'LineWidth',3);
        end
        if plotRow == nTets
            xlabel('Frequency (Hz)')
        end
        if k==1
            ylabel({'Mean Power (dB/Hz)',sprintf('Tet #%02i',tetNum)})
        end
        if plotRow == 1
            title(stateOrder{k})
            if ~isempty(S)
                legend(mh,sprintf('n=%i',size(S,1)))
            end
        end
    end
