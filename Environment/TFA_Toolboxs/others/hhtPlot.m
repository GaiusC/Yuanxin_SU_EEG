function [this] = hhtPlot(insFreq, insEnergy, tvec, FRange, MinThres, FreqLoc, isNormFreq)
    isYAxis = strcmpi(FreqLoc,'yaxis');

    % patch currently do not support datetime/duration, convert to
    % double
    if(isduration(tvec) || isdatetime(tvec))
        tvec = seconds(tvec);
    end

    frequnitstrs = getfrequnitstrs;
    if isNormFreq
        tvec = 2*pi*tvec; % Convert time axis to samples
        freqlbl = frequnitstrs{1};
        timelbl = getString(message('shared_signalwavelet:hht:general:TimeUnitSamples'));
    else
        freqlbl = frequnitstrs{2};
        timelbl = 'Time (s)';
    end

    if isYAxis
        xlbl = timelbl;
        ylbl = freqlbl;
        xyrange = [0,tvec(end),FRange(1),FRange(2)];
    else
        xlbl = freqlbl;
        ylbl = timelbl;
        xyrange = [FRange(1),FRange(2),0,tvec(end)];
    end

    this.hfig = newplot;
    numIMF = size(insFreq,2);
    for i = 1:numIMF
        % Plot each IMF
        insfi = insFreq(:,i);
        insei = insEnergy(:,i);
        insei((10*log10(insei))<MinThres) = nan;
        if(isYAxis)
            patch([tvec(1);tvec;tvec(end)], [0;insfi;0], [nan;insei;nan], ...
                'EdgeColor','interp','EdgeAlpha','interp',...
                'FaceColor', 'none', 'FaceVertexAlphaData',[nan;insei;nan],...
                'LineWidth', 2, 'FaceAlpha', 'interp');
        else
            patch([0;insfi;0], [tvec(1);tvec;tvec(end)], [nan;insei;nan], ...
                'EdgeColor','interp','EdgeAlpha','interp',...
                'FaceColor', 'none', 'FaceVertexAlphaData',[nan;insei;nan],...
                'LineWidth', 2, 'FaceAlpha', 'interp');
        end
    end
    axis(xyrange);
    xlabel(xlbl);
    ylabel(ylbl);
    title('Ï£¶û²®ÌØÍ¼');
    colormap parula
    colorbar
end