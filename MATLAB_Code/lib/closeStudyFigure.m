function closeStudyFigure(figureHandle)
%CLOSESTUDYFIGURE Release only the figure owned by a failed or finished export.
if isgraphics(figureHandle), close(figureHandle); end
end
