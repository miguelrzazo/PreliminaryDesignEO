function cleanup = studyGraphicsDefaults()
%STUDYGRAPHICSDEFAULTS Scope the original LaTeX defaults to one study call.
textInterpreter = get(groot,'defaultTextInterpreter');
tickInterpreter = get(groot,'defaultAxesTickLabelInterpreter');
legendInterpreter = get(groot,'defaultLegendInterpreter');
cleanup = onCleanup(@() set(groot,'defaultTextInterpreter',textInterpreter, ...
    'defaultAxesTickLabelInterpreter',tickInterpreter,'defaultLegendInterpreter',legendInterpreter));
set(groot,'defaultTextInterpreter','latex','defaultAxesTickLabelInterpreter','latex', ...
    'defaultLegendInterpreter','latex');
end
