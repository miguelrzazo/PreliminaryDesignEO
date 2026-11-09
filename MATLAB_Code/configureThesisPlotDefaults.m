function cleanup = configureThesisPlotDefaults()
%CONFIGURETHESISPLOTDEFAULTS Apply LaTeX to every MATLAB plot annotation.
% Return a cleanup object so callers restore the user's MATLAB defaults.
root = groot;
properties = {'defaultTextInterpreter', 'defaultAxesTickLabelInterpreter', ...
    'defaultLegendInterpreter'};
previous = cell(size(properties));
for k = 1:numel(properties)
    previous{k} = get(root, properties{k});
end
for k = 1:numel(properties)
    set(root, properties{k}, 'latex');
end
cleanup = onCleanup(@() restoreDefaults(root, properties, previous));
end

function restoreDefaults(root, properties, previous)
for k = 1:numel(properties)
    set(root, properties{k}, previous{k});
end
end
