function results = runWeatherLinkStudy(config)
%RUNWEATHERLINKSTUDY Preserve the historical seasonal probabilities and rate ranges.
% Two interleaved uniform draws per day reproduce the original seeded sequence.
previousPath = path; pathCleanup = onCleanup(@() path(previousPath));
addpath(fullfile(fileparts(mfilename('fullpath')),'lib'));
validateattributes(config.dayCount,{'numeric'},{'scalar','integer','finite','>=',7});
validateattributes(config.seasonProbabilities_percent,{'numeric'},{'real','finite','nonnegative','<=',100});
validateattributes(config.rateRanges_Mbps,{'numeric'},{'size',[4 2],'real','finite','nonnegative'});
assert(all(config.rateRanges_Mbps(:,2)>=config.rateRanges_Mbps(:,1)),'eo:WeatherRates','Rate ranges must have increasing endpoints.');
assert(isequal(size(config.seasonProbabilities_percent),[4 4]) && ...
    all(abs(sum(config.seasonProbabilities_percent,2)-100)<1e-12), ...
    'eo:WeatherProbabilities','Each seasonal probability row must sum to 100.');
previousRng = rng; rngCleanup = onCleanup(@() rng(previousRng));
rng(config.randomSeed);
draws = rand(2,config.dayCount);
days = (1:config.dayCount).'; seasons = 4*ones(config.dayCount,1);
seasons(days>=80 & days<=172)=1; seasons(days>=173 & days<=265)=2;
seasons(days>=266 & days<=355)=3;
probabilities = config.seasonProbabilities_percent(seasons,:);
conditionIndices = 1 + sum(draws(1,:).'*100 > cumsum(probabilities(:,1:3),2),2);
ranges = config.rateRanges_Mbps(conditionIndices,:);
dailyRates_Mbps = ranges(:,1) + draws(2,:).' .* (ranges(:,2)-ranges(:,1));
weekCount = floor(config.dayCount/7);
weeklyRates_Mbps = mean(reshape(dailyRates_Mbps(1:weekCount*7),7,weekCount),1);
results.config=config; results.conditionIndices=conditionIndices;
results.dailyConditions=reshape(config.conditionNames(conditionIndices),[],1);
results.dailyRates_Mbps=dailyRates_Mbps; results.weekIndices=1:weekCount;
results.weeklyRates_Mbps=weeklyRates_Mbps;
% Preserve the old annual statistic: mean of complete weekly means (364 days).
results.annualWeeklyMean_Mbps=mean(weeklyRates_Mbps);
results.annualDailyMean_Mbps=mean(dailyRates_Mbps);
if ~config.writeOutputs, return; end
if ~exist(config.outputDirectory,'dir'), mkdir(config.outputDirectory); end
graphicsCleanup = studyGraphicsDefaults();
weekIndices=results.weekIndices; weeklyRates_Mbps=results.weeklyRates_Mbps;
annualWeeklyMean_Mbps=results.annualWeeklyMean_Mbps; conditionNames=config.conditionNames;
dailyConditions=results.dailyConditions; dailyRates_Mbps=results.dailyRates_Mbps;
figure('Position', [100, 100, 1000, 600]);
hold on;
primaryColor = [0, 0.4470, 0.7410];
meanColor = [0.8500, 0.3250, 0.0980];
plot(weekIndices, weeklyRates_Mbps, 'o-', 'Color', primaryColor, ...
    'LineWidth', 2.5, 'MarkerSize', 8, 'MarkerFaceColor', primaryColor);
yline(annualWeeklyMean_Mbps, '--', 'Color', meanColor, 'LineWidth', 3, ...
    'DisplayName', sprintf('Media anual (%.1f Mbps)', annualWeeklyMean_Mbps));
xlim([1 max(2,numel(weekIndices))]);
xlabel('Semana del a\~no', 'FontSize', 12);
ylabel('Velocidad promedio (Mbps)', 'FontSize', 12);
title('Velocidad promedio semanal en banda X durante un a\~no', 'FontSize', 14);
grid on; grid minor;
legend({'Velocidad promedio semanal', sprintf('Media anual (%.1f Mbps)', annualWeeklyMean_Mbps)}, ...
    'Location', 'best', 'FontSize', 11);
hold off;
outputDirectory = config.outputDirectory;
exportgraphics(gcf, fullfile(outputDirectory, 'velmediasemanal.png'), 'Resolution', 150);
figure('Position', [150, 150, 1000, 600]);
hold on;
conditionColors = [
    0, 0.4470, 0.7410;
    0.8500, 0.3250, 0.0980;
    0.6350, 0.0780, 0.1840;
    0.4660, 0.6740, 0.1880
];
colorMap = containers.Map(conditionNames, ...
    {conditionColors(1,:), conditionColors(2,:), ...
     conditionColors(3,:), conditionColors(4,:)});
for conditionIndex = 1:length(conditionNames)
    conditionName = conditionNames{conditionIndex};
    conditionMask = strcmp(dailyConditions, conditionName);
    conditionDays = find(conditionMask);
    conditionRates_Mbps = dailyRates_Mbps(conditionMask);
    scatter(conditionDays, conditionRates_Mbps, 50, ...
        'MarkerFaceColor', colorMap(conditionName), ...
        'MarkerEdgeColor', 'none', 'MarkerFaceAlpha', 0.8, ...
        'DisplayName', conditionName);
end
annualDailyMean_Mbps = mean(dailyRates_Mbps);
yline(annualDailyMean_Mbps, '--', 'Color', 'k', 'LineWidth', 3, ...
    'DisplayName', sprintf('Media anual (%.1f Mbps)', annualDailyMean_Mbps));
xlim([1 max(2,config.dayCount)]);
xlabel('D\''ia del a\~no', 'FontSize', 12);
ylabel('Velocidad de descarga (Mbps)', 'FontSize', 12);
title('Velocidad diaria en banda X durante un a\~no', 'FontSize', 14);
grid on; grid minor;
legend('Location', 'best', 'FontSize', 11);
hold off;
exportgraphics(gcf, fullfile(outputDirectory, 'veldiariabandax.png'), 'Resolution', 150);
end
