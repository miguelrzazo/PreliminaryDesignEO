function publishThesisFigures(src)
%PUBLISHTHESISFIGURES Update the original thesis figures from current results.
root = fileparts(fileparts(mfilename('fullpath')));
if nargin < 1, src = fullfile(root,'Compiled data'); end
previousPath=path; pathCleanup=onCleanup(@() path(previousPath));
addpath(fullfile(root,'MATLAB_Code','lib'));
data = fullfile(src,'Mission');
missionFigures = fullfile(root,'Latex_Code','5.Mission');
payloadFigures = fullfile(root,'Latex_Code','4.Payload');
r = jsondecode(fileread(fullfile(data,'summary.json')));
plotDefaultsCleanup = configureThesisPlotDefaults();
set(groot,'defaultFigureVisible','off','defaultAxesFontSize',10);
s = readtable(fullfile(src,'Pipeline','OptimumConfigs','reporte_completo_configuraciones.csv'));
v = s(logical(s.Es_Viable),:);
writetable(v,fullfile(data,'viable_designs.csv'));
q = s(s.Num_Satelites==r.selected.Num_Satelites & s.Num_Telescopios==r.selected.Num_Telescopios & s.ID_Detector==r.selected.ID_Detector & strcmp(s.Tipo_Telescopio,r.selected.Tipo_Telescopio),:);
previousDirectory=pwd; directoryCleanup=onCleanup(@() cd(previousDirectory));
cd(fullfile(src,'Pipeline')); writeMassPlots(q); cd(previousDirectory);
prefix=sprintf('%dsat_%dtel_Det%d_%s',r.selected.Num_Satelites,r.selected.Num_Telescopios,r.selected.ID_Detector,r.selected.Tipo_Telescopio);
copyfile(fullfile(src,'Pipeline','Masa_seca',['MasaSeca_Comparativa_' prefix '.jpg']),missionFigures,'f');
copyfile(fullfile(src,'Pipeline','Masa_total',['MasaDesglose_' prefix '.pdf']),missionFigures,'f');
for item = {'MTF','SNR'}
    name = item{1};
    if strcmp(name,'MTF'), band=2; detector=r.selected.ID_Detector; threshold=0.25; else, band=3; detector=r.selected.ID_Detector+3; threshold=400; end
    file=sprintf('%s_Lambda%d_Detector%d_Telescopio%d_resultados.csv',name,band,detector,r.telescope_index);
    t=readtable(fullfile(src,'Pipeline',name,file),'ReadRowNames',true);
    h=str2double(erase(t.Properties.RowNames,'Alt_'));
    d=str2double(erase(t.Properties.VariableNames,'Diam_'));
    figure('Position',[100 100 900 700]);
    imageHandle=imagesc(d,h,t{:,:}); set(imageHandle,'AlphaData',~isnan(t{:,:}));
    set(gca,'Color','w','YDir','normal'); colormap(parula);
    if strcmp(name,'MTF'), wavelength=2.01; else, wavelength=0.76; end
    if strcmp(name,'MTF'), clim([threshold 0.35]); else, clim([threshold 2000]); end
    cb=colorbar; cb.Label.String=name; cb.Label.Interpreter='latex'; cb.TickLabelInterpreter='latex';
    title(sprintf('%s - Refractivo - Detector %d\n($\\lambda$ = %.2f $\\mu$m, Blanco: %s $< %g$)',name,r.selected.ID_Detector,wavelength,name,threshold),'Interpreter','latex','FontSize',14);
    xlabel('Di\''ametro de Pupila (mm)','Interpreter','latex'); ylabel('Altura Orbital (km)','Interpreter','latex');
    markDesign(r.selected.Diametro_Pupila_mm,r.selected.Altura_km);
    saveFigure(fullfile(missionFigures,['design_point_' lower(name) '.pdf']));
end

% Retain the original coverage heatmap, with the current design point.
file=sprintf('%d Satelite(s); %d Telescopio(s): %s; Detector %d',r.selected.Num_Satelites,r.selected.Num_Telescopios,r.selected.Tipo_Telescopio,r.selected.ID_Detector);
coverage=readmatrix(fullfile(src,'Pipeline','coverage',['coverage_' file '.csv']));
h=(200:10:1000).'; swath=20:10:500;
assert(isequal(size(coverage),[numel(h),numel(swath)]),'eo:CoverageAxes','Unexpected coverage axes.');
figure('Position',[100 100 1000 600]);
hm=imagesc(swath,h,coverage); set(hm,'AlphaData',isfinite(coverage));
set(gca,'Color','w','YDir','normal'); colormap(parula); clim([1 r.revisit_requirement_days]);
cb=colorbar; cb.Label.String='Tiempo revisita (dias)'; cb.Label.Interpreter='latex'; hold on;
plot(2*h*tand(10),h,'r--','LineWidth',1.5,'DisplayName','Limite FoV');
for detectorCount=1:3
    xline(detectorCount*2048*r.selected.GSD/1000,'k--',sprintf('%d Det',detectorCount),'HandleVisibility','off');
end
xlabel('Swath (km)','Interpreter','latex'); ylabel('Altura (km)','Interpreter','latex');
title(sprintf('Tiempo Revisita - %s\n(Blanco: $> %g$ dias)',file,r.revisit_requirement_days),'Interpreter','latex');
xlim([min(swath) max(swath)]); ylim([min(h) max(h)]); legend('Location','southeast','Interpreter','latex');
markDesign(r.selected.Swath_Requerido_km,r.selected.Altura_km);
saveFigure(fullfile(missionFigures,'design_point_coverage.pdf'));

% Use the original two h-D layouts, with the regional viability filter.
writeRegionalSelection(v,payloadFigures,r.selected.Num_Satelites,r.selected.Num_Telescopios);
profile=readtable(fullfile(data,'Frozen','acquisition_arcs.csv'));
[~,~,~,circular]=computeAcquisitionHeightEnvelope(0,pi/2,profile.Latitude_deg, ...
    r.selected.Inclinacion_SSO_deg,6378.137e3+r.selected.Altura_km*1000,6378.137e3);
figure('Position',[100 100 1000 600],'Theme','light'); hold on;
plot(profile.Latitude_deg,profile.AscendingHeight_km,'b-','LineWidth',1.5);
plot(profile.Latitude_deg,profile.DescendingHeight_km,'b--','LineWidth',1.5);
plot(profile.Latitude_deg,circular(1:height(profile)),'r-','LineWidth',1.5,'HandleVisibility','off');
plot(profile.Latitude_deg,fliplr(circular(height(profile)+1:end)),'r--','LineWidth',1.5,'HandleVisibility','off');
yline(r.selected.Altura_km,'k:'); grid on;
xlabel('Latitud del arco de adquisicion [grados]','Interpreter','latex'); ylabel('Altura sobre el elipsoide [km]','Interpreter','latex');
title('Variacion de altura sobre la zona de adquisicion','Interpreter','latex');
legend('Frozen, ascendente','Frozen, descendente','Altura nominal','Location','eastoutside','Interpreter','latex');
saveFigure(fullfile(missionFigures,'frozen_coverage_height_analysis.pdf'));
copyfile(fullfile(data,'Contacts','memoria.jpg'),fullfile(root,'Latex_Code','7.Segmento_Tierra','memoria.jpg'),'f');
copyfile(fullfile(data,'EndOfLife','deorbit.jpg'),missionFigures,'f');

% Existing chapter figures use copies of the same pipeline exports.
copyfile(fullfile(src,'Pipeline','coverage','*.pdf'),fullfile(root,'Latex_Code','4.Payload','Coverage'),'f');
copyfile(fullfile(src,'Pipeline','MTF','*.pdf'),fullfile(root,'Latex_Code','4.Payload','MTF'),'f');
copyfile(fullfile(src,'Pipeline','SNR','*.pdf'),fullfile(root,'Latex_Code','4.Payload','SNR'),'f');
% A compact annex keeps all viable rows without retaining old console output.
fid = fopen(fullfile(root,'Latex_Code','zzAnnex','TablasAnnex','viable_designs.tex'),'w');
fprintf(fid,'\\begin{longtable}{rrrlrrrr}\n\\caption{Soluciones viables tras el filtro regional. Masa de selección sin reserva de fin de misión.}\\\\\n\\toprule\nSat. & Tel. & Det. & Óptica & $h$ [km] & $D$ [mm] & $T_{eff}$ [d] & $m_{const}$ [kg] \\\\\n\\midrule\n\\endfirsthead\n\\toprule\nSat. & Tel. & Det. & Óptica & $h$ [km] & $D$ [mm] & $T_{eff}$ [d] & $m_{const}$ [kg] \\\\\n\\midrule\n\\endhead\n');
for k=1:height(v)
    fprintf(fid,'%d & %d & %d & %s & %.0f & %.0f & %.3f & %.3f \\\\\n',v.Num_Satelites(k),v.Num_Telescopios(k),v.ID_Detector(k),v.Tipo_Telescopio{k},v.Altura_km(k),v.Diametro_Pupila_mm(k),v.Dias_Cobertura(k),v.Masa_Total_Constelacion_kg(k));
end
fprintf(fid,'\\bottomrule\n\\end{longtable}\n'); fclose(fid);
fprintf('THESIS_FIGURES_COMPLETE: %d viable designs\n',height(v));
end

function markDesign(x,y)
hold on; plot(x,y,'ro','MarkerFaceColor','r','MarkerSize',6,'HandleVisibility','off');
text(x+8,y+60,'Punto de dise\~no','Interpreter','latex','FontSize',16, ...
    'BackgroundColor','white','EdgeColor','black');
plot([x+8,x],[y+50,y],'k-','LineWidth',1.5,'HandleVisibility','off');
end
function saveFigure(file)
set(gcf,'Color','w'); exportgraphics(gcf,file,'ContentType','vector'); close(gcf);
fprintf('Updated %s\n',file);
end
function writeRegionalSelection(v,out,satellites,telescopes)
v=v(v.Num_Satelites==satellites & v.Num_Telescopios==telescopes,:);
names={'Refractivo','Korsch','Cassegrain','TMA'};
figure('Position',[100 100 900 700]);
for detector=1:3
    subplot(3,1,detector); plotDetector(v,detector,names);
    title(sprintf('$h$ vs $D_{min}$ viable - Detector %d - Configuracion: %d satelite(s), %d telescopio(s) por satelite',detector,satellites,telescopes),'Interpreter','latex','FontSize',10);
end
saveFigure(fullfile(out,sprintf('HvsDmin_Viable_%dsat_%dtel.pdf',satellites,telescopes)));
figure('Position',[100 100 700 500]); hold on;
styles={'r--','g:','b-.'};
for detector=1:3
    q=sortrows(v(v.ID_Detector==detector & strcmp(v.Tipo_Telescopio,'Refractivo'),:),'Altura_km');
    plot(q.Altura_km,q.Diametro_Pupila_mm,styles{detector},'LineWidth',2,'DisplayName',sprintf('Detector %d',detector));
end
title(sprintf('$h$ vs $D_{min}$ - Refractivo - Configuracion: %d satelite(s), %d telescopio(s) por satelite',satellites,telescopes),'Interpreter','latex','FontSize',10);
xlabel('Altura orbital (km)','Interpreter','latex'); ylabel('$D_{min}$ Refractivo (mm)','Interpreter','latex');
grid on; legend('Location','best','Interpreter','latex'); axis tight;
saveFigure(fullfile(out,sprintf('HvsDmin_Viable_Refractivo_Config_%dsat_%dtel.pdf',satellites,telescopes)));
end
function plotDetector(v,detector,names)
hold on;
for telescope=1:numel(names)
    q=sortrows(v(v.ID_Detector==detector & strcmp(v.Tipo_Telescopio,names{telescope}),:),'Altura_km');
    if isempty(q), continue; end
    plot(q.Altura_km,q.Diametro_Pupila_mm,'LineWidth',2,'DisplayName',names{telescope});
end
xlabel('Altura orbital (km)','Interpreter','latex','FontSize',9); ylabel('$D_{min}$ (mm)','Interpreter','latex','FontSize',9);
legend('Location','northwest','Interpreter','latex'); grid on; axis tight;
end
