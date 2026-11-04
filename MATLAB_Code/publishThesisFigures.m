function publishThesisFigures(src)
% Render figures from one completed analytical run.
root = fileparts(fileparts(mfilename('fullpath')));
if nargin < 1, src = fullfile(root,'Compiled data'); end
out = fullfile(root,'Latex_Code','figures');
data = fullfile(src,'Mission');
if ~exist(out,'dir'), mkdir(out); end
r = jsondecode(fileread(fullfile(src,'Mission','summary.json')));
assert(isfinite(r.selected.Altura_km) && isfinite(r.selected.Diametro_Pupila_mm));
set(groot,'defaultFigureVisible','off','defaultAxesFontSize',11);

s = readtable(fullfile(src,'Pipeline','OptimumConfigs','reporte_completo_configuraciones.csv'));
v = s(logical(s.Es_Viable),:);
writetable(v,fullfile(data,'viable_designs.csv'));
q = s(s.Num_Satelites==r.selected.Num_Satelites & s.Num_Telescopios==r.selected.Num_Telescopios & s.ID_Detector==r.selected.ID_Detector & strcmp(s.Tipo_Telescopio,r.selected.Tipo_Telescopio),:);
q = sortrows(q,'Altura_km');
figure('Position',[100 100 850 450]);
plot(q.Altura_km,q.Masa_Seca_TM_kg,'--'); hold on;
plot(q.Altura_km,q.Masa_Seca_SEOSAT_kg,':');
plot(q.Altura_km,(q.Masa_Seca_TM_kg+q.Masa_Seca_SEOSAT_kg)/2,'-.');
plot(q.Altura_km,q.Masa_Seca_Satelite_kg,'LineWidth',1.5);
xlabel('Altura nominal [km]'); ylabel('Masa seca por satélite [kg]'); grid on;
legend('Thematic Mapper','SEOSAT','Promedio','Dos telescopios','Location','best');
finish('dry_mass');
figure('Position',[100 100 850 450]);
plot(q.Altura_km,q.Masa_Total_Satelite_kg,'--','Color',[.6 .6 .6]); hold on;
q = q(logical(q.Es_Viable),:);
plot(q.Altura_km,q.Masa_Seca_Satelite_kg,'o-');
plot(q.Altura_km,q.Masa_Combustible_kg,'s-');
plot(q.Altura_km,q.Masa_Total_Satelite_kg,'d-');
xlim([min(q.Altura_km)-20,max(q.Altura_km)+20]);
ylim([0,max(q.Masa_Total_Satelite_kg)*1.12]);
xlabel('Altura nominal [km]'); ylabel('Masa por satélite [kg]'); grid on;
legend('Óptica sin filtro regional','Seca viable','Propulsante de mantenimiento','Selección viable','Location','best');
finish('mass');

figure('Position',[100 100 850 450]); hold on;
labels = {};
for d = unique(v.ID_Detector).'
    z = v(v.ID_Detector==d,:);
    plot(z.Altura_km,z.Diametro_Pupila_mm,'o','LineWidth',1.3);
    labels{end+1} = sprintf('Detector %d, configuraciones viables',d); %#ok<AGROW>
end
xlabel('Altura nominal [km]'); ylabel('Diámetro de pupila [mm]'); grid on;
legend(labels,'Location','best'); finish('viable_apertures');

for item = {'MTF','SNR'}
    name = item{1};
    if strcmp(name,'MTF'), band=2; detector=r.selected.ID_Detector; else, band=3; detector=r.selected.ID_Detector+3; end
    file=sprintf('%s_Lambda%d_Detector%d_Telescopio%d_resultados.csv',name,band,detector,r.telescope_index);
    t = readtable(fullfile(src,'Pipeline',name,file),'ReadRowNames',true);
    h = str2double(erase(t.Properties.RowNames,'Alt_'));
    diam = str2double(erase(t.Properties.VariableNames,'Diam_'));
    figure('Position',[100 100 850 450]);
    hm = imagesc(diam,h,t{:,:});
    set(hm,'AlphaData',~isnan(t{:,:}));
    set(gca,'Color','w'); axis xy; colorbar; hold on;
    plot(r.selected.Diametro_Pupila_mm,r.selected.Altura_km,'rp','MarkerFaceColor','w','MarkerSize',14,'LineWidth',1.5);
    xlabel('Diámetro de pupila [mm]'); ylabel('Altura nominal [km]'); title(name);
    finish(lower(name));
end

figure('Position',[100 100 850 450]);
plot(r.latitudes_deg,r.revisit_days*(1-r.cloud_fraction),'o-'); hold on;
plot(r.latitudes_deg,r.revisit_days,'s-'); yline(r.revisit_requirement_days,'k--');
xlabel('Latitud [grados N]'); ylabel('Tiempo de revisita [días]'); grid on;
legend('Geométrica','Efectiva: corrección media por nubes','Requisito','Location','best'); finish('revisit');


profile = readtable(fullfile(src,'Mission','Frozen','acquisition_arcs.csv'));
figure('Position',[100 100 850 450]);
plot(profile.Latitude_deg,profile.AscendingHeight_km,'-','LineWidth',1.5); hold on;
plot(profile.Latitude_deg,profile.DescendingHeight_km,'--','LineWidth',1.5);
xlabel('Latitud geodésica [grados N]'); ylabel('Altura sobre el elipsoide [km]');
legend('Ascendente','Descendente','Location','best'); grid on;
finish('height_arcs');

figure('Position',[100 100 850 550]);
for c = 1:2
    if c==1, name='frozen'; else, name='circular'; end
    t = readtable(fullfile(src,'Mission','Frozen',[name '_evolution.csv']));

    subplot(2,1,1); plot(t.MissionDay/365.25,t.Eccentricity); hold on;
    subplot(2,1,2); plot(t.MissionDay/365.25,max(t.AscendingArcVariation_km,t.DescendingArcVariation_km)); hold on;
end
subplot(2,1,1); ylabel('Excentricidad media'); grid on; legend('Frozen','Inicialmente circular');
subplot(2,1,2); xlabel('Tiempo de misión [años]'); ylabel('Variación por arco [km]'); grid on;
finish('frozen_evolution');
copyfile(fullfile(src,'Mission','Contacts','memoria.jpg'),fullfile(out,'memory.jpg'));
copyfile(fullfile(src,'Mission','EndOfLife','deorbit.jpg'),fullfile(out,'deorbit.jpg'));



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

    function finish(name)
        set(gcf,'Color','w');
        exportgraphics(gcf,fullfile(out,[name '.pdf']),'ContentType','vector');
        close(gcf);
    end
end
