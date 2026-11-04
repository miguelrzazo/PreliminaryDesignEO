function publish_thesis_figures
% Render the approved MATLAB run; do not repeat or change its calculations.
root = fileparts(fileparts(mfilename('fullpath')));
src = fullfile(root,'review-results');
out = fullfile(root,'Latex_Code','review_figures');
data = fullfile(root,'Compiled data','TutorReview');
if ~exist(out,'dir'), mkdir(out); end
r = jsondecode(fileread(fullfile(src,'summary.json')));
assert(r.selected.Altura_km == 690 && r.selected.Diametro_Pupila_mm == 38);
set(groot,'defaultFigureVisible','off','defaultAxesFontSize',11);

s = readtable(fullfile(src,'Pipeline','OptimumConfigs','reporte_completo_configuraciones.csv'));
v = s(logical(s.Es_Viable),:);
writetable(v,fullfile(data,'viable_designs.csv'));
q = s(s.Num_Satelites==2 & s.Num_Telescopios==2 & s.ID_Detector==2 & strcmp(s.Tipo_Telescopio,'Refractivo'),:);
q = sortrows(q,'Altura_km');
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
    if strcmp(name,'MTF'), file='MTF_Lambda2_Detector2_Telescopio1_resultados.csv'; else, file='SNR_Lambda3_Detector5_Telescopio1_resultados.csv'; end
    t = readtable(fullfile(src,'Pipeline',name,file),'ReadRowNames',true);
    h = str2double(erase(t.Properties.RowNames,'Alt_'));
    diam = str2double(erase(t.Properties.VariableNames,'Diam_'));
    figure('Position',[100 100 850 450]);
    hm = imagesc(diam,h,t{:,:});
    set(hm,'AlphaData',~isnan(t{:,:}));
    set(gca,'Color','w'); axis xy; colorbar; hold on;
    plot(38,690,'rp','MarkerFaceColor','w','MarkerSize',14,'LineWidth',1.5);
    xlabel('Diámetro de pupila [mm]'); ylabel('Altura nominal [km]'); title(name);
    finish(lower(name));
end

figure('Position',[100 100 850 450]);
plot(r.latitudes_deg,r.revisit_days/1.2,'o-'); hold on;
plot(r.latitudes_deg,r.revisit_days,'s-'); yline(7,'k--');
xlabel('Latitud [grados N]'); ylabel('Tiempo de revisita [días]'); grid on;
legend('Geométrica','Efectiva: factor 6/5','Requisito','Location','best'); finish('revisit');
copyfile(fullfile(src,'latitude_revisit.csv'),data);

figure('Position',[100 100 850 550]);
for c = 1:2
    if c==1, name='frozen'; else, name='circular'; end
    t = readtable(fullfile(src,'FrozenSelected',[name '_evolution.csv']));
    copyfile(fullfile(src,'FrozenSelected',[name '_evolution.csv']),data);
    subplot(2,1,1); plot(t.MissionDay/365.25,t.Eccentricity); hold on;
    subplot(2,1,2); plot(t.MissionDay/365.25,max(t.AscendingArcVariation_km,t.DescendingArcVariation_km)); hold on;
end
subplot(2,1,1); ylabel('Excentricidad media'); grid on; legend('Frozen','Inicialmente circular');
subplot(2,1,2); xlabel('Tiempo de misión [años]'); ylabel('Variación por arco [km]'); grid on;
finish('frozen_evolution');
copyfile(fullfile(src,'ContactsSGP4','memoria.jpg'),fullfile(out,'memory.jpg'));
copyfile(fullfile(src,'ContactsSGP4','contacts_sgp4.csv'),data);
copyfile(fullfile(src,'EndOfLife','end_of_life.csv'),data);

% A compact annex keeps all viable rows without retaining old console output.
fid = fopen(fullfile(root,'Latex_Code','zzAnnex','TablasAnnex','review_designs.tex'),'w');
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
