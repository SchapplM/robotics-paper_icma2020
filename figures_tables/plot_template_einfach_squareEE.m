% Template for figure formatting
% (C) Institut fuer Mechatronische Systeme, Leibniz Universitaet Hannover

% clc;
clear;

% Load Data
addpath(genpath('../'));
l = [640;250;630;196;805;100]/1000; % Kinematikparameter
load('Experimental Validation/Joint Data Process 2/Process2_desired.mat');
EE_soll = posEE(l,deg2rad(deg_q))*1000;

load('Experimental Validation/Joint Data Process 2/Process2_withoutModel.mat');
EE_ohneModell = posEE(l,deg2rad(deg_q))*1000;

load('Experimental Validation/Joint Data Process 2/Process2_processModel.mat');
EE_mitModell = posEE(l,deg2rad(deg_q))*1000;

%% Farben definieren

Farben.rot     =  [255 0 0 ]/255; %rot
Farben.blau    =  [0 0 255 ]/255; %blau
Farben.gruen   =  [0 128 0 ]/255; %grün
Farben.schwarz =  [0 0 0 ]/255; %schwarz
Farben.magenta =  [255 0 255 ]/255; %magenta
Farben.cyan    =  [0 255 255 ]/255; %cyan
Farben.orange  =  [255 180 0 ]/255; %orange
Farben.grau    =  [136 138 142 ]/255; %grau
Farben.hellrot =  [255 150 150 ]/255; %hellrot

Farben.imesblau   = [0 80 155 ]/255; %imesblau
Farben.imesorange = [231 123 41 ]/255; %imesorange
Farben.imesgruen  = [200 211 23 ]/255; %imesgrün

% Schriftart und -größe setzen
schrift.art     = 'Times'; %um die Schriftart zu ändern muss der string des Texts in z.B. \texfsf{...} geschrieben werden
schrift.groesse = 11;

% Figuregröße und Offsets für Ränder bzw. Abstand zwischen den Plots setzen
laengen.breite_fig           = 12;
laengen.hoehe_fig            = 8;

laengen.offset_breite_links  = 2;
laengen.offset_breite_rechts = 0.5;
laengen.offset_oben          = 1.5;
laengen.Abstand_plots        = 0.1;
laengen.hoehe_label          = 1.2;
% Figure aufrufen
f = figure(1);
clf(1);

set(f,'DefaultAxesUnit','centimeters')
set(f,'DefaultAxesFontName',schrift.art)
set(f,'DefaultAxesFontSize',schrift.groesse)
set(f,'DefaultAxesTickLabelInterpreter', 'latex')
set(f,'DefaultLegendInterpreter', 'latex')
set(f,'defaultTextInterpreter','latex')
set(f,'DefaultTextFontSize',schrift.groesse)


f.Units             = 'centimeters';
f.OuterPosition  	= [30 5 laengen.breite_fig+0.4 laengen.hoehe_fig];
f.Color             = [1 1 1];
f.PaperSize         = [laengen.breite_fig laengen.hoehe_fig];
f.PaperPosition     = [0 0 0 0];
f.PaperPositionMode = 'auto';
f.ToolBar           = 'none';
f.MenuBar           = 'none';


%Anzahl der Zeilen
n=1;

%Anzahl labels in x 
m=1;

%Figureposition setzen
laengen.breite_axes    = laengen.breite_fig - (laengen.offset_breite_links+laengen.offset_breite_rechts);
laengen.hoehe_axes     = (laengen.hoehe_fig-laengen.offset_oben-(n-1)*laengen.Abstand_plots-m*laengen.hoehe_label)/n;

positionen(1).pos      = [laengen.offset_breite_links   laengen.hoehe_fig-(laengen.offset_oben + laengen.hoehe_axes)  laengen.breite_axes     laengen.hoehe_axes];

%% Plot 1

achsen(1).a = axes;
achsen(1).a.Position   = positionen(1).pos; 

plots(1).p1            = plot(EE_soll(:,1), EE_soll(:,2));
plots(1).p1.LineStyle  = '-';
plots(1).p1.LineWidth  = 1.5;
plots(1).p1.Color      = Farben.imesblau;
hold on
plots(1).p2            = plot(EE_ohneModell(:,1), EE_ohneModell(:,2));
plots(1).p2.LineStyle  = '-';
plots(1).p2.LineWidth  = 1.5;
plots(1).p2.Color      = Farben.imesorange;

plots(1).p4            = plot(EE_mitModell(:,1), EE_mitModell(:,2));
plots(1).p4.LineStyle  = '-';
plots(1).p4.LineWidth  = 1.5;
plots(1).p4.Color      = Farben.gruen;

% plots(1).p3            = line([0, i],[e_opt.axes(1),e_opt.axes(1)], 'LineStyle', '--');
% plots(1).p3.LineStyle  = '--';
% plots(1).p3.LineWidth  = 1.5;
% plots(1).p3.Color      = Farben.rot;

% plots(1).p5            = plot(0:i-1, eV.axes_rel(1,:));
% plots(1).p5.LineStyle  = '-';
% plots(1).p5.LineWidth  = 1.5;
% plots(1).p5.Color      = Farben.magenta;
% 
% plots(1).p6            = plot(0:i-1, eV.axes_rel(1,:));
% plots(1).p6.LineStyle  = '-';
% plots(1).p6.LineWidth  = 1.5;
% plots(1).p6.Color      = Farben.grau;


legende(1).l             = legend('reference trajectory', 'without model', 'process model', 'Location','southwest', 'Interpreter','latex');
legende(1).l.Units       = 'centimeters';
legende(1).l.Box         = 'on';


limits.ymin            = 235;
limits.ymax            = 250;    
limits.xmin            = 1341;
limits.xmax            = 1343;

achsen(1).a.XLim        = [limits.xmin limits.xmax];
achsen(1).a.YLim        = [limits.ymin limits.ymax];
achsen(1).a.XGrid       = 'on';
achsen(1).a.YGrid       = 'on';

label(1).y              = ylabel('y-axis (mm)', 'Interpreter','latex');
label(1).x              = xlabel('x-axis (mm)', 'Interpreter','latex');



% legenden(1).l            = legend([plots(1).p1,plots(1).p2],'Modell 1','Modell 2');
% legenden(1).l.Units      = 'centimeters';
% legenden(1).l.Position   = [10.2 1.8 0 0];
% 
% legenden(1).l.Box        = 'off';

f.Renderer = 'painters';
Dateiname = 'EE_square';
saveas(f,strcat(Dateiname,'.emf'))
saveas(f,strcat(Dateiname,'.pdf'))

%% Speichern und exportieren
% filename = 'mean_e_axes_Präsi_PnP';
% set(gcf,'PaperPositionMode','auto');
% print('-dmeta',[filename,'_tmp.emf']); %Ausgabe als .eps
% system(['gswin64c.exe -dNOPAUSE -dBATCH -dNOCACHE -dEPSCrop -sDEVICE=epswrite -sOutputFile=',filename,'.eps ',filename,'_tmp.eps']);
% system(['gswin64c.exe -dNOPAUSE -dBATCH -dNOCACHE -dEPSCrop -sDEVICE=pdfwrite -sOutputFile=',filename,'.pdf ',filename,'_tmp.eps']);
% system(['gswin64c.exe -dNOPAUSE -dBATCH -dNOCACHE -r300 -dEPSCrop -sDEVICE=png16m   -sOutputFile=',filename,'.png ',filename,'_tmp.eps']);

rmpath(genpath('./'));
