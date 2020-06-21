
# Kinematic constraints for Cloos QRC 350
# Description
# The robot has a coupling at the hand axes. This worksheet is necessary to fit this into the dynamics toolbox.
# The constraints can also be declared solely in the robot_env. This however leads to some functions not working as expected
# Sources
# [Diekmeyer2018_S678] Diekmeyer, Jonas: Identifikation der inversen Dynamik eines seriellen Roboters im geschlossenen Regelkreis; Studienarbeit; Leibniz Universität Hannover
# 
# Author
# Moritz Schappler, moritz.schappler@imes.uni-hannover.de, 2020-06
# (C) Institut für Mechatronische Systeme, Leibniz Universität Hannover
# Initialisierung
interface(warnlevel=0): # Unterdrücke die folgende Warnung.
restart: # Gibt eine Warnung, wenn über Terminal-Maple mit read gestartet wird.
interface(warnlevel=3):
kin_constraints_exist := true: # Für Speicherung
;
with(StringTools): # Für Zeitausgabe
with(LinearAlgebra):
with(codegen):
with(CodeGeneration):
codegen_act := true:
codegen_opt := 2: # Hoher Optimierungsgrad.
;
read "../helper/proc_MatlabExport":
read "../helper/proc_convert_s_t":
read "../helper/proc_convert_t_s":
with(RealDomain): # Schränkt alle Funktionen auf den reellen Bereich ein. Muss nach Definition von MatlabExport kommen. Sonst geht dieses nicht.
;
read "../robot_codegen_definitions/robot_env":
read sprintf("../codeexport/%s/tmp/tree_floatb_definitions", robot_name):
# Variable mit Winkeln der Zwangsbedingungen nur in Abhängigkeit der verallgemeinerten Koordinaten
kintmp_qs := Matrix(RowDimension(kintmp_s),1):
kintmp_qt := Matrix(RowDimension(kintmp_s),1):
# Ersetzungsausdrücke definieren.
# Variable zum speichern des Sinus und Cosinus der Winkel. Für dieses System ist das eigentlich nicht notwendig. Erstelle Variable, da sie von den anderen Skripten erwartet wird
kintmp_subsexp := Matrix(2*RowDimension(kintmp_s),2):
# Kinematic Constraints
# First axis is defined in negative direction
rho1_qt := -qJ_t(1,1):
# Hand axes 5 and 6 have a coupling via differential gear.
# [Diekmeyer2018_S678], p. 24
rho6_qt := qJ_t(6,1) - kDG*qJ_t(5,1):
kintmp_qt(1,1) := rho1_qt:
kintmp_qt(2,1) := rho6_qt:
# Nachverarbeitung
# Umrechnung in Substitutionsvariablen
kintmp_qs := convert_t_s(kintmp_qt):
# Exportiere Code für folgende Skripte
# Speichere Maple-Ausdruck (Eingabe-Format und internes Format)
save kintmp_qs, kintmp_qt, kin_constraints_exist, kintmp_subsexp, sprintf("../codeexport/%s/tmp/kinematic_constraints_maple_inert", robot_name):
save kintmp_qs, kintmp_qt, kin_constraints_exist, kintmp_subsexp, sprintf("../codeexport/%s/tmp/kinematic_constraints_maple_inert.m", robot_name):
# Liste mit abhängigen konstanten Kinematikparametern erstellen (wichtig für Matlab-Funktionsgenerierung)
read "../helper/proc_list_constant_expressions";
kc_symbols := Matrix(list_constant_expressions( kintmp_qs )):
save kc_symbols, sprintf("../codeexport/%s/tmp/kinematic_constraints_symbols_list_maple", robot_name):
MatlabExport(Transpose(kc_symbols), sprintf("../codeexport/%s/tmp/kinematic_constraints_symbols_list_matlab.m", robot_name), 2):

