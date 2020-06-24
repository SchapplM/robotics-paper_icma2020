% Return reduced dynamics parameters of the robot
% CloosQRC350DE
%
% The dynamics parameters can be modified by user input
% into the Maple code generation.
%
% Input: Full arrays of dynamics parameters
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[L1,L2,L3,L4,L5,L6,kDG]';
% m [7x1]
%   mass of all robot links (including the base)
% rSges [7x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% Icges [7x6]
%   inertia of all robot links about their respective center of mass, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertiavector2matrix.m)
%
% Output: Modified dynamics parameters, consistent with dynamics
% m,rSges,Icges
% (same meaning as input)

% Quelle: HybrDyn-Toolbox
% Datum: 2020-06-23 21:15
% Revision: 3f22bf868ffa24e21e77a0fe3b46e78b2d6fdc1f (2020-06-23)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function [m,rSges,Icges] = CloosQRC350DE_dynamics_parameters_modification(pkin,m,rSges,Icges)

%% Coder Information
%#codegen
%$cgargs {zeros(7,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'CloosQRC350DE_dynamics_parameters_modification: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'CloosQRC350DE_dynamics_parameters_modification: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'CloosQRC350DE_dynamics_parameters_modification: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'CloosQRC350DE_dynamics_parameters_modification: Icges has to be [7x6] (double)'); 

%% Variable Initialization
% Complete set of dynamics parameters that can be reduced by the user.
L1 = pkin(1);
L2 = pkin(2);
L3 = pkin(3);
L4 = pkin(4);
L5 = pkin(5);
L6 = pkin(6);
kDG = pkin(7);

M0 = m(1);
M1 = m(2);
M2 = m(3);
M3 = m(4);
M4 = m(5);
M5 = m(6);
M6 = m(7);

SX0 = rSges(1,1);
SY0 = rSges(1,2);
SZ0 = rSges(1,3);
SX1 = rSges(2,1);
SY1 = rSges(2,2);
SZ1 = rSges(2,3);
SX2 = rSges(3,1);
SY2 = rSges(3,2);
SZ2 = rSges(3,3);
SX3 = rSges(4,1);
SY3 = rSges(4,2);
SZ3 = rSges(4,3);
SX4 = rSges(5,1);
SY4 = rSges(5,2);
SZ4 = rSges(5,3);
SX5 = rSges(6,1);
SY5 = rSges(6,2);
SZ5 = rSges(6,3);
SX6 = rSges(7,1);
SY6 = rSges(7,2);
SZ6 = rSges(7,3);

XXC0 = Icges(1,1);
XYC0 = Icges(1,4);
XZC0 = Icges(1,5);
YYC0 = Icges(1,2);
YZC0 = Icges(1,6);
ZZC0 = Icges(1,3);
XXC1 = Icges(2,1);
XYC1 = Icges(2,4);
XZC1 = Icges(2,5);
YYC1 = Icges(2,2);
YZC1 = Icges(2,6);
ZZC1 = Icges(2,3);
XXC2 = Icges(3,1);
XYC2 = Icges(3,4);
XZC2 = Icges(3,5);
YYC2 = Icges(3,2);
YZC2 = Icges(3,6);
ZZC2 = Icges(3,3);
XXC3 = Icges(4,1);
XYC3 = Icges(4,4);
XZC3 = Icges(4,5);
YYC3 = Icges(4,2);
YZC3 = Icges(4,6);
ZZC3 = Icges(4,3);
XXC4 = Icges(5,1);
XYC4 = Icges(5,4);
XZC4 = Icges(5,5);
YYC4 = Icges(5,2);
YZC4 = Icges(5,6);
ZZC4 = Icges(5,3);
XXC5 = Icges(6,1);
XYC5 = Icges(6,4);
XZC5 = Icges(6,5);
YYC5 = Icges(6,2);
YZC5 = Icges(6,6);
ZZC5 = Icges(6,3);
XXC6 = Icges(7,1);
XYC6 = Icges(7,4);
XZC6 = Icges(7,5);
YYC6 = Icges(7,2);
YZC6 = Icges(7,6);
ZZC6 = Icges(7,3);

%% Parameter Postprocessing / Set Output
% Create the reduced set of dynamics parameters that is also used for generation of the dynamics equations.
% Aus parameters_dyn_mges_matlab.m
t1 = [M0; M1; M2; M3; M4; M5; M6;];
m = t1;

% Aus parameters_dyn_rSges_matlab.m
t1 = [0, 0, 0; SX1, SY1, SZ1; SX2, 0, SZ2; SX3, SY3, SZ3; 0, 0, SZ4; 0, SY5, 0; 0, 0, SZ6;];
rSges = t1;

% Aus parameters_dyn_Icges_matlab.m
t1 = [0, 0, 0, 0, 0, 0; XXC1, YYC1, ZZC1, XYC1, XZC1, YZC1; XXC2, YYC2, ZZC2, 0, XZC2, 0; XXC3, YYC3, ZZC3, XYC3, XZC3, YZC3; XXC4, YYC4, ZZC4, 0, 0, 0; XXC5, YYC5, ZZC5, 0, 0, 0; XXC6, YYC6, ZZC6, 0, 0, 0;];
Icges = t1;
