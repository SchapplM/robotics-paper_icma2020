function freturn = Minimalparameter_fcn(pkin, m, mrSges, Ifges, JAges, fc, fv)

%% Initialization

l2 = pkin(2);
l3 = pkin(3);
l4 = pkin(4);
l5 = pkin(5);
l6 = pkin(6);

J2xx = Ifges(2,1); J3xx = Ifges(3,1); J4xx = Ifges(4,1);
J5xx = Ifges(5,1); J6xx = Ifges(6,1);
J1yy = Ifges(1,2); J2yy = Ifges(2,2); J3yy = Ifges(3,2);
J4yy = Ifges(4,2); J5yy = Ifges(5,2); J6yy = Ifges(6,2);
J2zz = Ifges(2,3); J3zz = Ifges(3,3); J4zz = Ifges(4,3);
J5zz = Ifges(5,3); J6zz = Ifges(6,3);
J2xy = Ifges(2,4); J3xy = Ifges(3,4); J4xy = Ifges(4,4);
J5xy = Ifges(5,4); J6xy = Ifges(6,4);
J2xz = Ifges(2,5); J3xz = Ifges(3,5); J4xz = Ifges(4,5);
J5xz = Ifges(5,5); J6xz = Ifges(6,5); 
J2yz = Ifges(2,6); J3yz = Ifges(3,6); J4yz = Ifges(4,6);
J5yz = Ifges(5,6); J6yz = Ifges(6,6); 

JA1 = JAges(1);
JA2 = JAges(2);
JA3 = JAges(3);
JA4 = JAges(4);
JA5 = JAges(5);
JA6 = JAges(6);

rC1x = mrSges(1,1); rC2x = mrSges(2,1); rC3x = mrSges(3,1);
rC4x = mrSges(4,1); rC5x = mrSges(5,1); rC6x = mrSges(6,1);
rC2y = mrSges(2,2); rC3y = mrSges(3,2); rC4y = mrSges(4,2);
rC5y = mrSges(5,2); rC6y = mrSges(6,2);
rC1z = mrSges(1,3); rC2z = mrSges(2,3); rC3z = mrSges(3,3);
rC4z = mrSges(4,3); rC5z = mrSges(5,3); rC6z = mrSges(6,3);

m1 = m(1);
m2 = m(2);
m3 = m(3);
m4 = m(4);
m5 = m(5);
m6 = m(6);

fc1 = fc(1);
fc2 = fc(2);
fc3 = fc(3);
fc4 = fc(4);
fc5 = fc(5);
fc6 = fc(6);

fv1 = fv(1);
fv2 = fv(2);
fv3 = fv(3);
fv4 = fv(4);
fv5 = fv(5);
fv6 = fv(6);
%% Symbolic Code from Maple (RegressorMatrix.m)
t14 = (m4 + m5 + m6);
t12 = (m3 + t14);
t11 = (m2 + t12);
t9 = l3 .^ 2;
t20 = (t11 .* t9);
t8 = (l4 .^ 2);
t19 = (t8 .* m3);
t10 = (l2 .^ 2);
t18 = (-t9 + t10);
t17 = (m3 .* rC3y);
t16 = (m4 .* rC4y);
t15 = 2 .* l5 .* t16 + J4zz;
t13 = J6yy + (l6 + 2 .* rC6z) .* m6 .* l6;
t7 = l5 .^ 2;
freturn = [J1yy + JA1 + J2yy + t18 .* m2 + J3zz + (2 .* l2 .* rC1x + t10) .* m1 + t12 .* (-t8 + t18) J2xx - J2yy + t20 J2xz + (-m2 .* rC2z + t17) .* l3 J2zz + JA2 - t20 t11 .* l3 + m2 .* rC2x t19 + J3xx - J3zz + t14 .* (t7 + t8) + t15 -l4 .* t17 + J3xy J3xz -t19 + J3yy + t14 .* (t7 - t8) + t15 J3yz t12 .* l4 + m3 .* rC3x t14 .* l5 + m3 .* rC3z + t16 JA3 J4xx - J4zz + J5zz J4yy + J5zz JA4 J5xx - J5zz + t13 J5yy + t13 (rC6z + l6) .* m6 + m5 .* rC5z JA5 J6xx - J6yy J6zz JA6 fc1 fc2 fc3 fc4 fc5 fc6 fv1 fv2 fv3 fv4 fv5 fv6];
