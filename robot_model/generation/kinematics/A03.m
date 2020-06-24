function freturn = A03(q1, q2, q3)
  ddq1 = NaN;
	t3311 = cos(q2);
	t3310 = cos(q3);
	t3309 = cos(q1);
	t3308 = sin(q1);
	t3307 = sin(q2);
	t3306 = sin(q3);
	t3304 = t3311 .* t3306 + t3307 .* t3310;
	t3303 = t3307 .* t3306 - t3311 .* t3310;
	freturn = [t3309 .* t3304 -t3308 -t3309 .* t3303; -t3308 .* t3304 -t3309 t3308 .* t3303; -t3303 zeros(size(ddq1)) -t3304;];

end
