function freturn = A02(q1, q2)
  ddq1 = NaN;
	t3302 = cos(q1);
	t3301 = cos(q2);
	t3300 = sin(q1);
	t3299 = sin(q2);
	freturn = [t3302 .* t3299 t3302 .* t3301 t3300; -t3300 .* t3299 -t3300 .* t3301 t3302; t3301 -t3299 zeros(size(ddq1));];

end
