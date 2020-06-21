function freturn = A01(q1)
  ddq1 = NaN;
	t3298 = cos(q1);
	t3297 = sin(q1);
	freturn = [t3298 zeros(size(ddq1)) t3297; -t3297 zeros(size(ddq1)) t3298; zeros(size(ddq1)) -1 zeros(size(ddq1));];

end
