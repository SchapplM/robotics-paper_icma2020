function freturn = r03(L1, L2, L3, L4, q1, q2, q3)

	t3380 = sin(q3) .* L4;
	t3379 = cos(q2);
	t3378 = sin(q2);
	t3376 = cos(q3) .* L4 + L3;
	t3375 = t3376 .* t3378 + t3379 .* t3380 + L2;
	freturn = [cos(q1) .* t3375 -sin(q1) .* t3375 t3376 .* t3379 - t3378 .* t3380 + L1];

end
