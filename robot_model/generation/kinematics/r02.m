function freturn = r02(L1, L2, L3, q1, q2)

	t3374 = sin(q2) .* L3 + L2;
	freturn = [cos(q1) .* t3374 -sin(q1) .* t3374 cos(q2) .* L3 + L1];

end
