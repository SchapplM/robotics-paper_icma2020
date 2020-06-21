function freturn = r04(L1, L2, L3, L4, L5, q1, q2, q3)

	t3383 = sin(q3);
	t3385 = cos(q3);
	t3387 = t3383 .* L4 + t3385 .* L5;
	t3386 = cos(q2);
	t3384 = sin(q2);
	t3382 = t3385 .* L4 - t3383 .* L5 + L3;
	t3381 = t3382 .* t3384 + t3387 .* t3386 + L2;
	freturn = [t3381 .* cos(q1) -t3381 .* sin(q1) t3382 .* t3386 - t3387 .* t3384 + L1];

end
