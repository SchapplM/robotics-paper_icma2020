function freturn = r05(L1, L2, L3, L4, L5, q1, q2, q3)

	t3390 = sin(q3);
	t3392 = cos(q3);
	t3394 = t3390 .* L4 + t3392 .* L5;
	t3393 = cos(q2);
	t3391 = sin(q2);
	t3389 = t3392 .* L4 - t3390 .* L5 + L3;
	t3388 = t3389 .* t3391 + t3394 .* t3393 + L2;
	freturn = [t3388 .* cos(q1) -t3388 .* sin(q1) t3389 .* t3393 - t3394 .* t3391 + L1];

end
