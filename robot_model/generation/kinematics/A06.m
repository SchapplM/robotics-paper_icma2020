function freturn = A06(q1, q2, q3, q4, q5, q6)

	t3360 = sin(q3);
	t3361 = sin(q2);
	t3365 = cos(q3);
	t3366 = cos(q2);
	t3353 = t3361 .* t3360 - t3366 .* t3365;
	t3359 = sin(q4);
	t3373 = t3353 .* t3359;
	t3358 = sin(q5);
	t3364 = cos(q4);
	t3372 = t3358 .* t3364;
	t3362 = sin(q1);
	t3371 = t3359 .* t3362;
	t3367 = cos(q1);
	t3370 = t3359 .* t3367;
	t3363 = cos(q5);
	t3369 = t3363 .* t3364;
	t3351 = t3358 .* t3365 + t3360 .* t3369;
	t3352 = -t3360 .* t3358 + t3365 .* t3369;
	t3368 = t3351 .* t3366 + t3361 .* t3352;
	t3357 = -q6 + 0.227e3 ./ 0.1200e4 .* q5;
	t3356 = cos(t3357);
	t3355 = sin(t3357);
	t3354 = t3366 .* t3360 + t3361 .* t3365;
	t3350 = t3354 .* t3370 + t3362 .* t3364;
	t3349 = t3354 .* t3371 - t3367 .* t3364;
	t3348 = t3363 .* t3353 + t3354 .* t3372;
	t3347 = t3361 .* t3351 - t3352 .* t3366;
	t3346 = t3363 .* t3371 - t3368 .* t3367;
	t3345 = t3368 .* t3362 + t3363 .* t3370;
	freturn = [t3346 .* t3356 - t3355 .* t3350 t3346 .* t3355 + t3356 .* t3350 -t3348 .* t3367 + t3358 .* t3371; t3345 .* t3356 + t3349 .* t3355 t3345 .* t3355 - t3349 .* t3356 t3348 .* t3362 + t3358 .* t3370; t3347 .* t3356 + t3355 .* t3373 t3347 .* t3355 - t3356 .* t3373 t3353 .* t3372 - t3363 .* t3354;];

end