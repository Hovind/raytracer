void
set_colour(float colour[], float r, float g, float b)
{

	colour[0] = r;
	colour[1] = g;
	colour[2] = b;
}

void
scale(float vec[], float s)
{
	size_t i;
	for (i = 0; i < N; ++i)
		vec[i] *= s;
}

float
mixf(float a, float b, float mix)
{
    return b * mix + a * (1.0 - mix);
}

void
mix(float res[], float lhs[], float rhs[], float mix)
{
	size_t i;
	for (i = 0; i < 3; ++i)
	        res[i] = mixf(lhs[i], rhs[i], mix);
}

float
fresnel(float facing_ratio, float a)
{
	return mixf(1.0, pow(1.0 - facing_ratio, 3), a);
}

void
mul_colours(float res[], float lhs[], float rhs[])
{
	size_t i;
	for (i = 0; i < 3; ++i)
		res[i] = lhs[i] * rhs[i];
}

int
is_light(float colour[])
{
	int i;
	float sum = 0;
	for (i = 0; i < 3; ++i)
		sum += colour[i];
	return sum > 0;
}

float
randomf(void)
{
	return rand() / (float) RAND_MAX;
}

float
randomf_in_range(float min, float max)
{
	return randomf() * (max - min) + min;
}

float
min(float lhs, float rhs)
{
	if (lhs < rhs)
		return lhs;
	else
		return  rhs;
}

float
max(float lhs, float rhs)
{
	if (lhs > rhs)
		return lhs;
	else
		return  rhs;
}


