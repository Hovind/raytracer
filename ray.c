#include <stdio.h>
#include <string.h>

#define N 3

struct ray {
	float origin[N];
	float dir[N];
	float length;
}

void
set_vec(float vec[], float x, float y, float z)
{
	vec[0] = x;
	vec[1] = y;
	vec[2] = z;
}

void
copy_vec(float to[], float src[])
{
        memcpy(to, src, N*sizeof(src[0]));
}

void
print_vec(float src[])
{
        size_t i;
        for (i = 0; i < N; ++i)
                printf("%f ", src[i]);
}

void
add_vec(float res[], float lhs[], float rhs[])
{
	size_t i;
	for (i = 0; i < N; ++i)
		res[i] = lhs[i] + rhs[i];
}


void
sub_vec(float res[], float lhs[], float rhs[])
{
	size_t i;
	for (i = 0; i < N; ++i)
		res[i] = lhs[i] - rhs[i];
}

float
dot(float lhs[], float rhs[])
{
	size_t i;
	float ret = 0;
	for (i = 0; i < N; ++i)
		ret += lhs[i] * rhs[i];
	return ret;
}

unsigned char
dotb(float lhs[], float rhs[])
{
	return 255.0 * dot(lhs, rhs);
}

float
length_squared(float vec[])
{
	return dot_vec(vec, vec);
}

float
length(float vec[])
{
	return sqrt(length_squared(vec));
}

void
scale(float vec[], float s)
{
	size_t i;
	for (i = 0; i < N; ++i)
		vec[i] *= s;
}

void
normalize(float vec[])
{
	scale(vec, 1.0 / length(vec));
}

void
construct(float res[], float origin[], float dir[], float length)
{
	float tmp[N];
	copy(tmp, dir);
	scale(tmp, length);
	add(res, origin, tmp);
}

void
construct_refraction_dir(float refraction_dir[], float ray_dir[], float hit_normal[], int inside)
{
	float cosi;
	float k;
	float eta = 1.1;

	if (!inside)
		eta = 1 / eta;

	cosi = -1.0 * dot(hit_normal, ray_dir);
	k = 1 - eta * eta * (1 - cosi * cosi);

	add_weighted(refraction_dir, ray_dir, eta, hit_normal, eta * cosi - sqrt(k));
	normalize(refraction_dir);
}

float
ray_dot(struct ray *lhs, struct ray *rhs)
{
	return lhs->length * rhs->length * dot(lhs->dir, rhs->dir);
}