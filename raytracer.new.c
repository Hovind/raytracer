#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <mpi.h>

#define N 3
#define MAX_RAY_DEPTH 5

#define SET_COLOUR(colour, r, g, b) { colour[0] = r; colour[1] = g; colour[2] = b; }

void
copy(float to[], float src[])
{
        memcpy(to, src, N*sizeof(src[0]));
}

void
print(float src[])
{
        unsigned int i;
        for (i = 0; i < N; ++i)
                printf("%f ", src[i]);
        printf("\n");
}

void
add(float res[], float lhs[], float rhs[])
{
	int i;
	for (i = 0; i < N; ++i)
		res[i] = lhs[i] + rhs[i];
}


void
sub(float res[], float lhs[], float rhs[])
{
	int i;
	for (i = 0; i < N; ++i)
		res[i] = lhs[i] - rhs[i];
}

float
dot(float lhs[], float rhs[])
{
	int i;
	float ret;
	for (i = 0; i < N; ++i)
		ret += lhs[i] * rhs[i];
	return ret;
}

float
length_squared(float vec[])
{
	return dot(vec, vec);
}


float
length(float vec[])
{
	return sqrt(length_squared(vec));
}

void
scale(float vec[], float s)
{
		int i;
		for (i = 0; i < N; ++i)
				vec[i] *= s;
}

void
normalize(float vec[])
{
		scale(vec, length(vec));
}

float
mixf(float a, float b, float mix)
{
    return b * mix + a * (1.0 - mix);
}

void
mix(float res[], float lhs[], float rhs[], float mix)
{
		int i;
    for (i = 0; i < N; ++i)
        res[i] = mixf(lhs[i], rhs[i], mix);
}

float
fresnel(float hit_normal[], float ray_dir[], float a)
{
	return mixf(pow(1.0 + dot(ray_dir, hit_normal), 3), 1.0, a);
}

void
add_weighted(float res[], float lhs[], float a, float rhs[], float b)
{
	int i;
	for (i = 0; i < N; ++i)
		res[i] = a * lhs[i] + b * rhs[i];
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

struct sphere {
	float center[N];
	float radius2;
	float surface_colour[3];
	float emission_colour[3];
	float transparency;
	float reflection;
};
	
int
intersect(float origin[], float dir[], struct sphere *sphere, float *entry, float *exit)
{
	float origin_to_center[N];
	float origin_to_t_length;
	float t_radius2;
	float surface_to_t_length;

	sub(origin_to_center, sphere->center, origin);
	origin_to_t_length = dot(origin_to_center, dir);
	
	if (origin_to_t_length < 0)
		 return 0;

	/* Pythagoras */
	t_radius2 = length_squared(origin_to_center) - origin_to_t_length * origin_to_t_length;
	if (t_radius2 > sphere->radius2)
		 return 0;

	surface_to_t_length = sqrt(sphere->radius2 - t_radius2);

	if (entry)
		*entry = origin_to_t_length - surface_to_t_length;
	if (exit)
		*exit = origin_to_t_length + surface_to_t_length;
	
	return 1;
}

/*
int
find_first_intersection(struct sphere *spheres, unsgined int nspheres, float origin[], float dir[], struct sphere *intersected_sphere, float *intersected_entry) {	
	float entry;
	float exit;
	int ret = 0;

	*intersected_entry = INFINITY;
	for (i = 0; i < nspheres; ++i) {
		if (intersect(spheres[i], origin, dir, &entry, &exit)) {
			if (entry < 0)
				entry = exit;
	
			if (entry < *intersected_entry) {
				// We have found a "new first intersection" 
				ret = 1;
				*intersected_entry = entry;
				intersected_sphere = &spheres[i];
			}
		}
	}
	return ret;
}*/

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

void
generate_scene(struct sphere *spheres, int nspheres, unsigned int width,
unsigned int height, float width_inverse, float height_inverse, float
aspect_ratio, float angle)
{
	unsigned int i;
	float radius;
	
	/* Generate world sphere */
	spheres[0].center[0] = 0.0;
	spheres[0].center[0] = -10000.0;
	spheres[0].center[0] = -20.0;
	spheres[0].radius2 = 10000.0*10000.0;

	spheres[0].surface_colour[0] = 0.8;
	spheres[0].surface_colour[1] = 0.2;
	spheres[0].surface_colour[2] = 0.2;
	spheres[0].emission_colour[0] = 0.0;
	spheres[0].emission_colour[1] = 0.0;
	spheres[0].emission_colour[2] = 0.0;

	spheres[0].transparency = 0.0;
	spheres[0].reflection = 0.0;

	for (i = 0; i < nspheres - 1; ++i) {
		radius = randomf_in_range(0.1, 1.0);
		spheres[i].center[0] = randomf_in_range(-10.0, 10.0);
		spheres[i].center[1] = randomf_in_range(-10.0, 10.0);
		spheres[i].center[2] = randomf_in_range(-10.0, 10.0);
		spheres[i].radius2 = radius*radius;

		spheres[i].surface_colour[0] = randomf();
		spheres[i].surface_colour[1] = randomf();
		spheres[i].surface_colour[2] = randomf();

		spheres[i].emission_colour[0] = 0.0;
		spheres[i].emission_colour[1] = 0.0;
		spheres[i].emission_colour[2] = 0.0;
		/* SET_COLOUR(spheres[i].surface_colour, r, g, b);
		SET_COLOUR(spheres[i].emission_colour, 0.0, 0.0, 0.0); */

		spheres[i].transparency = randomf_in_range(0.1, 1.0);
		spheres[i].reflection = 1.0;
	}
	/* Add light source */
	spheres[nspheres - 1].center[0] = 0.0;
	spheres[nspheres - 1].center[1] = 20.0;
	spheres[nspheres - 1].center[2] = -30.0;
	spheres[nspheres - 1].radius2 = 3.0*3.0;

	spheres[nspheres - 1].surface_colour[0] = 0.0;
	spheres[nspheres - 1].surface_colour[1] = 0.0;
	spheres[nspheres - 1].surface_colour[2] = 0.0;
	spheres[nspheres - 1].emission_colour[0] = 3.0;
	spheres[nspheres - 1].emission_colour[1] = 3.0;
	spheres[nspheres - 1].emission_colour[2] = 3.0;

	spheres[nspheres - 1].transparency = 0.0;
	spheres[nspheres - 1].reflection = 0.0;
}

float c2cworld(unsigned int c, float measure_inverse)
{
	return 2.0 * (c + 0.5) * measure_inverse - 1.0;
}

float
y2yworld(unsigned int y, float height_inverse, float angle)
{
	return -1.0 * angle * c2cworld(y, height_inverse);
}

float
x2xworld(unsigned int x, float width_inverse, float angle, float aspect_ratio)
{
	return angle * aspect_ratio * c2cworld(x, width_inverse);
}

void
calculate_line(float *row, struct sphere *spheres, unsigned int nspheres, unsigned int y, unsigned int width, unsigned int height, float width_inverse, float height_inverse, float angle, float aspect_ratio)
{
	unsigned int x;
	float xworld;
	float yworld = y2yworld(y, height_inverse, angle);
	float origin[N] = {0.0, 0.0, 0.0};
	float dir[N];

	for (x = 0; x < width; ++x, row += 3) {
		/* We now have a set of coordinates (x, y),  we  project onto
		 * [0, angle*aspect_ratio] x [0, angle] */
		xworld = x2xworld(x, width_inverse, angle, aspect_ratio);
			
		/* Get direction of ray from camera */
		dir[0] = xworld;
		dir[1] = yworld;
		dir[2] = -1.0;
		normalize(dir);

		/* Trace ray */	
		trace(row, origin, dir, spheres, nspheres, 0);
	}
}

int
main(int argc, char **argv)
{
	/* Get provided support for threads, world size and world rank */
	int provided;
	int size;
	int rank;
	MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &provided);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	/* Seed random generator with lucky number */
	srand(13);

	unsigned int width = 1280;
	unsigned int height = 1024;
	float width_inverse = 1 / (float) width;
	float height_inverse = 1 / (float) height;
	float fov = 30;
	float aspect_ratio = width / (float) height;
	float angle = tan(M_PI * 0.5 * fov / 180.0);

	unsigned int nspheres = 100;
	struct sphere spheres[nspheres];
	MPI_Status status;

	/*
	std::vector<Sphere> spheres;
	if (rank == 0) {
			generate_scene(spheres, nbSpheres, width, height, invWidth, invHeight,
										 aspect_ratio, angle);
	}
	
	MPI_Bcast(spheres.data(), 3*nbSpheres, MPI_FLOAT, 0, MPI_COMM_WORLD);
	*/	
	generate_scene(spheres, nspheres, width, height, width_inverse, height_inverse, aspect_ratio, angle);
	int line;
	float row[3 * width];
	if (rank == 0) {
		/* Thy bidding, master? */
		unsigned int slave;
		unsigned int x;
		float image[3 * width * height];

		/* Send initial tasks */
		for (line = 0; line < size - 1 && line < height; ++line) {
			MPI_Send(&line, 1, MPI_INT, line + 1, 0, MPI_COMM_WORLD); 
		}
		for (; line < height; ++line) {
			MPI_Recv(row, 3*width, MPI_FLOAT, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);

			/* Memcpy the received stuff */
			for (x = 0; x < width; ++x)
				image[status.MPI_TAG * width + x] = row[x];

			/* Send more work */
			MPI_Send(&line, 1, MPI_INT, status.MPI_SOURCE, 0, MPI_COMM_WORLD);
		}

		/* Signal work done */

		for (slave = 1; slave < size; ++slave) {
			MPI_Recv(row, 3*width, MPI_FLOAT, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);

			/* Memcpy the received stuff */
			for (x = 0; x < width; ++x)
				image[status.MPI_TAG * width + x] = row[x];

			/* Send height to put slave to rest */
			MPI_Send(&height, 1, MPI_INT, slave, status.MPI_SOURCE, MPI_COMM_WORLD);
		}
		/* Save picture */
		/*std::ofstream ofs("./untitled.ppm", std::ios::out | std::ios::binary);
		ofs << "P6\n" << width << " " << height << "\n255\n";
		for (unsigned i = 0; i < width * height; ++i) {
				ofs << (unsigned char)(std::min(float(1), image[i].x) * 255) <<
							 (unsigned char)(std::min(float(1), image[i].y) * 255) <<
							 (unsigned char)(std::min(float(1), image[i].z) * 255);
		}
		ofs.close();*/

		/* Deallocate */

	} else {
		/* Work, work */
		do {
			MPI_Recv(&line, 1, MPI_INT, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
			//std::cout << "Got line " << line << std::endl;
			if (line < height) {
				calculate_line(row, spheres, nspheres, line, width, height, width_inverse, height_inverse, angle, aspect_ratio);
				MPI_Send(row, 3*width, MPI_FLOAT, 0, line, MPI_COMM_WORLD); 
			}
		} while(line < height);
	}
	MPI_Finalize();
	return 0;
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
max(float lhs, float rhs)
{
	if (lhs > rhs)
		return lhs;
	else
		return  rhs;
}

int
trace(float colour[], float ray_origin[], float ray_dir[], struct sphere *spheres, unsigned int nspheres, int depth)
{
	int i;
	struct sphere *sphere = NULL;
	float distance = INFINITY;
	float hit_point[N];
	float hit_normal[N];
	float bias = 1e-4;

	colour[0] = 0.0;
	colour[1] = 0.0;
	colour[2] = 0.0;

	for (i = 0; i < nspheres; ++i) {
		float entry;
		float exit;
		if (intersect(ray_origin, ray_dir, &spheres[i], &entry, &exit)) {
			if (entry < 0)
				entry = exit;

			if (entry < distance) {
				sphere = &sphere[i];
				distance = entry;
			}
		}
	}
	if (!sphere)
		return 0;
	/* Construct hit point and hit normal */
	construct(hit_point, ray_origin, ray_dir, distance);
	sub(hit_normal, hit_point, sphere->center);

	// if the object material is glass, split the ray into a reflection
	// and a refraction ray.
	if ((sphere->transparency > 0 || sphere->reflection > 0) && depth < MAX_RAY_DEPTH) {
		/* Compute reflection */
		float refraction_colour[3];
		float reflection_colour[3];

		float reflection_origin[N];
		float reflection_dir[N];

		float fresnel_effect;

		construct(reflection_origin, hit_point, hit_normal, bias);
		construct(reflection_dir, ray_dir, hit_normal, -2.0 * dot(ray_dir, hit_normal));

		trace(reflection_colour, reflection_origin, reflection_dir, spheres, nspheres, depth + 1);

		/* Compute refraction */
		if (sphere->transparency) {
			/* Refraction origin is the hit point of ray */
			float refraction_origin[N];
			float refraction_dir[N];
			int inside = 0;

			construct(refraction_origin, hit_point, hit_normal, -1.0 * bias);
			construct_refraction_dir(refraction_dir, ray_dir, hit_normal, inside);
			trace(refraction_colour, refraction_origin, refraction_dir, spheres, nspheres, depth + 1);
			scale(refraction_colour, sphere->transparency);
		}

		/* Calculate fresnel effect */
		fresnel_effect = fresnel(hit_normal, ray_dir, 0.1);
		mix(colour, reflection_colour, refraction_colour, fresnel_effect);
	} else {
		// object is a diffuse opaque object
		// compute illumination

		/* Find lights */
		for (i = 0; i < nspheres; ++i) {
			if (is_light(spheres[i].emission_colour)) {
				int j;
				float transmission_factor;
				float light_origin[N];
				float light_dir[N];
				
				construct(light_origin, hit_point, hit_normal, bias);
				sub(light_dir, spheres[i].center, hit_point);
				normalize(light_dir);
				transmission_factor = dot(hit_normal, light_dir);
				if (transmission_factor < 0)
					continue;

				for (j = 0; j < nspheres; ++j) {
					if (i != j && intersect(light_origin, light_dir, &sphere[j], NULL, NULL)) {
						transmission_factor = 0;
						break;
					}
				}
				construct(colour, colour, spheres[i].emission_colour, transmission_factor);
			}
		} 
	}
} 

