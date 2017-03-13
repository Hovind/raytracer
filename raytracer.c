#include <pthread.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <mpi.h>

#define N 3
#define MAX_RAY_DEPTH 5
#define M_PI 3.14159265358979323846
#define BIAS 1e-4

struct screen_config {
	float width_inverse;
	float height_inverse;
	float fov;
	float aspect_ratio;
	float angle;
};
	
struct sphere {
	float center[N];
	float radius2;
	float surface_colour[3];
	float emission_colour[3];
	float transparency;
	float reflection;
};

typedef struct {
	int provided;
	int size;
	int rank;
} MPI_Context;

void
MPI_Init_context(int *argc, char ***argv, MPI_Context *context)
{	
	MPI_Init_thread(argc, argv, MPI_THREAD_FUNNELED, &context->provided);
	MPI_Comm_size(MPI_COMM_WORLD, &context->size);
	MPI_Comm_rank(MPI_COMM_WORLD, &context->rank);
}


void
construct_refraction_unit_ray(struct ray refraction_ray, r_dir);
	k = 1 - eta * eta * (1 - cosi * cosi);

	add_weighted(refraction_dir, ray_dir, eta, hit_normal, eta * cosi - sqrt(k));
	normalize(refraction_dir);
}

int
intersect(struct ray *r, struct sphere *sphere, float *entry, float *exit)
{
	float origin_to_center[N];
	float origin_to_t_length;
	float t_radius2;
	float surface_to_t_length;

	sub(origin_to_center, sphere->center, r->origin);
	origin_to_t_length = dot(origin_to_center, r->dir);
	
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

struct sphere *
generate_scene(unsigned int nspheres, unsigned int width, unsigned int height, struct screen_config *screen)
{
	unsigned int i;
	float radius;
	struct sphere *spheres = malloc(nspheres * sizeof(*spheres));

	screen->width_inverse = 1 / (float) width;
	screen->height_inverse = 1 / (float) height;
	screen->fov = 30.0;
	screen->aspect_ratio = width / (float) height;
	screen->angle = tanf(M_PI * 0.5 * screen->fov / 180.0);
	
	/* Generate world sphere */
	spheres[0].center[0] = 0.0;
	spheres[0].center[1] = -10004.0;
	spheres[0].center[2] = -20.0;
	spheres[0].radius2 = 10000.0*10000.0;

	set_colour(spheres[0].surface_colour, 0.8, 0.2, 0.2);
	set_colour(spheres[0].emission_colour, 0.0, 0.0, 0.0);

	spheres[0].transparency = 0.0;
	spheres[0].reflection = 0.0;

	for (i = 1; i + 1 < nspheres; ++i) {
		spheres[i].center[0] = randomf_in_range(-10.0, 10.0);
		spheres[i].center[1] = randomf_in_range(-1.0, 1.0);
		spheres[i].center[2] = randomf_in_range(-25.0, -15.0);
		radius = randomf_in_range(0.9, 1.0);
		spheres[i].radius2 = radius*radius;

		set_colour(spheres[i].surface_colour, randomf(), randomf(), randomf());
		set_colour(spheres[i].emission_colour, 0.0, 0.0, 0.0);

		spheres[i].transparency = randomf_in_range(0.0, 0.5);
		spheres[i].reflection = 1.0;
	}
	/* Add light source */
	spheres[nspheres - 1].center[0] = 0.0;
	spheres[nspheres - 1].center[1] = 20.0;
	spheres[nspheres - 1].center[2] = -30.0;
	spheres[nspheres - 1].radius2 = 3.0*3.0;

	set_colour(spheres[nspheres - 1].surface_colour, 0.0, 0.0, 0.0);
	set_colour(spheres[nspheres - 1].emission_colour, 3.0, 3.0, 3.0);

	spheres[nspheres - 1].transparency = 0.0;
	spheres[nspheres - 1].reflection = 0.0;

	return spheres;
}

float c2cworld(unsigned int c, float measure_inverse)
{
	return 2.0 * (c + 0.5) * measure_inverse - 1.0;
}

float
y2yworld(unsigned int y, struct screen_config screen)
{
	return -1.0 * screen.angle * c2cworld(y, screen.height_inverse);
}

float
x2xworld(unsigned int x, struct screen_config screen)
{
	return screen.angle * screen.aspect_ratio * c2cworld(x, screen.width_inverse);
}

int
ray_reaches(struct ray *r, unsigned int i, struct sphere *spheres, unsigned int nspheres)
{
	unsigned int j;
	for (j = 0; j < nspheres; ++j)
		if (i != j && intersect(r, spheres + j, NULL, NULL))
			return 0;
	return 1;
}

int
find_first_intersection(float *distance, struct sphere *sphere, struct ray *camera_ray, struct sphere *spheres, unsigned int nspheres)
{ 
	for (i = 0, *distance = INFINITY; i < nspheres; ++i) {
		float entry;
		float exit;
		if (intersect(camera_ray, spheres + i, &entry, &exit)) {
			if (entry < 0)
				entry = exit;
	
			if (entry < *distance) {
				sphere = spheres + i;
				*distance = entry;
			}
		}
	}
}

int
trace(float colour[], struct ray *r, struct sphere *spheres, unsigned int nspheres, int depth)
{
	unsigned int i;
	int inside; 
	struct sphere *sphere;
	struct ray hit_normal_ray;
	float distance;

	if (find_first_intersection(&distance, sphere, camera_ray, spheres, nspheres)) {
		set_colour(colour, 2.0, 2.0, 2.0);
		return 0;
	}

	set_colour(colour, 2.0, 2.0, 2.0);

	/* Construct hit point and hit normal */
	construct_hit_normal_unit_ray(hit_normal_ray, distance);
	inside = flip_ray_if_inside(&hit_normal_ray, r);

	if ((sphere->transparency > 0 || sphere->reflection > 0) && depth < MAX_RAY_DEPTH) {
		/* Compute reflection */
		unsigned char refraction_colour[3];
		unsigned char reflection_colour[3];
		struct ray reflection_ray;
		struct ray reflection_ray;
		float facing_ratio;
		float fresnel_effect;

		set_colour(refraction_colour, 0.0, 0.0, 0.0);
		set_colour(reflection_colour, 0.0, 0.0, 0.0);

		facing_ratio = - 1.0 * dot(ra, hit_normal);
		construct_reflection_ray(reflection_ray, hit_point, hit_normal, bias); ...
		trace(reflection_ray);
		scale_colour(relfection_colour, fresnel_effect);

		/* Compute refraction */
		if (sphere->transparency) {
			/* Refraction origin is the hit point of ray */

			struct ray reflection_ray;
			construct_reflection_ray(reflection_ray);
			trace();
			scale_colour(refraction_colour, (1.0 - fresnel_effect) * sphere->transparency);
		}

		/* Calculate fresnel effect */
		add_colours(colour, reflection_colour, refraction_colour);
		mul_colours(colour, colour, sphere->surface_colour);
	} else {
		/* Find lights */
		for (i = 0; i < nspheres; ++i) {
			if (is_light(spheres[i].emission_colour)) {
				struct ray light_ray;
				float transmission_factor;
				
				construct_light_unit_ray();
				transmission_factor = dot(hit_normal_ray.dir, light_ray.dir);
				if (transmission_factor > 0 && ray_reaches(&light_ray, i, spheres, nspheres)) {
					float product_colour[N];

					mul_colours(product_colour, sphere->surface_colour, spheres[i].emission_colour);
					construct(colour, colour, product_colour, transmission_factor);
				}
			}
		} 
	}
	add(colour, colour, sphere->emission_colour);
	return 1;
}

struct segment_args {
	float *row;
	size_t i;
	pthread_mutex_t m;
	size_t segment_length;
	float yworld;
	struct sphere *spheres;
	unsigned int nspheres;
	struct screen_config screen;
};

void *
calculate_segment(void *vargs)
{
	struct segment_args *args;
	struct ray camera_ray;
	size_t j;
	size_t x;
	float xworld;

	args = vargs;

	pthread_mutex_lock(&args->m);
	j = args->i++;
	pthread_mutex_unlock(&args->m);

	
	for (x = j * args->segment_length; x < (j + 1) * args->segment_length; ++x) {
		xworld = x2xworld(x, args->screen);
		construct_camera_ray(&camera_ray, xworld, yworld);
		trace(args->row + 3 * x, &camera_ray, args->spheres, args->nspheres, 0);
	}
	return NULL;
}
		
void
calculate_line(float *row, size_t y, size_t width, size_t nsegments, struct sphere *spheres, unsigned int nspheres, struct screen_config screen)
{
	size_t i;

	pthread_t *threads = malloc(nsegments * sizeof(*threads));

	float yworld = y2yworld(y, screen);
	struct segment_args args = {
		.row = row,
		.i = 0,
		.segment_length = width / nsegments,
		.yworld = yworld,
		.spheres = spheres,
		.nspheres = nspheres,
		.screen = screen,
	};

	for (i = 1; i < nsegments; ++i)
		pthread_create(threads + i, NULL, calculate_segment, &args);

	/* Do some work on master thread as well, lead by example */
	calculate_segment(&args);

	for (i = 1; i < nsegments; ++i)
		pthread_join(threads[i], NULL);

	free(threads);
}

int
main(int argc, char **argv)
{
	/* Get provided support for threads, world size and world rank */
	MPI_Status status;
	MPI_Context context;

	size_t line;
	float *row;
	struct sphere *spheres;	
	struct screen_config screen;

	size_t width = 1280;
	size_t height = 1024;
	unsigned int nspheres = 100;
	size_t nsegments = 2;

	srand(13);
	MPI_Init_context(&argc, &argv, &context);

	spheres = generate_scene(nspheres, width, height, &screen);
	row = malloc(3 * width * sizeof(*row));

	if (context.rank == 0) {
		int slave;
		float *image = malloc(3 * width * height * sizeof(*image));

		for (line = 0; line < context.size - 1 && line < height; ++line)
			MPI_Send(&line, 1, MPI_INT, line + 1, 0, MPI_COMM_WORLD); 

		for (; line < height; ++line) {
			MPI_Recv(row, 3 * width, MPI_FLOAT, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
			memcpy(image + 3 * status.MPI_TAG * width, row, 3 * width * sizeof(*image));
			MPI_Send(&line, 1, MPI_INT, status.MPI_SOURCE, 0, MPI_COMM_WORLD);
		}

		for (slave = 1; slave < context.size; ++slave) {
			MPI_Recv(row, 3 * width, MPI_FLOAT, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
			memcpy(image + 3 * status.MPI_TAG * width, row, 3 * width * sizeof(*image));
			MPI_Send(&height, 1, MPI_INT, slave, status.MPI_SOURCE, MPI_COMM_WORLD);
		}
		save_ppm("untitled.ppm", image, width, height);
		free(image);

	} else {
		/* Work, work */
		while (MPI_Recv(&line, 1, MPI_INT, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &status), line < height) {
			calculate_line(row, line, width, nsegments, spheres, nspheres, screen);
			MPI_Send(row, 3*width, MPI_FLOAT, 0, line, MPI_COMM_WORLD); 
		}
	}
	/* Deallocate */
	free(row);
	free(spheres);

	/* Deinit MPI */
	MPI_Finalize();
	return 0;
}
