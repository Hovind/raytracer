#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <mpi.h>

#define N 3

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
	scale(1 / length(vec), vec);
}


struct sphere {
	float center[N];
	float radius;
	float surface_colour[3];
	float emission_colour[3];
	float transparency;
	float reflection;
};
	
int
intersect(struct sphere *sphere, float origin[], float dir[], float *entry, float *exit)
{
	float origin_to_center[N];
	float origin_to_t_length;
	float t_radius2;
	float surface_to_t_length;

	sub(origin_to_center, sphere->center, origin);
	origin_to_t_length = dot(origin_to_center, dir);
	
	if (origin_to_t_length < 0)
		 return false;

	/* Pythagoras */
	t_radius2 = length_square(origin_to_center) - origin_to_t_length * origin_to_t_length;
	if (t_radius2 > sphere->radius2)
		 return 0;

	surface_to_t_length = sqrt(sphere->radius2 - t_radius2);
	*entry = origin_to_t_length - surface_to_t_length;
	*exit = origin_to_t_length + surface_to_t_length;
	
	return 1;
}

#define MAX_RAY_DEPTH 5

float
mix(float a, float b, float mix)
{
    return b * mix + a * (1.0 - mix);
}

//[comment]
// This is the main trace function. It takes a ray as argument (defined by its origin
// and direction). We test if this ray intersects any of the geometry in the scene.
// If the ray intersects an object, we compute the intersection point, the normal
// at the intersection point, and shade this point using this information.
// Shading depends on the surface property (is it transparent, reflective, diffuse).
// The function returns a color for the ray. If the ray intersects an object that
// is the color of the object at the intersection point, otherwise it returns
// the background color.
//[/comment]

int
trace(float origin[], float dir[], struct sphere *spheres, unsigned int nspheres, int depth, float colour[])
{
	//if (raydir.length() != 1) std::cerr << "Error " << raydir << std::endl;
	float entry_min = INFINITY;
	struct sphere *sphere = NULL;
	float entry = INFINITY;
	float exit = INFINITY;
	unsigned int i;
	
	// find intersection of this ray with the sphere in the scene
	for (i = 0; i < nspheres; ++i) {
		if (intersect(spheres[i], origin, dir, &entry, &exit)) {
			if (entry < 0)
				t0 = t1;
	
			if (entry < entry_min) {
			        entry_min = entry;
			        sphere = &spheres[i];
			}
		}
	}
	// if there's no intersection return black or background color
	if (!sphere)
		return 0;
	surface_colour[3] = {0.0, 0.0, 0.0}; // color of the ray/surfaceof the object intersected by the ray
	
	copy(origin_to_hit, dir);
	scale(origin_to_hit, entry);
	
	add(hit_point, origin, origin_to_hit);
	sub(hit_normal, hit_point, sphere->center);
	normalize(hit_normal);

	/* Add some bias to the point from which we will be tracing */
	bias = 1e-4;
	inside = 0;

	if (dot(dir, nhit) > 0) {
		nhit = -nhit;
		inside = 1;
	}

	if ((sphere->transparency > 0 || sphere->reflection > 0) && depth < MAX_RAY_DEPTH) {
		facing_ratio = dot(dir, hit_normal);

		/* Change the mix value to tweak the effect */
		fresnel_effect = mix(pow(1 - facingratio, 3), 1.0, 0.1);

		/* Compute reflection direction, all vectors are already normalized */
		copy(reflection_dir, hit_normal);
		scale(reflection_dir, 2 * dot(dir, reflection_dir));
		sub(reflection_dir, dir, reflection_dir);
		normalize(reflection_dir);

		copy(tmp, hit_normal);
		scale(tmp, bias);
		sum(tmp, hit_point, tmp);
		trace(reflection, tmp, reflection_dir, spheres, depth + 1);
		refraction = 0;

		Vec3f refraction = 0;
		// if the sphere is also transparent compute refraction ray (transmission)
		if (sphere->transparency) {
			float ior = 1.1, eta = (inside) ? ior : 1 / ior; // are we inside or outside the surface?
			float cosi = -dot(hit_normal, dir);
			float k = 1 - eta * eta * (1 - cosi * cosi);
			Vec3f refrdir = raydir * eta + nhit * (eta *  cosi - sqrt(k));
			normalize(refraction_dir);
			copy(tmp, hit_normal);
			scale(tmp, bias);
			sub(tmp, hit_point, tmp);
			trace(refraction, tmp, refraction_dir, spheres, depth + 1);
		}
		// the result is a mix of reflection and refraction (if the sphere is transparent)
		surface_colour = sphere->surface_colour * mix(reflection, refraction * sphere->transparency, fresnel_effect);
		//surfaceColor = (reflection * fresneleffect + refraction * (1 - fresneleffect) * sphere->transparency) * sphere->surfaceColor;
	} else {
	        // it's a diffuse object, no need to raytrace any further
		for (i = 0; i < nspheres; ++i) {
			/* Check if light */
			if (spheres[i].emission_colour[0] > 0) {
				transmission = {1.0, 1.0, 1.0};
				sub(light_dir, spheres[i].center, hit_point);
				normalize(light_dir);
				//Vec3f lightDirection = spheres[i].center - phit;
				//lightDirection.normalize();
				for (unsigned j = 0; j < spheres.size(); ++j) {
					if (i != j) {
						if (intersect(spheres[j], phit + nhit * bias, lightDirection, NULL, NULL)) {
							transmission = 0;
							break;
						}
					}
				}
				surface_colour += sphere->surface_colour * transmission * spheres[i].emission_colour * MAX(0.0, dot(nhit, light_dir));
			}
		}
	}
    add(surface_colour, sphere->emission_colour, colour);
    return 1;
}

float
randomf(void) {
	return rand() / float(RAND_MAX);
}

float
randomf_in_range(float min, float max)
{
	return randomf() * (max - min) + min;
}

void
generate_scene(struct sphere *spheres, int nspheres)
{
	unsigned int i;
	float x;
	float y;
	float z;
	float radius;
	float r;
	float g;
	float b;
	float t;

	/* Seed random generator with lucky number */
	srand48(13);
	
	/* Generate world sphere */
	spheres[0] = {
		.center = {0, -10000.0, -20.0},
		.radius = 10000,
		.surface_colour = {0.80, 0.20, 0.20},
		.emission_colour = {0.0, 0.0, 0.0},
		.transparency = 0.0,
		.reflection = 0.0,
	};

	for (i = 0; i < nspheres - 1; ++i) {
		x = randomf_in_range(-10.0, 10.0);
		y = randomf_in_range(-10.0, 10.0);
		z = randomf_in_range(-10.0, 10.0);

		radius = randomf_in_range(0.1, 1.0);
		r = randomf();
		g = randomf();
		b = randomf();

		t = random_in_range(0.0, 0.5);

		spheres[i] = {
			.center = {x, y, z},
			.radius = radius,
			.surface_colour = {r, g, b},
			.emission_colour = {0.0, 0.0, 0.0},
			.transparency = t,
			.reflection = 1.0,
		};
	}
	/* Add light source */
	spheres[nspheres - 1] = {
		.center = {0.0, 20.0, -30),
		.radius = 3,
		.surface_colour = {0.0, 0.0, 0.0},
		.emisson_colour = {3.0, 3.0, 3.0},
		.transparency = 0.0,
		.reflection = 0.0,
	};
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
	return angle * aspect_ratio * c2cworld(y, height_inverse);
}

{
	

void calculate_line(float *row, struct sphere *spheres, unsigned int nspheres, unsigned int y, unsigned int width, unsigned int height, float width_inverse, float height_inverse, float angle, float aspect_ratio)
{
	unsigned int x;
	float xworld;
	float yworld = y2yworld(y, height_inverse, angle);
	float zeros[N] = {0.0, 0.0, 0.0};

	for (x = 0; x < width; ++x, row += 3) {
		/* We now have a set of coordinates (x, y),  we  project onto
		 * [0, angle*aspectratio] x [0, angle] */
		xworld = x2xworld(x, width_inverse, angle, aspect_ratio);
			
		/* Get direction of ray from camera */
		dir = {xworld, yworld, -1.0};
		normalize(dir);

		/* Trace ray */	
		trace(row, origin, dir, spheres, 0);
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

	unsigned int width = 1280;
	unsigned int height = 1024;
	float width_inverse = 1 / float(width);
	float height_inverse = 1 / float(height);
	float fov = 30;
	float aspect_ratio = width / float(height);
	float angle = tan(M_PI * 0.5 * fov / 180.0);

	unsigned int nbSpheres = 100;
	struct sphere spheres[nspheres];
	MPI_Status status;

	/*
	std::vector<Sphere> spheres;
	if (rank == 0) {
			generate_scene(spheres, nbSpheres, width, height, invWidth, invHeight,
										 aspectratio, angle);
	}
	
	MPI_Bcast(spheres.data(), 3*nbSpheres, MPI_FLOAT, 0, MPI_COMM_WORLD);
	*/	
	generate_scene(spheres, nspheres, width, height, width_inverse, height_inverse, aspect_ratio, angle);
	int line;
	float row[3 * width];
	if (rank == 0) {
		/* Thy bidding, master? */
		float image[3 * width * height];

		/* Send initial tasks */
		for (line = 0; line < size - 1 && line < height; ++line) {
			MPI_Send(&line, 1, MPI_INT, line + 1, 0, MPI_COMM_WORLD); 
		}
		for (; line < height; ++line) {
			MPI_Recv(row, 3*width, MPI_FLOAT, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);

			/* Memcpy the received stuff */
			for (int x = 0; x < width; ++x)
				image[status.MPI_TAG * width + x] = row[x];

			/* Send more work */
			MPI_Send(&line, 1, MPI_INT, status.MPI_SOURCE, 0, MPI_COMM_WORLD);
		}

		/* Signal work done */

		for (int slave = 1; slave < size; ++slave) {
			MPI_Recv(row, 3*width, MPI_FLOAT, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);

			/* Memcpy the received stuff */
			for (int x = 0; x < width; ++x)
				image[status.MPI_TAG * width + x] = row[x];

			/* Send height to put slave to rest */
			MPI_Send(&height, 1, MPI_INT, slave, status.MPI_SOURCE, MPI_COMM_WORLD);
		}
		/* Save picture */
		std::ofstream ofs("./untitled.ppm", std::ios::out | std::ios::binary);
		ofs << "P6\n" << width << " " << height << "\n255\n";
		for (unsigned i = 0; i < width * height; ++i) {
				ofs << (unsigned char)(std::min(float(1), image[i].x) * 255) <<
							 (unsigned char)(std::min(float(1), image[i].y) * 255) <<
							 (unsigned char)(std::min(float(1), image[i].z) * 255);
		}
		ofs.close();

		/* Deallocate */
		delete [] image;

	} else {
		/* Work, work */
		do {
			MPI_Recv(&line, 1, MPI_INT, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
			//std::cout << "Got line " << line << std::endl;
			if (line < height) {
				calculate_line(row, spheres, line, width, height, width_inverse, height_inverse, angle, aspect_ratio);
				MPI_Send(row, 3*width, MPI_FLOAT, 0, line, MPI_COMM_WORLD); 
			}
		} while(line < height);
	}
	delete [] row;
	MPI_Finalize();
	return 0;
}


