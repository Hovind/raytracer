#include <time.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#define N 5

void
copy(float src[], float to[])
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


float
randomf(void)
{
        return rand() / float(RAND_MAX);
}


int
main()
{
	srand(time(NULL));
	float src[N] = {1.0, 2.0, 3.0, 4.0, 5.0};
	float to[N];
	print(to);
	copy(src, to);
	print(to);
	printf("Float: %f\n", randomf());
	return 0;
}
