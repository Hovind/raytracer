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

int
main()
{
	float src[N] = {1.0, 2.0, 3.0, 4.0, 5.0};
	float to[N];
	print(to);
	copy(src, to);
	print(to);
	return 0;
}
