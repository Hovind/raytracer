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

void
float2bytes(unsigned char bytes[], float floats[], size_t n)
{
	size_t i;
	for (i = 0; i < n; ++i)
		bytes[i] = min(1.0, floats[i]) * 255.0;
}

void
save_ppm(char file_name[], float image[], unsigned int width, unsigned int height)
{
	size_t i;
	FILE *fp = fopen(file_name, "wb");
	fprintf(fp, "P6\n%u %u\n255\n", width, height);
	for (i = 0; i < height * width; ++i) {
		unsigned char bytes[3];
		float2bytes(bytes, image + 3 * i, 3);
		fwrite(bytes, sizeof(bytes[0]), 3, fp);
	}
	fclose(fp);
}
