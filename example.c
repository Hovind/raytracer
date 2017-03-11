__global__ void
kernel(unsigned char *ptr)
{
	// map from threadIdx/BlockIdx to pixel position
	int x = blockIdx.x;
	int y = blockIdx.y;
	int offset = x + y * gridDim.x;
	// now calculate the value at that position
	int juliaValue = julia( x, y );
	ptr[offset*4 + 0] = 255 * juliaValue;
	ptr[offset*4 + 1] = 0;
	ptr[offset*4 + 2] = 0;
	ptr[offset*4 + 3] = 255;
}

int
main(void)
{
	CPUBitmap bitmap(DIM, DIM);
	unsigned char *dev_bitmap;
	HANDLE_ERROR(cudaMalloc((void **) &dev_bitmap, bitmap.image_size()));
	dim3
	grid(DIM,DIM);
	kernel<<<grid,1>>>(dev_bitmap);
	HANDLE_ERROR(cudaMemcpy(bitmap.get_ptr(), dev_bitmap, bitmap.image_size(), cudaMemcpyDeviceToHost));
	bitmap.display_and_exit();
	HANDLE_ERROR(cudaFree(dev_bitmap));
}
