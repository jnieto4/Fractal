#include <cstdio>
#include <algorithm>
#include <sys/time.h>
#include <mpi.h>
#include "BMP24.h"


static void fractal(const int width, unsigned char* const pic, int row_start, int row_end)
{
  	const double delta = 0.004737716;
  	const double xMid = -0.663889302;
  	const double yMid =  0.353461972;

  	// compute pixels of image
  	const double xMin = xMid - delta;
  	const double yMin = yMid - delta;
  	const double dw = 2.0 * delta / width;
  	double cy = yMin + dw * row_start;

	for (int row = row_start; row < row_end; row++) {  // rows
    		double cx = xMin;
    		for (int col = 0; col < width; col++) {  // columns
      			double x = cx;
      			double y = cy;
      			double x2, y2;
      			int count = 256;
      			do {
        			x2 = x * x;
        			y2 = y * y;
        			y = 2.0 * x * y + cy;
        			x = x2 - y2 + cx;
        			count--;
      			} while ((count > 0) && ((x2 + y2) < 5.0));
      			pic[row * width + col] = (unsigned char)count;
      			cx += dw;
    		}
    	cy = yMin + dw + row * dw;
  	}
}


int main(int argc, char* argv [])
{	
	// MIP Initaliation
	MPI_Init(NULL,NULL);
	int rank, size;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	
	if(rank == 0) {
  	printf("Fractal v2.4\n");

  	// check command line
  	if (argc != 2) {
		fprintf(stderr, "USAGE: %s image_width\n", argv[0]);
	       	exit(-1);
	}
	}
  	const int width = atoi(argv[1]);
	if(rank == 0) {
  	if (width < 12) {
		fprintf(stderr, "ERROR: image_width must be at least 12 pixels\n");
	       	exit(-1);
	}
  	printf("image width: %d\n", width);
}
  	// allocate image memory
  	unsigned char* pic = new unsigned char [width * width];
	unsigned char* final_pic = new unsigned char [width * width];

  	// start time
  	timeval beg, end;
	MPI_Barrier(MPI_COMM_WORLD);
  	gettimeofday(&beg, NULL);

	int row_start = rank * width / size;
	int row_end = row_start + (width / size);
	int chunk_size = (row_end - row_start) * width;
  	// execute timed code
  	fractal(width, pic, row_start, row_end);
	
	MPI_Gather(     &pic[(rank*chunk_size)], 
			chunk_size, 
			MPI_UNSIGNED_CHAR, 
			final_pic, 
			chunk_size, 
			MPI_UNSIGNED_CHAR,
		       	0, 
			MPI_COMM_WORLD);

	// end time
	gettimeofday(&end, NULL);

	if(rank == 0){
  		const double runtime = end.tv_sec - beg.tv_sec + (end.tv_usec - beg.tv_usec) / 1000000.0;
  		printf("compute time: %.6f s\n", runtime);
	
  		// write image to BMP file
  		if (width <= 1024) {
    			BMP24 bmp(0, 0, width, width);
    			for (int y = 0; y < width; y++) {
      				for (int x = 0; x < width; x++) {
        				bmp.dot(x, y, 0x0000ff - final_pic[y * width + x] * 0x000001 + final_pic[y * width + x] * 0x010100);
      				}
    			}
    			bmp.save("fractal.bmp");
  		}
	}

  	// clean up
  	delete [] pic;
	delete [] final_pic;
	MPI_Finalize();
  	return 0;
}
