#include <stdbool.h> 

typedef struct {
   unsigned char red;
   unsigned char green;
   unsigned char blue;
} pixel;

typedef struct {
    int red;
    int green;
    int blue;
    // int num;
} pixel_sum;

/*
* General Optimization we do:
1) x++ -> ++x
*/

/* Compute min and max of two integers, respectively */
int min(int a, int b) { return (a < b ? a : b); }
int max(int a, int b) { return (a > b ? a : b); }

int calcIndex(int i, int j, int n) {
    return ((i)*(n)+(j));
}

/*
 * initialize_pixel_sum - Initializes all fields of sum to 0
 */
void initialize_pixel_sum(pixel_sum *sum) {
	//Optimize: not go to memory to read:
    sum->red = 0;
	sum->green = 0;
	sum->blue = 0;
    // sum->num = 0;
    return;
}

/*
 * assign_sum_to_pixel - Truncates pixel's new value to match the range [0,255]
 */
static void assign_sum_to_pixel(pixel *current_pixel, pixel_sum sum, int kernelScale) {

    // divide by kernel's weight
    sum.red /= kernelScale;
    sum.green /= kernelScale;
    sum.blue /= kernelScale;

    // truncate each pixel's color values to match the range [0,255]
    current_pixel->red = (unsigned char) (min(max(sum.red, 0), 255));
    current_pixel->green = (unsigned char) (min(max(sum.green, 0), 255));
    current_pixel->blue = (unsigned char) (min(max(sum.blue, 0), 255));
    return;
}

/*
* sum_pixels_by_weight - Sums pixel values, scaled by given weight
*/
static void sum_pixels_by_weight(pixel_sum *sum, pixel p, int weight) {
    sum->red += ((int) p.red) * weight;
    sum->green += ((int) p.green) * weight;
    sum->blue += ((int) p.blue) * weight;
    // sum->num++;
    return;
}

/*
 *  Applies kernel for pixel at (i,j)
 */
static pixel applyKernel(int dim, int i, int j, pixel *src, int kernelSize, int kernel[kernelSize][kernelSize], int kernelScale, bool filter) {

    int ii, jj;
    int currRow, currCol;
    pixel_sum sum;
    pixel current_pixel;
    int min_intensity = 766; // arbitrary value that is higher than maximum possible intensity, which is 255*3=765
    int max_intensity = -1; // arbitrary value that is lower than minimum possible intensity, which is 0
    int min_index, max_index;// more efficient
    pixel loop_pixel;

	//this pixel is in the range of change - func smooth takes care of that
	//so, we don't need min or max

    initialize_pixel_sum(&sum);

	//calculating lengths:
	int halfKernelSize = kernelSize >> 1; //=kernelSize/2
	int length1 = i + halfKernelSize;
	int length2 = j + halfKernelSize;

	//Optimize by calculate index efficiently in loop:
	int kIndex = 0;
	ii = i - halfKernelSize;
	int indexPixel = ii * dim;

	//best to check here and calculate to filter in his loop for minimize the memory we read.
	if(!filter){
		for(; ii <= length1; ++ii) {
			jj = j - halfKernelSize;
			for(; jj <= length2; ++jj) {
				// apply kernel on pixel at [ii,jj]
				sum_pixels_by_weight(&sum, src[indexPixel + jj], *((int *)kernel +kIndex));
				++kIndex;
			}
			indexPixel += dim;
		}

		// assign kernel's result to pixel at [i,j]
		assign_sum_to_pixel(&current_pixel, sum, kernelScale);
		return current_pixel;
	}

    for(; ii <= length1; ++ii) {
		jj = j - halfKernelSize;
        for(; jj <= length2; ++jj) {
			// check if smaller than min or higher than max and update
			int finalIndex = indexPixel + jj;
            loop_pixel = src[finalIndex];

			// apply kernel on pixel at [ii,jj]
			sum_pixels_by_weight(&sum, src[finalIndex], *((int *)kernel +kIndex));
			++kIndex;
			
			//calculating once this arg:
			int loop_pixel_intensity = ((int) loop_pixel.red) + ((int) loop_pixel.green) + ((int) loop_pixel.blue);

            if (loop_pixel_intensity <= min_intensity) {
                min_intensity = loop_pixel_intensity;

				//more efficient to remember the loc:
                min_index = finalIndex;
            }

            if (loop_pixel_intensity > max_intensity) {
                max_intensity = loop_pixel_intensity;

				//more efficient to remember the loc:
                max_index = finalIndex;
            }
		}
		indexPixel += dim;
	}

    // filter out min and max
    sum_pixels_by_weight(&sum, src[min_index], -1);
    sum_pixels_by_weight(&sum, src[max_index], -1);

    // assign kernel's result to pixel at [i,j]
    assign_sum_to_pixel(&current_pixel, sum, kernelScale);
    return current_pixel;
}

/*
* Apply the kernel over each pixel.
* Ignore pixels where the kernel exceeds bounds. These are pixels with row index smaller than kernelSize/2 and/or
* column index smaller than kernelSize/2
*/
void smooth(int dim, pixel *src, pixel *dst, int kernelSize, int kernel[kernelSize][kernelSize], int kernelScale, bool filter) {

    int i, j;
	//removing calculations from loop:
	int halfKernelSize = kernelSize << 1;
	int length = dim - halfKernelSize;
	//calculating index efficiently
	int dstIndex = halfKernelSize * dim;
    for (i = halfKernelSize ; i < length; ++i) {
        for (j =  halfKernelSize ; j < length; ++j) {
            dst[dstIndex + j] = applyKernel(dim, i, j, src, kernelSize, kernel, kernelScale, filter);
        }

		dstIndex += dim;
    }
}

void charsToPixels(Image *charsImg, pixel* pixels) {
//calculating index efficiently
    int row = 0;
	int index = 0;
	int indexMult3 = 0;
    for (row = 0 ; row < m ; ++row) {
		int col = 0;
        for (; col < n ; ++col) {

            pixels[index].red = image->data[indexMult3];
            pixels[index].green = image->data[indexMult3 + 1];
            pixels[index].blue = image->data[indexMult3 + 2];

			++index;
			indexMult3 += 3;
        }
    }
}

void pixelsToChars(pixel* pixels, Image *charsImg) {
	//calculating index efficiently
    int row = 0;
	int index = 0;
	int indexMult3 = 0;
    for (row = 0 ; row < m ; ++row) {
		int col = 0;
        for (; col < n ; ++col) {

            image->data[indexMult3] = pixels[index].red;
            image->data[indexMult3 + 1] = pixels[index].green;
            image->data[indexMult3 + 2] = pixels[index].blue;

			++index;
			indexMult3 += 3;
        }
    }
}

void copyPixels(pixel* src, pixel* dst) {
	//calculating index efficiently
    int row = 0;
	int index = 0;
    for (row = 0 ; row < m ; ++row) {
		int col = 0;
        for (; col < n ; ++col) {

            dst[index].red = src[index].red;
            dst[index].green = src[index].green;
            dst[index].blue = src[index].blue;

			++index;
        }
    }
}

void doConvolution(Image *image, int kernelSize, int kernel[kernelSize][kernelSize], int kernelScale, bool filter) {
	//caculating length once:
	size_t size = m*n*sizeof(pixel);
    pixel* pixelsImg = malloc(size);
    pixel* backupOrg = malloc(size);

    charsToPixels(image, pixelsImg);
    copyPixels(pixelsImg, backupOrg);

    smooth(m, backupOrg, pixelsImg, kernelSize, kernel, kernelScale, filter);

    pixelsToChars(pixelsImg, image);

    free(pixelsImg);
    free(backupOrg);
}

void myfunction(Image *image, char* srcImgpName, char* blurRsltImgName, char* sharpRsltImgName, char* filteredBlurRsltImgName, char* filteredSharpRsltImgName, char flag) {

    /*
    * [1, 1, 1]
    * [1, 1, 1]
    * [1, 1, 1]
    */
    int blurKernel[3][3] = {{1, 1, 1}, {1, 1, 1}, {1, 1, 1}};

    /*
    * [-1, -1, -1]
    * [-1, 9, -1]
    * [-1, -1, -1]
    */
    int sharpKernel[3][3] = {{-1,-1,-1},{-1,9,-1},{-1,-1,-1}};

    if (flag == '1') {  
        // blur image
        doConvolution(image, 3, blurKernel, 9, false);

        // write result image to file
        writeBMP(image, srcImgpName, blurRsltImgName);  

        // sharpen the resulting image
        doConvolution(image, 3, sharpKernel, 1, false);
        
        // write result image to file
        writeBMP(image, srcImgpName, sharpRsltImgName); 
    } else {
        // apply extermum filtered kernel to blur image
        doConvolution(image, 3, blurKernel, 7, true);

        // write result image to file
        writeBMP(image, srcImgpName, filteredBlurRsltImgName);

        // sharpen the resulting image
        doConvolution(image, 3, sharpKernel, 1, false);

        // write result image to file
        writeBMP(image, srcImgpName, filteredSharpRsltImgName); 
    }
}


