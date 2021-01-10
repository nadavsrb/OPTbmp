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
    if(kernelScale != 1) {  //often the scale is 1.
        // divide by kernel's weight
        sum.red /= kernelScale;
        sum.green /= kernelScale;
        sum.blue /= kernelScale;
    }

    // truncate each pixel's color values to match the range [0,255]
    // this func is called a lot we won't want she also calls funcs.
    if (sum.red < 0){sum.red = 0;} else if(sum.red > 255) {sum.red = 255;}
    current_pixel->red = (unsigned char) sum.red;
    if (sum.green < 0){sum.green = 0;} else if(sum.green > 255) {sum.green = 255;}
    current_pixel->green = (unsigned char) sum.green;
    if (sum.blue < 0){sum.blue = 0;} else if(sum.blue > 255) {sum.blue = 255;}
    current_pixel->blue = (unsigned char) sum.blue;
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
    pixel_sum sum;
    pixel current_pixel;
    int min_intensity = 766; // arbitrary value that is higher than maximum possible intensity, which is 255*3=765
    int max_intensity = -1; // arbitrary value that is lower than minimum possible intensity, which is 0
    pixel min_p, max_p;// more efficient
    pixel loop_pixel;

	// we check that this pixel is in the range of change (in the "red square as defined in ex5")
	int halfKernelSize = kernelSize >> 1; //=kernelSize/2
	int endRange = dim - halfKernelSize;
	if( i < halfKernelSize ||
		i >= endRange ||
		j < halfKernelSize ||
		j >= endRange){
			return *(src+ i*dim +j);
	}
	//so, we don't need min or max

    // func for this little thing isn't effcient:
    // initialize_pixel_sum(&sum);
    sum.red = 0;
	sum.green = 0;
	sum.blue = 0;

	//calculating lengths:
	int length1 = i + halfKernelSize;
	int length2 = j + halfKernelSize;

	//Optimize by calculate index efficiently in loop:
	int* kIndex = (int *) kernel;
	ii = i - halfKernelSize;
	int indexPixel = ii * dim;
    int startJJ = j - halfKernelSize;
	//best to check here and calculate to filter in his loop for minimize the memory we read.
	if(!filter){
		for(; ii <= length1; ++ii) {
			jj = startJJ;
            int finalIndex = indexPixel + jj;
			for(; jj <= length2; ++jj) {
				// apply kernel on pixel at [ii,jj]

				// Optimization don't use func in this loop instead (many things in stack):
				pixel p = src[finalIndex];
				int weight = *(kIndex);
				sum.red += ((int) p.red) * weight;
    			sum.green += ((int) p.green) * weight;
    			sum.blue += ((int) p.blue) * weight;

				++kIndex;
                ++finalIndex;
			}
			indexPixel += dim;
		}
	} else {
        for(; ii <= length1; ++ii) {
            jj = startJJ;
            int finalIndex = indexPixel + jj;
            for(; jj <= length2; ++jj) {
                // check if smaller than min or higher than max and update
                loop_pixel = src[finalIndex];
                int r = (int) loop_pixel.red;
                int g = (int) loop_pixel.green;
                int b = (int) loop_pixel.blue;
                // apply kernel on pixel at [ii,jj]
                // Optimization don't use func in this loop instead (many things in stack):
                int weight = *(kIndex);
                sum.red += r * weight;
                sum.green += g * weight;
                sum.blue += b * weight;
                
                ++kIndex;
                
                //calculating once this arg:
                int loop_pixel_intensity = r + b + g;

                if (loop_pixel_intensity <= min_intensity) {
                    min_intensity = loop_pixel_intensity;

                    //more efficient to remember the loc:
                    min_p = loop_pixel;
                }

                if (loop_pixel_intensity > max_intensity) {
                    max_intensity = loop_pixel_intensity;

                    //more efficient to remember the loc:
                    max_p = loop_pixel;
                }

                ++finalIndex;
            }
            indexPixel += dim;
        }

        // filter out min and max
        // func are not effecient:
        // sum_pixels_by_weight(&sum, min_p, -1);
        // sum_pixels_by_weight(&sum, max_p, -1);
        sum.red -= min_p.red + max_p.red;
        sum.green -= min_p.green + max_p.green;
        sum.blue -= min_p.blue + max_p.blue;
    }

    // this func is called a lot so we won't use as func - heavy on stack
	// assign_sum_to_pixel(&current_pixel, sum, kernelScale);
    if(kernelScale != 1) {  //often the scale is 1.
        // divide by kernel's weight
        sum.red /= kernelScale;
        sum.green /= kernelScale;
        sum.blue /= kernelScale;
    }

    // truncate each pixel's color values to match the range [0,255]
    // this func is called a lot we won't want she also calls funcs.
    if (sum.red < 0){sum.red = 0;} else if(sum.red > 255) {sum.red = 255;}
    current_pixel.red = (unsigned char) sum.red;
    if (sum.green < 0){sum.green = 0;} else if(sum.green > 255) {sum.green = 255;}
    current_pixel.green = (unsigned char) sum.green;
    if (sum.blue < 0){sum.blue = 0;} else if(sum.blue > 255) {sum.blue = 255;}
    current_pixel.blue = (unsigned char) sum.blue;

    return current_pixel;
}

/*
* Apply the kernel over each pixel.
* doesn't Ignore pixels where the kernel exceeds bounds-> applyKernel ignores.
*/
void smooth(int dim, pixel *src, pixel *dst, int kernelSize, int kernel[kernelSize][kernelSize], int kernelScale, bool filter) {

    int i, j;
	//calculating index efficiently
	int dstIndex = dim * dim - 1;
    int start = dim - 1;
    //go backwards is more efficient

    //lets move out inner calculations
    int ii, jj;
    pixel_sum sum;
    pixel current_pixel;
    int min_intensity = 766; // arbitrary value that is higher than maximum possible intensity, which is 255*3=765
    int max_intensity = -1; // arbitrary value that is lower than minimum possible intensity, which is 0
    pixel min_p, max_p;// more efficient
    pixel loop_pixel;

    // we check that this pixel is in the range of change (in the "red square as defined in ex5")
    int halfKernelSize = kernelSize >> 1; //=kernelSize/2
    int endRange = dim - halfKernelSize;

    //calculating lengths:
    int startLength = start + halfKernelSize;
    int length1 = startLength;
    int startIIAndJJ = start - halfKernelSize;
    int startII = startIIAndJJ;
    int startStartIndexPixel = startIIAndJJ * (dim + 1);
    for (i = start ; i >= 0; --i) {
        int length2 = startLength;
        int startJJ = startIIAndJJ;
        int startIndexPixel = startStartIndexPixel;
        for (j =  start ; j >= 0; --j) {
            if( i < halfKernelSize ||
                i >= endRange ||
                j < halfKernelSize ||
                j >= endRange){
                    dst[dstIndex] = *(src+ i*dim +j);
			        --dstIndex;
                    --length2;
                    --startJJ;
                    --startIndexPixel;
                    continue;
            }
            //so, we don't need min or max

            // func for this little thing isn't effcient:
            // initialize_pixel_sum(&sum);
            sum.red = 0;
            sum.green = 0;
            sum.blue = 0;

            //Optimize by calculate index efficiently in loop:
            int* kIndex = (int *) kernel;
            ii = startII;
            int indexPixel = startIndexPixel;
            //best to check here and calculate to filter in his loop for minimize the memory we read.
            if(!filter){
                for(; ii <= length1; ++ii) {
                    jj = startJJ;
                    int finalIndex = indexPixel;
                    for(; jj <= length2; ++jj) {
                        // apply kernel on pixel at [ii,jj]

                        // Optimization don't use func in this loop instead (many things in stack):
                        pixel p = src[finalIndex];
                        int weight = *(kIndex);
                        sum.red += ((int) p.red) * weight;
                        sum.green += ((int) p.green) * weight;
                        sum.blue += ((int) p.blue) * weight;

                        ++kIndex;
                        ++finalIndex;
                    }
                    indexPixel += dim;
                }
            } else {
                for(; ii <= length1; ++ii) {
                    jj = startJJ;
                    int finalIndex = indexPixel;
                    for(; jj <= length2; ++jj) {
                        // check if smaller than min or higher than max and update
                        loop_pixel = src[finalIndex];
                        int r = (int) loop_pixel.red;
                        int g = (int) loop_pixel.green;
                        int b = (int) loop_pixel.blue;
                        // apply kernel on pixel at [ii,jj]
                        // Optimization don't use func in this loop instead (many things in stack):
                        int weight = *(kIndex);
                        sum.red += r * weight;
                        sum.green += g * weight;
                        sum.blue += b * weight;
                        
                        ++kIndex;
                        
                        //calculating once this arg:
                        int loop_pixel_intensity = r + b + g;

                        if (loop_pixel_intensity <= min_intensity) {
                            min_intensity = loop_pixel_intensity;

                            //more efficient to remember the loc:
                            min_p = loop_pixel;
                        }

                        if (loop_pixel_intensity > max_intensity) {
                            max_intensity = loop_pixel_intensity;

                            //more efficient to remember the loc:
                            max_p = loop_pixel;
                        }

                        ++finalIndex;
                    }
                    indexPixel += dim;
                }

                // filter out min and max
                // func are not effecient:
                // sum_pixels_by_weight(&sum, min_p, -1);
                // sum_pixels_by_weight(&sum, max_p, -1);
                sum.red -= min_p.red + max_p.red;
                sum.green -= min_p.green + max_p.green;
                sum.blue -= min_p.blue + max_p.blue;

                min_intensity = 766; // arbitrary value that is higher than maximum possible intensity, which is 255*3=765
                max_intensity = -1; // arbitrary value that is lower than minimum possible intensity, which is 0
            }

            // this func is called a lot so we won't use as func - heavy on stack
            // assign_sum_to_pixel(&current_pixel, sum, kernelScale);
            if(kernelScale != 1) {  //often the scale is 1.
                // divide by kernel's weight
                sum.red /= kernelScale;
                sum.green /= kernelScale;
                sum.blue /= kernelScale;
            }

            // truncate each pixel's color values to match the range [0,255]
            // this func is called a lot we won't want she also calls funcs.
            if (sum.red < 0){sum.red = 0;} else if(sum.red > 255) {sum.red = 255;}
            current_pixel.red = (unsigned char) sum.red;
            if (sum.green < 0){sum.green = 0;} else if(sum.green > 255) {sum.green = 255;}
            current_pixel.green = (unsigned char) sum.green;
            if (sum.blue < 0){sum.blue = 0;} else if(sum.blue > 255) {sum.blue = 255;}
            current_pixel.blue = (unsigned char) sum.blue;

            dst[dstIndex] = current_pixel;
			--dstIndex;
            --length2;
            --startJJ;
            --startIndexPixel;
        }
        --length1;
        --startII;
        startStartIndexPixel -= dim;
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
			++indexMult3;
            pixels[index].green = image->data[indexMult3];
			++indexMult3;
            pixels[index].blue = image->data[indexMult3];
			++indexMult3;

			++index;
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
			++indexMult3;
            image->data[indexMult3] = pixels[index].green;
			++indexMult3;
            image->data[indexMult3] = pixels[index].blue;
			++indexMult3;

			++index;
        }
    }
}

void copyPixels(pixel *src, pixel *dst) {
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

    charsToPixels(image, backupOrg);

	// Optimization: we make sure smooth go over all the pixels instead (we changed smooth):
    // copyPixels(pixelsImg, backupOrg);

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


