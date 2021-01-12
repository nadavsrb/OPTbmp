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

// /*
// * Apply the kernel over each pixel.
// * doesn't Ignore pixels where the kernel exceeds bounds-> applyKernel ignores.
// */
// void smooth(int dim, pixel *src, pixel *dst, int kernelSize, int kernel[kernelSize][kernelSize], int kernelScale, bool filter) {

//     int i, j;
// 	//calculating index efficiently
// 	int Index = dim * dim - 1;
//     pixel* dstIndex = dst + Index;
//     pixel* srcIndex = src + Index;
//     int start = dim - 1;
//     //go backwards is more efficient

//     //lets move out inner calculations
//     int ii, jj;
//     pixel_sum sum;
//     pixel current_pixel;
//     int min_intensity = 766; // arbitrary value that is higher than maximum possible intensity, which is 255*3=765
//     int max_intensity = -1; // arbitrary value that is lower than minimum possible intensity, which is 0
//     pixel min_p, max_p;// more efficient
//     pixel loop_pixel;

//     // we check that this pixel is in the range of change (in the "red square as defined in ex5")
//     int halfKernelSize = kernelSize >> 1; //=kernelSize/2
//     int endRange = dim - halfKernelSize;
//     int sizeRedRaw = endRange - halfKernelSize;


//     //calculating lengths:
//     int length1 = start;
//     int startIIAndJJ = dim - kernelSize;
//     int startII = startIIAndJJ;
//     pixel* startStartIndexPixel = srcIndex - (dim + 1)*(kernelSize-1);
    
//     //LAST ROWS ISN'T IN THE RED CIRCLE:
//     for(i = start; i >= endRange; --i){
//         for (j =  start ; j >= 0; --j)  {
//             *dstIndex = *srcIndex;
//             --dstIndex;
//             --srcIndex;
//         }
//     }

//     for (; i >= halfKernelSize; --i) {
//         int length2 = start;
//         int startJJ = startIIAndJJ;
//         pixel* startIndexPixel = startStartIndexPixel;

//         for (j =  start ; j >= endRange; --j) {
//             *dstIndex = *srcIndex;
//             --dstIndex;
//             --srcIndex;
//         }

//         for (; j >= halfKernelSize; --j) {
//             //so, we don't need min or max

//             // func for this little thing isn't effcient:
//             // initialize_pixel_sum(&sum);
//             sum.red = 0;
//             sum.green = 0;
//             sum.blue = 0;

//             //Optimize by calculate index efficiently in loop:
//             int* kIndex = (int *) kernel;
//             ii = startII;
//             pixel* indexPixel = startIndexPixel;
//             //best to check here and calculate to filter in his loop for minimize the memory we read.
//             if(!filter){
//                 for(; ii <= length1; ++ii) {
//                     jj = startJJ;
//                     pixel* finalIndex = indexPixel;
//                     for(; jj <= length2; ++jj) {
//                         // apply kernel on pixel at [ii,jj]

//                         // Optimization don't use func in this loop instead (many things in stack):
//                         pixel p = *finalIndex;
//                         int weight = *(kIndex);
//                         sum.red += ((int) p.red) * weight;
//                         sum.green += ((int) p.green) * weight;
//                         sum.blue += ((int) p.blue) * weight;

//                         ++kIndex;
//                         ++finalIndex;
//                     }
//                     indexPixel += dim;
//                 }
//             } else {
//                 for(; ii <= length1; ++ii) {
//                     jj = startJJ;
//                     pixel* finalIndex = indexPixel;
//                     for(; jj <= length2; ++jj) {
//                         // check if smaller than min or higher than max and update
//                         loop_pixel = *finalIndex;
//                         int r = (int) loop_pixel.red;
//                         int g = (int) loop_pixel.green;
//                         int b = (int) loop_pixel.blue;
//                         // apply kernel on pixel at [ii,jj]
//                         // Optimization don't use func in this loop instead (many things in stack):
//                         int weight = *(kIndex);
//                         sum.red += r * weight;
//                         sum.green += g * weight;
//                         sum.blue += b * weight;
                        
//                         ++kIndex;
                        
//                         //calculating once this arg:
//                         int loop_pixel_intensity = r + b + g;

//                         if (loop_pixel_intensity <= min_intensity) {
//                             min_intensity = loop_pixel_intensity;

//                             //more efficient to remember the loc:
//                             min_p = loop_pixel;
//                         }

//                         if (loop_pixel_intensity > max_intensity) {
//                             max_intensity = loop_pixel_intensity;

//                             //more efficient to remember the loc:
//                             max_p = loop_pixel;
//                         }

//                         ++finalIndex;
//                     }
//                     indexPixel += dim;
//                 }

//                 // filter out min and max
//                 // func are not effecient:
//                 // sum_pixels_by_weight(&sum, min_p, -1);
//                 // sum_pixels_by_weight(&sum, max_p, -1);
//                 sum.red -= min_p.red + max_p.red;
//                 sum.green -= min_p.green + max_p.green;
//                 sum.blue -= min_p.blue + max_p.blue;

//                 min_intensity = 766; // arbitrary value that is higher than maximum possible intensity, which is 255*3=765
//                 max_intensity = -1; // arbitrary value that is lower than minimum possible intensity, which is 0
//             }

//             // this func is called a lot so we won't use as func - heavy on stack
//             // assign_sum_to_pixel(&current_pixel, sum, kernelScale);
//             if(kernelScale != 1) {  //often the scale is 1.
//                 // divide by kernel's weight
//                 sum.red /= kernelScale;
//                 sum.green /= kernelScale;
//                 sum.blue /= kernelScale;
//             }

//             // truncate each pixel's color values to match the range [0,255]
//             // this func is called a lot we won't want she also calls funcs.
//             if (sum.red < 0){sum.red = 0;} else if(sum.red > 255) {sum.red = 255;}
//             current_pixel.red = (unsigned char) sum.red;
//             if (sum.green < 0){sum.green = 0;} else if(sum.green > 255) {sum.green = 255;}
//             current_pixel.green = (unsigned char) sum.green;
//             if (sum.blue < 0){sum.blue = 0;} else if(sum.blue > 255) {sum.blue = 255;}
//             current_pixel.blue = (unsigned char) sum.blue;

//             *dstIndex = current_pixel;
// 			--dstIndex;
//             --length2;
//             --startJJ;
//             --startIndexPixel;
//         }
//         srcIndex -= sizeRedRaw;

//         for (; j >= 0; --j) {
//             *dstIndex = *srcIndex;
//             --dstIndex;
//             --srcIndex;
//         }

//         --length1;
//         --startII;
//         startStartIndexPixel -= dim;
//     }

//     //FIRST ROWS ISN'T IN THE RED CIRCLE:
//     for(; i >= 0; --i){
//         for (j =  start ; j >= 0; --j)  {
//             *dstIndex = *srcIndex;
//             --dstIndex;
//             --srcIndex;
//         }
//     }
// }

void smoothNoFilter(int dim, pixel *src, pixel *dst, int kernelSize, int kernel[kernelSize][kernelSize], int kernelScale) {

    int i, j;
	//calculating index efficiently
	int Index = dim * dim - 1;
    pixel* dstIndex = dst + Index;
    pixel* srcIndex = src + Index;
    int start = dim - 1;
    //go backwards is more efficient

    //lets move out inner calculations
    int ii, jj;
    pixel_sum sum;
    pixel current_pixel;
    pixel loop_pixel;

    // we check that this pixel is in the range of change (in the "red square as defined in ex5")
    int halfKernelSize = kernelSize >> 1; //=kernelSize/2
    int endRange = dim - halfKernelSize;
    int sizeRedRaw = endRange - halfKernelSize;


    //calculating lengths:
    int length1 = start;
    int startIIAndJJ = dim - kernelSize;
    int startII = startIIAndJJ;
    pixel* startStartIndexPixel = srcIndex - (dim + 1)*(kernelSize-1);
    
    //LAST ROWS ISN'T IN THE RED CIRCLE:
    for(i = start; i >= endRange; --i){
        for (j =  start ; j >= 0; --j)  {
            *dstIndex = *srcIndex;
            --dstIndex;
            --srcIndex;
        }
    }

    for (; i >= halfKernelSize; --i) {
        int length2 = start;
        int startJJ = startIIAndJJ;
        pixel* startIndexPixel = startStartIndexPixel;

        for (j =  start ; j >= endRange; --j) {
            *dstIndex = *srcIndex;
            --dstIndex;
            --srcIndex;
        }

        for (; j >= halfKernelSize; --j) {
            //so, we don't need min or max

            // func for this little thing isn't effcient:
            // initialize_pixel_sum(&sum);
            sum.red = 0;
            sum.green = 0;
            sum.blue = 0;

            //Optimize by calculate index efficiently in loop:
            int* kIndex = (int *) kernel;
            ii = startII;
            pixel* indexPixel = startIndexPixel;

            for(; ii <= length1; ++ii) {
                jj = startJJ;
                pixel* finalIndex = indexPixel;
                for(; jj <= length2; ++jj) {
                    // apply kernel on pixel at [ii,jj]

                    // Optimization don't use func in this loop instead (many things in stack):
                    pixel p = *finalIndex;
                    int weight = *(kIndex);
                    sum.red += ((int) p.red) * weight;
                    sum.green += ((int) p.green) * weight;
                    sum.blue += ((int) p.blue) * weight;

                    ++kIndex;
                    ++finalIndex;
                }
                indexPixel += dim;
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

            *dstIndex = current_pixel;
			--dstIndex;
            --length2;
            --startJJ;
            --startIndexPixel;
        }
        srcIndex -= sizeRedRaw;

        for (; j >= 0; --j) {
            *dstIndex = *srcIndex;
            --dstIndex;
            --srcIndex;
        }

        --length1;
        --startII;
        startStartIndexPixel -= dim;
    }

    //FIRST ROWS ISN'T IN THE RED CIRCLE:
    for(; i >= 0; --i){
        for (j =  start ; j >= 0; --j)  {
            *dstIndex = *srcIndex;
            --dstIndex;
            --srcIndex;
        }
    }
}

void smoothWithFilter(int dim, pixel *src, pixel *dst, int kernelSize, int kernel[kernelSize][kernelSize], int kernelScale) {

    int i, j;
	//calculating index efficiently
	int Index = dim * dim - 1;
    pixel* dstIndex = dst + Index;
    pixel* srcIndex = src + Index;
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
    int sizeRedRaw = endRange - halfKernelSize;


    //calculating lengths:
    int length1 = start;
    int startIIAndJJ = dim - kernelSize;
    int startII = startIIAndJJ;
    pixel* startStartIndexPixel = srcIndex - (dim + 1)*(kernelSize-1);
    
    //LAST ROWS ISN'T IN THE RED CIRCLE:
    for(i = start; i >= endRange; --i){
        for (j =  start ; j >= 0; --j)  {
            *dstIndex = *srcIndex;
            --dstIndex;
            --srcIndex;
        }
    }

    for (; i >= halfKernelSize; --i) {
        int length2 = start;
        int startJJ = startIIAndJJ;
        pixel* startIndexPixel = startStartIndexPixel;

        for (j =  start ; j >= endRange; --j) {
            *dstIndex = *srcIndex;
            --dstIndex;
            --srcIndex;
        }

        for (; j >= halfKernelSize; --j) {
            //so, we don't need min or max

            // func for this little thing isn't effcient:
            // initialize_pixel_sum(&sum);
            sum.red = 0;
            sum.green = 0;
            sum.blue = 0;

            //Optimize by calculate index efficiently in loop:
            int* kIndex = (int *) kernel;
            ii = startII;
            pixel* indexPixel = startIndexPixel;
            //best to check here and calculate to filter in his loop for minimize the memory we read.
            for(; ii <= length1; ++ii) {
                jj = startJJ;
                pixel* finalIndex = indexPixel;
                for(; jj <= length2; ++jj) {
                    // check if smaller than min or higher than max and update
                    loop_pixel = *finalIndex;
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

            *dstIndex = current_pixel;
			--dstIndex;
            --length2;
            --startJJ;
            --startIndexPixel;
        }
        srcIndex -= sizeRedRaw;

        for (; j >= 0; --j) {
            *dstIndex = *srcIndex;
            --dstIndex;
            --srcIndex;
        }

        --length1;
        --startII;
        startStartIndexPixel -= dim;
    }

    //FIRST ROWS ISN'T IN THE RED CIRCLE:
    for(; i >= 0; --i){
        for (j =  start ; j >= 0; --j)  {
            *dstIndex = *srcIndex;
            --dstIndex;
            --srcIndex;
        }
    }
}


pixel_sum* getSumMatrix(int dim, pixel *src){
    pixel_sum* sum = malloc((dim-2)*(dim-2)*sizeof(pixel_sum));
    pixel * srcKernelIndex = src + (dim -2) *dim -3;
    pixel_sum* sumIndex = sum + (dim-2)*(dim-2) -1;

    int endIndexKernel = dim*dim -1;
    for(int i =  dim - 3; i >= 0; --i){

        pixel_sum colSum[3] = {0};
        int colToFillIndex =  2;
        for(int a = 1; a >=0; --a){
            for(int b = 2; b >=0; --b){
                colSum[a].red += (int) src[endIndexKernel - a - b * dim].red;
                colSum[a].green += (int) src[endIndexKernel - a - b * dim].green;
                colSum[a].blue += (int) src[endIndexKernel - a - b * dim].blue;
            }
        }
        endIndexKernel -= dim;

        for(int j = dim - 3; j>=0 ; --j){
            colSum[colToFillIndex].red = (int) srcKernelIndex->red;
            colSum[colToFillIndex].green = (int) srcKernelIndex->green;
            colSum[colToFillIndex].blue = (int) srcKernelIndex->blue;

            for(int k = 1; k < 3; ++k){
                colSum[colToFillIndex].red += (int) (srcKernelIndex + dim*k)->red;
                colSum[colToFillIndex].green += (int) (srcKernelIndex + dim*k)->green;
                colSum[colToFillIndex].blue += (int) (srcKernelIndex + dim*k)->blue;
            }

            ++colToFillIndex;
            if(colToFillIndex == 3) {colToFillIndex = 0;}

            sumIndex->red = colSum[0].red + colSum[1].red + colSum[2].red;
            sumIndex->green = colSum[0].green + colSum[1].green + colSum[2].green;
            sumIndex->blue = colSum[0].blue + colSum[1].blue + colSum[2].blue;

            --sumIndex;
            --srcKernelIndex;
        }
        srcKernelIndex-=2;
    }

    return sum;
}


pixel_sum* getFilteredSumMatrix(int dim, pixel *src){
    pixel_sum* sum = malloc((dim-2)*(dim-2)*sizeof(pixel_sum));
    pixel * srcKernelIndex = src + (dim -2) *dim -3;
    pixel_sum* sumIndex = sum + (dim-2)*(dim-2) -1;

    int endIndexKernel = dim*dim -1;
    for(int i =  dim - 3; i >= 0; --i){

        pixel_sum colSum[3] = {0};
        int colToFillIndex =  2;
        int minIntensities[3] = {766, 766, 766};
        int minRaw[3];
        pixel_sum minPixels[3];
        int maxIntensities[3] = {-1, -1, -1};
        int maxRaw[3];
        pixel_sum maxPixels[3]; 
        pixel_sum pixelInSum;
        for(int a = 1; a >=0; --a){
            for(int b = 2; b >=0; --b){
                pixelInSum.red = (int) src[endIndexKernel - a - b * dim].red;
                pixelInSum.green = (int) src[endIndexKernel - a - b * dim].green;
                pixelInSum.blue = (int) src[endIndexKernel - a - b * dim].blue;

                colSum[a].red += pixelInSum.red;
                colSum[a].green += pixelInSum.green;
                colSum[a].blue += pixelInSum.blue;

                int intensity = pixelInSum.red + pixelInSum.green + pixelInSum.blue;

                if(intensity <= minIntensities[a]){
                    minIntensities[a] = intensity;
                    minRaw[a] = 2 - b;
                    minPixels[a] = pixelInSum;
                }

                if(intensity > maxIntensities[a]) {
                    maxIntensities[a] = intensity;
                    maxRaw[a] = 2 - b;
                    maxPixels[a] = pixelInSum;
                }

            }
        }
        endIndexKernel -= dim;

        for(int j = dim - 3; j>=0 ; --j){
            pixelInSum.red = (int) srcKernelIndex->red;
            pixelInSum.green = (int) srcKernelIndex->green;
            pixelInSum.blue = (int) srcKernelIndex->blue;

            colSum[colToFillIndex].red = pixelInSum.red;
            colSum[colToFillIndex].green = pixelInSum.green;
            colSum[colToFillIndex].blue = pixelInSum.blue;

            int intensity = pixelInSum.red + pixelInSum.green + pixelInSum.blue;

            if(intensity <= minIntensities[colToFillIndex]){
                minIntensities[colToFillIndex] = intensity;
                minPixels[colToFillIndex] = pixelInSum;
                minRaw[colToFillIndex] = 0;
            }

            if(intensity > maxIntensities[colToFillIndex]) {
                maxIntensities[colToFillIndex] = intensity;
                maxPixels[colToFillIndex] = pixelInSum;
                maxRaw[colToFillIndex] = 0;
            }

            for(int k = 1; k < 3; ++k){
                pixelInSum.red = (int) (srcKernelIndex + dim*k)->red;
                pixelInSum.green = (int) (srcKernelIndex + dim*k)->green;
                pixelInSum.blue = (int) (srcKernelIndex + dim*k)->blue;

                colSum[colToFillIndex].red += pixelInSum.red;
                colSum[colToFillIndex].green += pixelInSum.green;
                colSum[colToFillIndex].blue += pixelInSum.blue;

                intensity = pixelInSum.red + pixelInSum.green + pixelInSum.blue;

                if(intensity <= minIntensities[colToFillIndex]){
                    minIntensities[colToFillIndex] = intensity;
                    minPixels[colToFillIndex] = pixelInSum;
                    minRaw[colToFillIndex] = k;
                }
                

                if(intensity > maxIntensities[colToFillIndex]) {
                    maxIntensities[colToFillIndex] = intensity;
                    maxPixels[colToFillIndex] = pixelInSum;
                    maxRaw[colToFillIndex] = k;
                }
            }

            int minIndex  = colToFillIndex;
            int maxIndex  = colToFillIndex;
            int checkIndex = colToFillIndex - 1;
            for (int s = 0; s < 2; ++s){
                if(checkIndex == -1 ) {checkIndex = 2;}
                
                if(minIntensities[minIndex] > minIntensities[checkIndex] ||
                (minIntensities[minIndex] == minIntensities[checkIndex] &&
                minRaw[checkIndex] >= minRaw[minIndex])) {
                    minIndex = checkIndex;
                }

                if(maxIntensities[maxIndex] < maxIntensities[checkIndex] ||
                (maxIntensities[maxIndex] == maxIntensities[checkIndex] &&
                maxRaw[checkIndex] < maxRaw[maxIndex])) {
                    maxIndex = checkIndex;
                }

                --checkIndex;
            }

            sumIndex->red = colSum[0].red + colSum[1].red + colSum[2].red - minPixels[minIndex].red - maxPixels[maxIndex].red;
            sumIndex->green = colSum[0].green + colSum[1].green + colSum[2].green - minPixels[minIndex].green - maxPixels[maxIndex].green;
            sumIndex->blue = colSum[0].blue + colSum[1].blue + colSum[2].blue - minPixels[minIndex].blue - maxPixels[maxIndex].blue;

            ++colToFillIndex;
            if(colToFillIndex == 3) {colToFillIndex = 0;}
            minIntensities[colToFillIndex] = 766;
            maxIntensities[colToFillIndex] = -1;

            --sumIndex;
            --srcKernelIndex;
        }
        srcKernelIndex-=2;
    }
    return sum;
}

void blur(int dim, pixel *src, pixel *dst, pixel_sum* sumMat, int kernelScale) {
    int index = dim * dim - 1;
    for(int i = dim; i > 0; --i){
        dst[index] = src[index];
        --index;
    }
    pixel_sum sum;
    pixel current_pixel;
    pixel_sum* sumMatIndex = sumMat + (dim-2)*(dim-2) -1;
    for(int i = dim-2; i > 0; --i){
        dst[index] = src[index];
        --index;

        for(int j = dim-2; j > 0; --j) {
        sum = *sumMatIndex;
        // this func is called a lot so we won't use as func - heavy on stack
        // assign_sum_to_pixel(&current_pixel, sum, kernelScale);
        sum.red /= kernelScale;
        sum.green /= kernelScale;
        sum.blue /= kernelScale;

        // truncate each pixel's color values to match the range [0,255]
        // this func is called a lot we won't want she also calls funcs.
        if (sum.red < 0){sum.red = 0;} else if(sum.red > 255) {sum.red = 255;}
        current_pixel.red = (unsigned char) sum.red;
        if (sum.green < 0){sum.green = 0;} else if(sum.green > 255) {sum.green = 255;}
        current_pixel.green = (unsigned char) sum.green;
        if (sum.blue < 0){sum.blue = 0;} else if(sum.blue > 255) {sum.blue = 255;}
        current_pixel.blue = (unsigned char) sum.blue;

        dst[index] = current_pixel;

        --sumMatIndex;
        --index;
        }

        dst[index] = src[index];
        --index;
    }
    for(int i = dim; i > 0; --i){
        dst[index] = src[index];
        --index;
    }
}



void sharppen(int dim, pixel *src, pixel *dst, pixel_sum* sumMat) {
    int index = dim * dim - 1;

    for(int i = dim; i > 0; --i){
        dst[index] = src[index];
        --index;
    }

    pixel_sum sum;
    pixel current_pixel;
    pixel_sum* sumMatIndex = sumMat + (dim-2)*(dim-2) -1;
    for(int i = dim-2; i > 0; --i){
        dst[index] = src[index];
        --index;

        for(int j = dim-2; j > 0; --j) {
        sum = *sumMatIndex;
        // this func is called a lot so we won't use as func - heavy on stack
        // assign_sum_to_pixel(&current_pixel, sum, kernelScale);
        sum.red = 10 * src[index].red - sum.red;
        sum.green = 10 * src[index].green - sum.green;
        sum.blue = 10 * src[index].blue - sum.blue;

        // truncate each pixel's color values to match the range [0,255]
        // this func is called a lot we won't want she also calls funcs.
        if (sum.red < 0){sum.red = 0;} else if(sum.red > 255) {sum.red = 255;}
        current_pixel.red = (unsigned char) sum.red;
        if (sum.green < 0){sum.green = 0;} else if(sum.green > 255) {sum.green = 255;}
        current_pixel.green = (unsigned char) sum.green;
        if (sum.blue < 0){sum.blue = 0;} else if(sum.blue > 255) {sum.blue = 255;}
        current_pixel.blue = (unsigned char) sum.blue;

        dst[index] = current_pixel;

        --sumMatIndex;
        --index;
        }

        dst[index] = src[index];
        --index;
    }

    for(int i = dim; i > 0; --i){
        dst[index] = src[index];
        --index;
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
    if (filter) {
        smoothWithFilter(m, backupOrg, pixelsImg, kernelSize, kernel, kernelScale);
    } else {
        smoothNoFilter(m, backupOrg, pixelsImg, kernelSize, kernel, kernelScale);
    }

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

    size_t size = m*n*sizeof(pixel);
    pixel* pixelsImg = malloc(size);
    pixel* backupOrg = malloc(size);
    //caculating length once:
    charsToPixels(image, backupOrg);
    if (flag == '1') {  
        pixel_sum* sum = getSumMatrix(m, backupOrg);
        // Optimization: we make sure smooth go over all the pixels instead (we changed smooth):
        // copyPixels(pixelsImg, backupOrg);
        blur(m, backupOrg, pixelsImg, sum, 9);

        pixelsToChars(pixelsImg, image);

        // write result image to file
        writeBMP(image, srcImgpName, blurRsltImgName);

        sum = getSumMatrix(m, pixelsImg);

        // sharpen the resulting image
        sharppen(m, pixelsImg, backupOrg, sum);
        
        pixelsToChars(backupOrg, image);
        
        // write result image to file
        writeBMP(image, srcImgpName, sharpRsltImgName); 

        free(sum);
        free(pixelsImg);
        free(backupOrg);
    } else {
        pixel_sum* sum = getFilteredSumMatrix(m, backupOrg);
        // Optimization: we make sure smooth go over all the pixels instead (we changed smooth):
        // copyPixels(pixelsImg, backupOrg);
        blur(m, backupOrg, pixelsImg, sum, 7);
        
        pixelsToChars(pixelsImg, image);

        // write result image to file
        writeBMP(image, srcImgpName, filteredBlurRsltImgName);

        sum = getSumMatrix(m, pixelsImg);

        // sharpen the resulting image
        sharppen(m, pixelsImg, backupOrg, sum);
        
        pixelsToChars(backupOrg, image);

        // write result image to file
        writeBMP(image, srcImgpName, filteredSharpRsltImgName); 

        free(sum);
        free(pixelsImg);
        free(backupOrg);
    }
}
