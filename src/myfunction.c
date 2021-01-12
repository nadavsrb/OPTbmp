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
pixel_sum* getSumMatrix(int dim, pixel *src){
    int srcSize = dim*dim;
    int sumSize = srcSize-(dim<<2)+4;
    pixel_sum* sum = malloc(sumSize*sizeof(pixel_sum));
    pixel * srcKernelIndex = src + sumSize + (dim<<1) -7;
    pixel* srcKernelIndexIterNext = srcKernelIndex + dim;
    pixel* srcKernelIndexIterNextNext = srcKernelIndexIterNext + dim;
    pixel_sum* sumIndex = sum + sumSize -1;
    int startIAndJ = dim - 3;
    pixel* firstSrcIndex = src + srcSize -1 - (dim<<1);
    pixel* firstSrcIndexNext = firstSrcIndex + dim;
    pixel* firstSrcIndexNextNext = firstSrcIndexNext + dim;
    pixel* secondSrcIndex = firstSrcIndex - 1;
    pixel* secondSrcIndexNext = firstSrcIndexNext -1;
    pixel* secondSrcIndexNextNext = firstSrcIndexNextNext -1;
    pixel_sum colSum[3];
    pixel_sum* colSumIterNext = colSum + 1;
    pixel_sum* colSumIterNextNext = colSumIterNext + 1;
    pixel_sum* startColToFillIndex = colSum + 2;
    pixel_sum* endColToFillIndex = colSum + 3;

    int i =  startIAndJ;
    for(; i >= 0; --i){
        
        colSum->red = (int) firstSrcIndex->red;
        colSum->green = (int) firstSrcIndex->green;
        colSum->blue = (int) firstSrcIndex->blue;

        colSum->red += (int) firstSrcIndexNext->red;
        colSum->green += (int) firstSrcIndexNext->green;
        colSum->blue += (int) firstSrcIndexNext->blue;

        colSum->red += (int) firstSrcIndexNextNext->red;
        colSum->green += (int) firstSrcIndexNextNext->green;
        colSum->blue += (int) firstSrcIndexNextNext->blue;

        colSumIterNext->red = (int) secondSrcIndex->red;
        colSumIterNext->green = (int) secondSrcIndex->green;
        colSumIterNext->blue = (int) secondSrcIndex->blue;

        colSumIterNext->red += (int) secondSrcIndexNext->red;
        colSumIterNext->green += (int) secondSrcIndexNext->green;
        colSumIterNext->blue += (int) secondSrcIndexNext->blue;

        colSumIterNext->red += (int) secondSrcIndexNextNext->red;
        colSumIterNext->green += (int) secondSrcIndexNextNext->green;
        colSumIterNext->blue += (int) secondSrcIndexNextNext->blue;

        firstSrcIndex -= dim;
        firstSrcIndexNext -= dim;
        firstSrcIndexNextNext -= dim;
        secondSrcIndex = firstSrcIndex - 1;
        secondSrcIndexNext = firstSrcIndexNext - 1;
        secondSrcIndexNextNext = firstSrcIndexNextNext - 1;
        
        pixel_sum *colSumIndex = startColToFillIndex;
        int j = startIAndJ;
        for(; j>=0 ; --j){

            colSumIndex->red = (int) srcKernelIndex->red + (int) srcKernelIndexIterNext->red + (int) srcKernelIndexIterNextNext->red;
            colSumIndex->green = (int) srcKernelIndex->green + (int) srcKernelIndexIterNext->green + (int) srcKernelIndexIterNextNext->green;
            colSumIndex->blue = (int) srcKernelIndex->blue + (int) srcKernelIndexIterNext->blue + (int) srcKernelIndexIterNextNext->blue;

            ++colSumIndex;
            if(colSumIndex == endColToFillIndex) {colSumIndex = colSum;}

            sumIndex->red = colSum->red + colSumIterNext->red + colSumIterNextNext-> red;
            sumIndex->green = colSum->green + colSumIterNext->green + colSumIterNextNext-> green;
            sumIndex->blue = colSum->blue + colSumIterNext->blue + colSumIterNextNext-> blue;

            --sumIndex;
            --srcKernelIndex;
            --srcKernelIndexIterNext;
            --srcKernelIndexIterNextNext;
        }
        srcKernelIndex -= 2;
        srcKernelIndexIterNext -= 2;
        srcKernelIndexIterNextNext -= 2;
    }

    return sum;
}


pixel_sum* getFilteredSumMatrix(int dim, pixel *src){
    pixel_sum* sum = malloc((dim-2)*(dim-2)*sizeof(pixel_sum));
    pixel * srcKernelIndex = src + (dim -2) *dim -3;
    pixel_sum* sumIndex = sum + (dim-2)*(dim-2) -1;

    int endIndexKernel = dim*dim -1;
    int i =  dim - 3;
    for(; i >= 0; --i){

        pixel_sum colSum[3] = {0};
        int colToFillIndex =  2;
        int minIntensities[3] = {766, 766, 766};
        int minRaw[3];
        pixel_sum minPixels[3];
        int maxIntensities[3] = {-1, -1, -1};
        int maxRaw[3];
        pixel_sum maxPixels[3]; 
        pixel_sum pixelInSum;
        int a = 1;
        for(; a >=0; --a){
            int b = 2;
            for(; b >=0; --b){
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
        
        int j = dim - 3;
        for(; j>=0 ; --j){
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

            int k = 1;
            for(; k < 3; ++k){
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

            int s = 0;
            for (; s < 2; ++s){
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
    pixel* srcIndex = src+index;
    pixel* dstIndex = dst+index;
    int i = dim;
    int start = dim - 2;

    for(; i > 0; --i){
        *dstIndex = *srcIndex;
        --srcIndex;
        --dstIndex;
    }

    pixel_sum sum;
    pixel current_pixel;
    pixel_sum* sumMatIndex = sumMat + index - (dim<<2) +4;
    int maxkernelScale = 255 * kernelScale;
    for(i = start; i > 0; --i){
        *dstIndex = *srcIndex;
        --srcIndex;
        --dstIndex;

        int j = start;
        for(; j > 0; --j) {
        sum = *sumMatIndex;
        // this func is called a lot so we won't use as func - heavy on stack
        // assign_sum_to_pixel(&current_pixel, sum, kernelScale);

        // sum.red = srcIndex->red * 10 - sum.red;
        // sum.green = srcIndex->green * 10 - sum.green;
        // sum.blue = srcIndex->blue *10 - sum.blue;

        // truncate each pixel's color values to match the range [0,255]
        // this func is called a lot we won't want she also calls funcs.
        if (sum.red < 0){
            current_pixel.red = 0;
        } else if(sum.red > maxkernelScale) {
            current_pixel.red = 255;
        } else {
            current_pixel.red = (unsigned char) (sum.red / kernelScale);
        }

        if (sum.green < 0){
            current_pixel.green = 0;
        } else if(sum.green > maxkernelScale) {
            current_pixel.green = 255;
        } else {
            current_pixel.green = (unsigned char) (sum.green / kernelScale);
        }

        if (sum.blue < 0){
            current_pixel.blue = 0;
        } else if(sum.blue > maxkernelScale) {
            current_pixel.blue = 255;
        } else {
            current_pixel.blue = (unsigned char) (sum.blue / kernelScale);
        }
        

        *dstIndex = current_pixel;

        --sumMatIndex;
        --srcIndex;
        --dstIndex;
        }

        *dstIndex = *srcIndex;
        --srcIndex;
        --dstIndex;
    }

    for(i = dim; i > 0; --i){
         *dstIndex = *srcIndex;
        --srcIndex;
        --dstIndex;
    }
}



void sharppen(int dim, pixel *src, pixel *dst, pixel_sum* sumMat) {
    int index = dim * dim - 1;
    pixel* srcIndex = src+index;
    pixel* dstIndex = dst+index;
    int i = dim;
    int start = dim - 2;

    for(; i > 0; --i){
        *dstIndex = *srcIndex;
        --srcIndex;
        --dstIndex;
    }

    pixel_sum sum;
    pixel current_pixel;
    pixel_sum* sumMatIndex = sumMat + start*start -1;
    for(i = start; i > 0; --i){
        *dstIndex = *srcIndex;
        --srcIndex;
        --dstIndex;

        int j = start;
        for(; j > 0; --j) {
        sum = *sumMatIndex;
        // this func is called a lot so we won't use as func - heavy on stack
        // assign_sum_to_pixel(&current_pixel, sum, kernelScale);
        int r = srcIndex->red;
        int b = srcIndex->blue;
        int g = srcIndex->green;

        sum.red = (r<<3) + (r<<1) - sum.red;
        sum.green = (g<<3) + (g<<1) - sum.green;
        sum.blue = (b<<3) + (b<<1) - sum.blue;

        // sum.red = srcIndex->red * 10 - sum.red;
        // sum.green = srcIndex->green * 10 - sum.green;
        // sum.blue = srcIndex->blue *10 - sum.blue;

        // truncate each pixel's color values to match the range [0,255]
        // this func is called a lot we won't want she also calls funcs.
        if (sum.red < 0){sum.red = 0;} else if(sum.red > 255) {sum.red = 255;}
        current_pixel.red = (unsigned char) sum.red;
        if (sum.green < 0){sum.green = 0;} else if(sum.green > 255) {sum.green = 255;}
        current_pixel.green = (unsigned char) sum.green;
        if (sum.blue < 0){sum.blue = 0;} else if(sum.blue > 255) {sum.blue = 255;}
        current_pixel.blue = (unsigned char) sum.blue;

        *dstIndex = current_pixel;

        --sumMatIndex;
        --srcIndex;
        --dstIndex;
        }

        *dstIndex = *srcIndex;
        --srcIndex;
        --dstIndex;
    }

    for(i = dim; i > 0; --i){
         *dstIndex = *srcIndex;
        --srcIndex;
        --dstIndex;
    }

}

void charsToPixels(Image *charsImg, pixel* pixels) {
//calculating index efficiently
	unsigned char* indexImage = charsImg->data;
    pixel* indexPixels = pixels;

    int row;
    for (row = 0 ; row < m ; ++row) {
		int col = 0;
        for (; col < n ; ++col) {

            indexPixels->red = *indexImage;
			++indexImage;
            indexPixels->green = *indexImage;
			++indexImage;
            indexPixels->blue = *indexImage;
			++indexImage;

			++indexPixels;
        }
    }
}

void pixelsToChars(pixel* pixels, Image *charsImg) {
	//calculating index efficiently
	unsigned char* indexImage = charsImg->data;
    pixel* indexPixels = pixels;

    int row;
    for (row = 0 ; row < m ; ++row) {
		int col = 0;
        for (; col < n ; ++col) {

            *indexImage = indexPixels->red;
			++indexImage;
            *indexImage = indexPixels->green;
			++indexImage;
            *indexImage = indexPixels->blue;
			++indexImage;

			++indexPixels;
        }
    }
}

void myfunction(Image *image, char* srcImgpName, char* blurRsltImgName, char* sharpRsltImgName, char* filteredBlurRsltImgName, char* filteredSharpRsltImgName, char flag) {
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
