#pragma GCC optimize ("O3")
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

void getBlurNoFilter(int dim, pixel *src, pixel *imagePixels, Image *charsImg, int kernelScale){
    int srcSize = dim*dim;
    int index = srcSize - 1;
    int dimMult2 = dim<<1;
    int dimMult3 = dimMult2 + dim;
    int sumSize = srcSize-(dimMult2<<1)+4;
    int maxkernelScale = 255 * kernelScale;
    int startIAndJ = dim - 2;
    int startIAndJLoopUnrolling = startIAndJ - 1;
    unsigned char* indexImage = charsImg->data + (srcSize<<1) + srcSize -dimMult3 - 1;
    pixel* srcKernelIndex = src + sumSize + dimMult2 -7;
    pixel* srcKernelIndexIterNext = srcKernelIndex + dim;
    pixel* srcKernelIndexIterNextNext = srcKernelIndexIterNext + dim;
    pixel* imagePixelsIndex = imagePixels + index - dim;
    pixel* firstSrcIndex = src + srcSize -1 - dimMult2;
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
    // row 2:
    pixel_sum colSum2[3];
    pixel_sum* colSumIterNext2 = colSum2 + 1;
    pixel_sum* colSumIterNextNext2 = colSumIterNext2 + 1;
    pixel_sum* startColToFillIndex2 = colSum2 + 2;
    pixel* firstSrcIndexBefore = firstSrcIndex - dim;
    pixel* secondSrcIndexBefore = secondSrcIndex - dim;
    pixel* srcKernelIndexIterBefore = srcKernelIndex - dim;
    unsigned char* indexImage2 = indexImage -dimMult3;
    pixel* imagePixelsIndex2 = imagePixelsIndex - dim;

    int i = startIAndJLoopUnrolling;
    for(; i > 0; i -= 2){
        imagePixelsIndex->blue = *indexImage;
        --indexImage;
        imagePixelsIndex->green = *indexImage;
        --indexImage;
        imagePixelsIndex->red = *indexImage;
        --indexImage;

        --imagePixelsIndex;

        imagePixelsIndex2->blue = *indexImage2;
        --indexImage2;
        imagePixelsIndex2->green = *indexImage2;
        --indexImage2;
        imagePixelsIndex2->red = *indexImage2;
        --indexImage2;

        --imagePixelsIndex2;
        
        int helpSumR = (int) firstSrcIndex->red + (int) firstSrcIndexNext->red;
        int helpSumG = (int) firstSrcIndex->green + (int) firstSrcIndexNext->green;
        int helpSumB = (int) firstSrcIndex->blue + (int) firstSrcIndexNext->blue;
    
        colSum->red = helpSumR + (int) firstSrcIndexNextNext->red;
        colSum->green = helpSumG + (int) firstSrcIndexNextNext->green;
        colSum->blue = helpSumB + (int) firstSrcIndexNextNext->blue;

        colSum2->red = (int) firstSrcIndexBefore->red + helpSumR;
        colSum2->green = (int) firstSrcIndexBefore->green + helpSumG;
        colSum2->blue = (int) firstSrcIndexBefore->blue + helpSumB;

        helpSumR = (int) secondSrcIndex->red + (int) secondSrcIndexNext->red;
        helpSumG = (int) secondSrcIndex->green + (int) secondSrcIndexNext->green;
        helpSumB = (int) secondSrcIndex->blue + (int) secondSrcIndexNext->blue;

        colSumIterNext->red = helpSumR + (int) secondSrcIndexNextNext->red;
        colSumIterNext->green = helpSumG + (int) secondSrcIndexNextNext->green;
        colSumIterNext->blue = helpSumB + (int) secondSrcIndexNextNext->blue;

        colSumIterNext2->red = (int) secondSrcIndexBefore->red + helpSumR;
        colSumIterNext2->green = (int) secondSrcIndexBefore->green + helpSumG;
        colSumIterNext2->blue = (int) secondSrcIndexBefore->blue + helpSumB;

        firstSrcIndexNext = firstSrcIndexBefore;
        firstSrcIndexNextNext = firstSrcIndex;
        firstSrcIndexBefore -= dim;
        firstSrcIndex = firstSrcIndexBefore;
        firstSrcIndexBefore -= dim;

        secondSrcIndexNext = secondSrcIndexBefore;
        secondSrcIndexNextNext = secondSrcIndex;
        secondSrcIndex = firstSrcIndex - 1;
        secondSrcIndexBefore = firstSrcIndexBefore - 1;
        
        pixel_sum *colSumIndex = startColToFillIndex;
        pixel_sum *colSumIndex2 = startColToFillIndex2;
        int j = startIAndJ;
        for(; j>0 ; --j){
            helpSumR = (int) srcKernelIndex->red + (int) srcKernelIndexIterNext->red;
            helpSumG = (int) srcKernelIndex->green + (int) srcKernelIndexIterNext->green;
            helpSumB = (int) srcKernelIndex->blue + (int) srcKernelIndexIterNext->blue;

            colSumIndex->red = helpSumR + (int) srcKernelIndexIterNextNext->red;
            colSumIndex->green = helpSumG + (int) srcKernelIndexIterNextNext->green;
            colSumIndex->blue = helpSumB + (int) srcKernelIndexIterNextNext->blue;

            colSumIndex2->red = (int) srcKernelIndexIterBefore->red + helpSumR;
            colSumIndex2->green = (int) srcKernelIndexIterBefore->green + helpSumG;
            colSumIndex2->blue = (int) srcKernelIndexIterBefore->blue + helpSumB;

            ++colSumIndex;
            ++colSumIndex2;
            if(colSumIndex == endColToFillIndex) {
                colSumIndex = colSum;
                colSumIndex2 = colSum2;
            }

            pixel_sum sum = {colSum->red + colSumIterNext->red + colSumIterNextNext-> red,
            colSum->green + colSumIterNext->green + colSumIterNextNext-> green,
            colSum->blue + colSumIterNext->blue + colSumIterNextNext-> blue};

            pixel_sum sum2 = {colSum2->red + colSumIterNext2->red + colSumIterNextNext2-> red,
            colSum2->green + colSumIterNext2->green + colSumIterNextNext2-> green,
            colSum2->blue + colSumIterNext2->blue + colSumIterNextNext2-> blue};

            if (sum.blue < 0){
                *indexImage = imagePixelsIndex->blue = 0;
            } else if(sum.blue > maxkernelScale) {
                *indexImage = imagePixelsIndex->blue = 255;
            } else {
                *indexImage = imagePixelsIndex->blue = (unsigned char) (sum.blue / kernelScale);
            }

            --indexImage;

            if (sum.green < 0){
                *indexImage = imagePixelsIndex->green = 0;
            } else if(sum.green > maxkernelScale) {
                *indexImage = imagePixelsIndex->green = 255;
            } else {
                *indexImage = imagePixelsIndex->green = (unsigned char) (sum.green / kernelScale);
            }

            --indexImage;

            if (sum.red < 0){
                *indexImage = imagePixelsIndex->red = 0;
            } else if(sum.red > maxkernelScale) {
                *indexImage = imagePixelsIndex->red = 255;
            } else {
                *indexImage = imagePixelsIndex->red = (unsigned char) (sum.red / kernelScale);
            }

            --indexImage;

            --imagePixelsIndex;

            if (sum2.blue < 0){
                *indexImage2 = imagePixelsIndex2->blue = 0;
            } else if(sum2.blue > maxkernelScale) {
                *indexImage2 = imagePixelsIndex2->blue = 255;
            } else {
                *indexImage2 = imagePixelsIndex2->blue = (unsigned char) (sum2.blue / kernelScale);
            }

            --indexImage2;

            if (sum2.green < 0){
                *indexImage2 = imagePixelsIndex2->green = 0;
            } else if(sum2.green > maxkernelScale) {
                *indexImage2 = imagePixelsIndex2->green = 255;
            } else {
                *indexImage2 = imagePixelsIndex2->green = (unsigned char) (sum2.green / kernelScale);
            }

            --indexImage2;

            if (sum2.red < 0){
                *indexImage2 = imagePixelsIndex2->red = 0;
            } else if(sum2.red > maxkernelScale) {
                *indexImage2 = imagePixelsIndex2->red = 255;
            } else {
                *indexImage2 = imagePixelsIndex2->red = (unsigned char) (sum2.red / kernelScale);
            }

            --indexImage2;

            --imagePixelsIndex2;

            --srcKernelIndexIterBefore;
            --srcKernelIndex;
            --srcKernelIndexIterNext;
            --srcKernelIndexIterNextNext;
        }

        srcKernelIndexIterBefore -= 2;
        srcKernelIndex -= 2;
        srcKernelIndexIterNext -= 2;
        srcKernelIndexIterNextNext = srcKernelIndexIterNext;
        srcKernelIndexIterNext = srcKernelIndex;
        srcKernelIndex = srcKernelIndexIterBefore;
        srcKernelIndexIterBefore -= dim;

        imagePixelsIndex->blue = *indexImage;
        --indexImage;
        imagePixelsIndex->green = *indexImage;
        --indexImage;
        imagePixelsIndex->red = *indexImage;

        imagePixelsIndex2->blue = *indexImage2;
        --indexImage2;
        imagePixelsIndex2->green = *indexImage2;
        --indexImage2;
        imagePixelsIndex2->red = *indexImage2;
        --indexImage2;

        --imagePixelsIndex2;

        indexImage = indexImage2;
        indexImage2 -= dimMult3;
        imagePixelsIndex = imagePixelsIndex2;
        imagePixelsIndex2 -= dim;
    }

    for(++i; i > 0; --i){
        imagePixelsIndex->blue = *indexImage;
        --indexImage;
        imagePixelsIndex->green = *indexImage;
        --indexImage;
        imagePixelsIndex->red = *indexImage;
        --indexImage;

        --imagePixelsIndex;
        
        colSum->red = (int) firstSrcIndex->red + (int) firstSrcIndexNext->red + (int) firstSrcIndexNextNext->red;
        colSum->green = (int) firstSrcIndex->green + (int) firstSrcIndexNext->green + (int) firstSrcIndexNextNext->green;
        colSum->blue = (int) firstSrcIndex->blue + (int) firstSrcIndexNext->blue + (int) firstSrcIndexNextNext->blue;

        colSumIterNext->red = (int) secondSrcIndex->red + (int) secondSrcIndexNext->red + (int) secondSrcIndexNextNext->red;
        colSumIterNext->green = (int) secondSrcIndex->green + (int) secondSrcIndexNext->green + (int) secondSrcIndexNextNext->green;
        colSumIterNext->blue = (int) secondSrcIndex->blue + (int) secondSrcIndexNext->blue + (int) secondSrcIndexNextNext->blue;

        firstSrcIndexNext = firstSrcIndex;
        firstSrcIndexNextNext = firstSrcIndexNext;
        firstSrcIndex -= dim;
        secondSrcIndexNext = secondSrcIndex;
        secondSrcIndexNextNext = secondSrcIndexNext;
        secondSrcIndex = firstSrcIndex - 1;
        
        pixel_sum *colSumIndex = startColToFillIndex;
        int j = startIAndJ;
        for (; j > 0; --j) {
            colSumIndex->red = (int)srcKernelIndex->red + (int)srcKernelIndexIterNext->red + (int)srcKernelIndexIterNextNext->red;
            colSumIndex->green = (int)srcKernelIndex->green + (int)srcKernelIndexIterNext->green + (int)srcKernelIndexIterNextNext->green;
            colSumIndex->blue = (int)srcKernelIndex->blue + (int)srcKernelIndexIterNext->blue + (int)srcKernelIndexIterNextNext->blue;

            ++colSumIndex;
            if (colSumIndex == endColToFillIndex) { colSumIndex = colSum; }

            pixel_sum sum = { colSum->red + colSumIterNext->red + colSumIterNextNext->red,
            colSum->green + colSumIterNext->green + colSumIterNextNext->green,
            colSum->blue + colSumIterNext->blue + colSumIterNextNext->blue };

            if (sum.blue < 0) {
                *indexImage = imagePixelsIndex->blue = 0;
            }
            else if (sum.blue > maxkernelScale) {
                *indexImage = imagePixelsIndex->blue = 255;
            }
            else {
                *indexImage = imagePixelsIndex->blue = (unsigned char)(sum.blue / kernelScale);
            }

            --indexImage;

            if (sum.green < 0) {
                *indexImage = imagePixelsIndex->green = 0;
            }
            else if (sum.green > maxkernelScale) {
                *indexImage = imagePixelsIndex->green = 255;
            }
            else {
                *indexImage = imagePixelsIndex->green = (unsigned char)(sum.green / kernelScale);
            }

            --indexImage;

            if (sum.red < 0) {
                *indexImage = imagePixelsIndex->red = 0;
            }
            else if (sum.red > maxkernelScale) {
                *indexImage = imagePixelsIndex->red = 255;
            }
            else {
                *indexImage = imagePixelsIndex->red = (unsigned char)(sum.red / kernelScale);
            }

            --indexImage;

            --imagePixelsIndex;

            --srcKernelIndex;
            --srcKernelIndexIterNext;
            --srcKernelIndexIterNextNext;

        }

        srcKernelIndex -= 2;
        srcKernelIndexIterNext -= 2;
        srcKernelIndexIterNextNext -= 2;

        imagePixelsIndex->blue = *indexImage;
        --indexImage;
        imagePixelsIndex->green = *indexImage;
        --indexImage;
        imagePixelsIndex->red = *indexImage;
        --indexImage;

        --imagePixelsIndex;
    }
}

void getSharppenNoFilter(int dim, pixel *src, Image *charsImg){
    int srcSize = dim*dim;
    int index = srcSize - 1;
    int dimMult2 = dim<<1;
    int dimMult3 = dimMult2 + dim;
    int sumSize = srcSize-(dimMult2<<1)+4;
    int startIAndJ = dim - 2;
    int startIAndJLoopUnrolling = startIAndJ - 1;
    unsigned char* indexImage = charsImg->data + (srcSize<<1) + srcSize -dimMult3 - 4;
    pixel * srcKernelIndex = src + sumSize + dimMult2 -7;
    pixel* srcKernelIndexIterNext = srcKernelIndex + dim;
    pixel* srcKernelIndexIterNextNext = srcKernelIndexIterNext + dim;
    pixel* firstSrcIndex = src + srcSize -1 - dimMult2;
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
    // row 2:
    pixel_sum colSum2[3];
    pixel_sum* colSumIterNext2 = colSum2 + 1;
    pixel_sum* colSumIterNextNext2 = colSumIterNext2 + 1;
    pixel_sum* startColToFillIndex2 = colSum2 + 2;
    pixel* firstSrcIndexBefore = firstSrcIndex - dim;
    pixel* secondSrcIndexBefore = secondSrcIndex - dim;
    pixel* srcKernelIndexIterBefore = srcKernelIndex - dim;
    unsigned char* indexImage2 = indexImage -dimMult3;

    int i =  startIAndJLoopUnrolling;
    for(; i > 0; i-=2){
        int rNext2 = (int) secondSrcIndex->red;
        int gNext2 = (int) secondSrcIndex->green;
        int bNext2 = (int) secondSrcIndex->blue;

        int rNext = (int) secondSrcIndexNext->red;
        int gNext = (int) secondSrcIndexNext->green;
        int bNext = (int) secondSrcIndexNext->blue;

        int helpSumR = (int) firstSrcIndex->red + (int) firstSrcIndexNext->red;
        int helpSumG = (int) firstSrcIndex->green + (int) firstSrcIndexNext->green;
        int helpSumB = (int) firstSrcIndex->blue + (int) firstSrcIndexNext->blue;
    
        colSum->red = helpSumR + (int) firstSrcIndexNextNext->red;
        colSum->green = helpSumG + (int) firstSrcIndexNextNext->green;
        colSum->blue = helpSumB + (int) firstSrcIndexNextNext->blue;

        colSum2->red = (int) firstSrcIndexBefore->red + helpSumR;
        colSum2->green = (int) firstSrcIndexBefore->green + helpSumG;
        colSum2->blue = (int) firstSrcIndexBefore->blue + helpSumB;

        helpSumR = rNext2 + rNext;
        helpSumG = gNext2 + gNext;
        helpSumB = bNext2 + bNext;

        colSumIterNext->red = helpSumR + (int) secondSrcIndexNextNext->red;
        colSumIterNext->green = helpSumG + (int) secondSrcIndexNextNext->green;
        colSumIterNext->blue = helpSumB + (int) secondSrcIndexNextNext->blue;

        colSumIterNext2->red = (int) secondSrcIndexBefore->red + helpSumR;
        colSumIterNext2->green = (int) secondSrcIndexBefore->green + helpSumG;
        colSumIterNext2->blue = (int) secondSrcIndexBefore->blue + helpSumB;

        firstSrcIndexNext = firstSrcIndexBefore;
        firstSrcIndexNextNext = firstSrcIndex;
        firstSrcIndexBefore -= dim;
        firstSrcIndex = firstSrcIndexBefore;
        firstSrcIndexBefore -= dim;

        secondSrcIndexNext = secondSrcIndexBefore;
        secondSrcIndexNextNext = secondSrcIndex;
        secondSrcIndex = firstSrcIndex - 1;
        secondSrcIndexBefore = firstSrcIndexBefore - 1;
        
        pixel_sum *colSumIndex = startColToFillIndex;
        pixel_sum *colSumIndex2 = startColToFillIndex2;
        int j = startIAndJ;
        for(; j>0 ; --j){
            int r = rNext;
            int g = gNext;
            int b = bNext;

            int r2 = rNext2;
            int g2 = gNext2;
            int b2 = bNext2;

            rNext = (int) srcKernelIndexIterNext->red;
            gNext = (int) srcKernelIndexIterNext->green;
            bNext = (int) srcKernelIndexIterNext->blue;

            rNext2 = (int) srcKernelIndex->red;
            gNext2 = (int) srcKernelIndex->green;
            bNext2 = (int) srcKernelIndex->blue;

            helpSumR = rNext2 + rNext;
            helpSumG = gNext2 + gNext;
            helpSumB = bNext2 + bNext;

            colSumIndex->red = helpSumR + (int) srcKernelIndexIterNextNext->red;
            colSumIndex->green = helpSumG + (int) srcKernelIndexIterNextNext->green;
            colSumIndex->blue = helpSumB + (int) srcKernelIndexIterNextNext->blue;

            colSumIndex2->red = (int) srcKernelIndexIterBefore->red + helpSumR;
            colSumIndex2->green = (int) srcKernelIndexIterBefore->green + helpSumG;
            colSumIndex2->blue = (int) srcKernelIndexIterBefore->blue + helpSumB;

            ++colSumIndex;
            ++colSumIndex2;
            if(colSumIndex == endColToFillIndex) {
                colSumIndex = colSum;
                colSumIndex2 = colSum2;
            }

            pixel_sum sum = {colSum->red + colSumIterNext->red + colSumIterNextNext-> red,
            colSum->green + colSumIterNext->green + colSumIterNextNext-> green,
            colSum->blue + colSumIterNext->blue + colSumIterNextNext-> blue};

            pixel_sum sum2 = {colSum2->red + colSumIterNext2->red + colSumIterNextNext2-> red,
            colSum2->green + colSumIterNext2->green + colSumIterNextNext2-> green,
            colSum2->blue + colSumIterNext2->blue + colSumIterNextNext2-> blue};

            sum.red = (r<<3) + (r<<1) - sum.red;
            sum.green = (g<<3) + (g<<1) - sum.green;
            sum.blue = (b<<3) + (b<<1) - sum.blue;

            sum2.red = (r2<<3) + (r2<<1) - sum2.red;
            sum2.green = (g2<<3) + (g2<<1) - sum2.green;
            sum2.blue = (b2<<3) + (b2<<1) - sum2.blue;

            if (sum.blue < 0){sum.blue = 0;} else if(sum.blue > 255) {sum.blue = 255;}
            *indexImage = (unsigned char) sum.blue;
            --indexImage;
            if (sum.green < 0){sum.green = 0;} else if(sum.green > 255) {sum.green = 255;}
            *indexImage = (unsigned char) sum.green;
            --indexImage;
            if (sum.red < 0){sum.red = 0;} else if(sum.red > 255) {sum.red = 255;}
            *indexImage = (unsigned char) sum.red;
            --indexImage;

            if (sum2.blue < 0){sum2.blue = 0;} else if(sum2.blue > 255) {sum2.blue = 255;}
            *indexImage2 = (unsigned char) sum2.blue;
            --indexImage2;
            if (sum2.green < 0){sum2.green = 0;} else if(sum2.green > 255) {sum2.green = 255;}
            *indexImage2 = (unsigned char) sum2.green;
            --indexImage2;
            if (sum2.red < 0){sum2.red = 0;} else if(sum2.red > 255) {sum2.red = 255;}
            *indexImage2 = (unsigned char) sum2.red;
            --indexImage2;

            --srcKernelIndexIterBefore;
            --srcKernelIndex;
            --srcKernelIndexIterNext;
            --srcKernelIndexIterNextNext;
        }

        indexImage2 -= 6;
        indexImage = indexImage2;
        indexImage2 -= dimMult3;
        srcKernelIndexIterBefore -= 2;
        srcKernelIndex -= 2;
        srcKernelIndexIterNext -= 2;
        srcKernelIndexIterNextNext = srcKernelIndexIterNext;
        srcKernelIndexIterNext = srcKernelIndex;
        srcKernelIndex = srcKernelIndexIterBefore;
        srcKernelIndexIterBefore -= dim;
    }
    for(++i; i > 0; --i){
        int rNext = (int) secondSrcIndexNext->red;
        int gNext = (int) secondSrcIndexNext->green;
        int bNext = (int) secondSrcIndexNext->blue;
    
        colSum->red = (int) firstSrcIndex->red + (int) firstSrcIndexNext->red + (int) firstSrcIndexNextNext->red;
        colSum->green = (int) firstSrcIndex->green + (int) firstSrcIndexNext->green + (int) firstSrcIndexNextNext->green;
        colSum->blue = (int) firstSrcIndex->blue + (int) firstSrcIndexNext->blue + (int) firstSrcIndexNextNext->blue;

        colSumIterNext->red = (int) secondSrcIndex->red + rNext + (int) secondSrcIndexNextNext->red;
        colSumIterNext->green = (int) secondSrcIndex->green + gNext + (int) secondSrcIndexNextNext->green;
        colSumIterNext->blue = (int) secondSrcIndex->blue + bNext + (int) secondSrcIndexNextNext->blue;

        firstSrcIndexNext = firstSrcIndex;
        firstSrcIndexNextNext = firstSrcIndexNext;
        firstSrcIndex -= dim;
        secondSrcIndexNext = secondSrcIndex;
        secondSrcIndexNextNext = secondSrcIndexNext;
        secondSrcIndex = firstSrcIndex - 1;
        
        pixel_sum *colSumIndex = startColToFillIndex;
        int j = startIAndJ;
        for(; j>0 ; --j){
            int r = rNext;
            int g = gNext;
            int b = bNext;

            rNext = (int) srcKernelIndexIterNext->red;
            gNext = (int) srcKernelIndexIterNext->green;
            bNext = (int) srcKernelIndexIterNext->blue;

            colSumIndex->red = (int) srcKernelIndex->red + rNext + (int) srcKernelIndexIterNextNext->red;
            colSumIndex->green = (int) srcKernelIndex->green + gNext + (int) srcKernelIndexIterNextNext->green;
            colSumIndex->blue = (int) srcKernelIndex->blue + bNext + (int) srcKernelIndexIterNextNext->blue;

            ++colSumIndex;
            if(colSumIndex == endColToFillIndex) {colSumIndex = colSum;}

            pixel_sum sum = {colSum->red + colSumIterNext->red + colSumIterNextNext-> red,
            colSum->green + colSumIterNext->green + colSumIterNextNext-> green,
            colSum->blue + colSumIterNext->blue + colSumIterNextNext-> blue};

            sum.red = (r<<3) + (r<<1) - sum.red;
            sum.green = (g<<3) + (g<<1) - sum.green;
            sum.blue = (b<<3) + (b<<1) - sum.blue;

            if (sum.blue < 0){sum.blue = 0;} else if(sum.blue > 255) {sum.blue = 255;}
            *indexImage = (unsigned char) sum.blue;
            --indexImage;
            if (sum.green < 0){sum.green = 0;} else if(sum.green > 255) {sum.green = 255;}
            *indexImage = (unsigned char) sum.green;
            --indexImage;
            if (sum.red < 0){sum.red = 0;} else if(sum.red > 255) {sum.red = 255;}
            *indexImage = (unsigned char) sum.red;
            --indexImage;

            --srcKernelIndex;
            --srcKernelIndexIterNext;
            --srcKernelIndexIterNextNext;
        }

        indexImage -= 6;
        srcKernelIndex -= 2;
        srcKernelIndexIterNext -= 2;
        srcKernelIndexIterNextNext -= 2;
    }
}

void getBlurWithFilter(int dim, pixel *src, pixel *imagePixels, Image *charsImg, int kernelScale){
    int srcSize = dim*dim;
    int index = srcSize - 1;
    int dimMult2 = dim<<1;
    int sumSize = srcSize-(dimMult2<<1)+4;
    int maxkernelScale = 255 * kernelScale;
    int startIAndJ = dim - 2;
    unsigned char* indexImage = charsImg->data + (srcSize<<1) + srcSize -dimMult2 -dim - 1;
    pixel* srcKernelIndex = src + sumSize + dimMult2 -7;
    pixel* srcKernelIndexIterNext = srcKernelIndex + dim;
    pixel* srcKernelIndexIterNextNext = srcKernelIndexIterNext + dim;
    pixel* imagePixelsIndex = imagePixels + index - dim;
    pixel* firstSrcIndex = src + srcSize -1 - dimMult2;
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
    int minIntensities[3];
    int minRaw[3];
    pixel_sum minPixels[3];
    int maxIntensities[3];
    int maxRaw[3];
    pixel_sum maxPixels[3];

    int i = startIAndJ;
    for(; i > 0; --i){
        imagePixelsIndex->blue = *indexImage;
        --indexImage;
        imagePixelsIndex->green = *indexImage;
        --indexImage;
        imagePixelsIndex->red = *indexImage;
        --indexImage;

        --imagePixelsIndex;

        int* minIntensitiesIndex = minIntensities;
        int* minRawIndex = minRaw;
        pixel_sum* minPixelsIndex = minPixels;
        int* maxIntensitiesIndex = maxIntensities;
        int* maxRawIndex = maxRaw;
        pixel_sum* maxPixelsIndex = maxPixels;

        {
            pixel_sum pixel1 = {(int) firstSrcIndex->red, (int) firstSrcIndex->green, (int) firstSrcIndex->blue};
            pixel_sum pixel2 = {(int) firstSrcIndexNext->red, (int) firstSrcIndexNext->green, (int) firstSrcIndexNext->blue};
            pixel_sum pixel3 = {(int) firstSrcIndexNextNext->red, (int) firstSrcIndexNextNext->green, (int) firstSrcIndexNextNext->blue};
            
            colSum->red = pixel1.red + pixel2.red + pixel3.red;
            colSum->green = pixel1.green + pixel2.green + pixel3.green;
            colSum->blue = pixel1.blue + pixel2.blue + pixel3.blue;

            int intensity1 = pixel1.red + pixel1.green + pixel1.blue;
            int intensity2 = pixel2.red + pixel2.green + pixel2.blue;
            int intensity3 = pixel3.red + pixel3.green + pixel3.blue;

            if(intensity1 == intensity2 && intensity2 == intensity3){
                *maxIntensitiesIndex = intensity1;
                *maxRawIndex = 0;
                *maxPixelsIndex = pixel1;

                *minIntensitiesIndex = intensity3;
                *minRawIndex = 2;
                *minPixelsIndex = pixel3;
            }else if(intensity2 <= intensity1){
                if(intensity3 <= intensity2){
                    *maxIntensitiesIndex = intensity1;
                    *maxRawIndex = 0;
                    *maxPixelsIndex = pixel1;

                    *minIntensitiesIndex = intensity3;
                    *minRawIndex = 2;
                    *minPixelsIndex = pixel3;
                } else {
                    if(intensity3 <= intensity1){
                        *maxIntensitiesIndex = intensity1;
                        *maxRawIndex = 0;
                        *maxPixelsIndex = pixel1;
                    } else {
                        *maxIntensitiesIndex = intensity3;
                        *maxRawIndex = 2;
                        *maxPixelsIndex = pixel3;
                    }

                    *minIntensitiesIndex = intensity2;
                    *minRawIndex = 1;
                    *minPixelsIndex = pixel2;
                }
            }else if(intensity3 > intensity1){
                if (intensity3 > intensity2) {

                    *maxIntensitiesIndex = intensity3;
                    *maxRawIndex = 2;
                    *maxPixelsIndex = pixel3;
                } else {
                    *maxIntensitiesIndex = intensity2;
                    *maxRawIndex = 1;
                    *maxPixelsIndex = pixel2;
                }

                *minIntensitiesIndex = intensity1;
                *minRawIndex = 0;
                *minPixelsIndex = pixel1;
            }else{

                *maxIntensitiesIndex = intensity2;
                *maxRawIndex = 1;
                *maxPixelsIndex = pixel2;
                
                *minIntensitiesIndex = intensity3;
                *minRawIndex = 2;
                *minPixelsIndex = pixel3;
            }
        }
        ++minIntensitiesIndex;
        ++minRawIndex;
        ++minPixelsIndex;
        ++maxIntensitiesIndex;
        ++maxRawIndex;
        ++maxPixelsIndex;
        {
            pixel_sum pixel1 = {(int) secondSrcIndex->red, (int) secondSrcIndex->green, (int) secondSrcIndex->blue};
            pixel_sum pixel2 = {(int) secondSrcIndexNext->red, (int) secondSrcIndexNext->green, (int) secondSrcIndexNext->blue};
            pixel_sum pixel3 = {(int) secondSrcIndexNextNext->red, (int) secondSrcIndexNextNext->green, (int) secondSrcIndexNextNext->blue};
            
            colSumIterNext->red = pixel1.red + pixel2.red + pixel3.red;
            colSumIterNext->green = pixel1.green + pixel2.green + pixel3.green;
            colSumIterNext->blue = pixel1.blue + pixel2.blue + pixel3.blue;

            int intensity1 = pixel1.red + pixel1.green + pixel1.blue;
            int intensity2 = pixel2.red + pixel2.green + pixel2.blue;
            int intensity3 = pixel3.red + pixel3.green + pixel3.blue;

            if(intensity1 == intensity2 && intensity2 == intensity3){
                *maxIntensitiesIndex = intensity1;
                *maxRawIndex = 0;
                *maxPixelsIndex = pixel1;

                *minIntensitiesIndex = intensity3;
                *minRawIndex = 2;
                *minPixelsIndex = pixel3;
            }else if(intensity2 <= intensity1){
                if(intensity3 <= intensity2){
                    *maxIntensitiesIndex = intensity1;
                    *maxRawIndex = 0;
                    *maxPixelsIndex = pixel1;

                    *minIntensitiesIndex = intensity3;
                    *minRawIndex = 2;
                    *minPixelsIndex = pixel3;
                } else {
                    if(intensity3 <= intensity1){
                        *maxIntensitiesIndex = intensity1;
                        *maxRawIndex = 0;
                        *maxPixelsIndex = pixel1;
                    } else {
                        *maxIntensitiesIndex = intensity3;
                        *maxRawIndex = 2;
                        *maxPixelsIndex = pixel3;
                    }

                    *minIntensitiesIndex = intensity2;
                    *minRawIndex = 1;
                    *minPixelsIndex = pixel2;
                }
            }else if(intensity3 > intensity1){
                if (intensity3 > intensity2) {

                    *maxIntensitiesIndex = intensity3;
                    *maxRawIndex = 2;
                    *maxPixelsIndex = pixel3;
                } else {
                    *maxIntensitiesIndex = intensity2;
                    *maxRawIndex = 1;
                    *maxPixelsIndex = pixel2;
                }

                *minIntensitiesIndex = intensity1;
                *minRawIndex = 0;
                *minPixelsIndex = pixel1;
            }else{

                *maxIntensitiesIndex = intensity2;
                *maxRawIndex = 1;
                *maxPixelsIndex = pixel2;
                
                *minIntensitiesIndex = intensity3;
                *minRawIndex = 2;
                *minPixelsIndex = pixel3;
            }
        }
        ++minIntensitiesIndex;
        ++minRawIndex;
        ++minPixelsIndex;
        ++maxIntensitiesIndex;
        ++maxRawIndex;
        ++maxPixelsIndex;
        
        firstSrcIndexNextNext = firstSrcIndexNext;
        firstSrcIndexNext = firstSrcIndex; 
        firstSrcIndex -= dim;
        secondSrcIndexNextNext = secondSrcIndexNext;
        secondSrcIndexNext = secondSrcIndex;
        secondSrcIndex = firstSrcIndex - 1;
        
        pixel_sum *colSumIndex = startColToFillIndex;
        int j = startIAndJ;
        for(; j>0 ; --j){
            pixel_sum pixel1 = {(int) srcKernelIndex->red, (int) srcKernelIndex->green, (int) srcKernelIndex->blue};
            pixel_sum pixel2 = {(int) srcKernelIndexIterNext->red, (int) srcKernelIndexIterNext->green, (int) srcKernelIndexIterNext->blue};
            pixel_sum pixel3 = {(int) srcKernelIndexIterNextNext->red, (int) srcKernelIndexIterNextNext->green, (int) srcKernelIndexIterNextNext->blue};

            colSumIndex->red = pixel1.red + pixel2.red + pixel3.red;
            colSumIndex->green = pixel1.green + pixel2.green + pixel3.green;
            colSumIndex->blue = pixel1.blue + pixel2.blue + pixel3.blue;

            int intensity1 = pixel1.red + pixel1.green + pixel1.blue;
            int intensity2 = pixel2.red + pixel2.green + pixel2.blue;
            int intensity3 = pixel3.red + pixel3.green + pixel3.blue;

            if(intensity1 == intensity2 && intensity2 == intensity3){
                *maxIntensitiesIndex = intensity1;
                *maxRawIndex = 0;
                *maxPixelsIndex = pixel1;

                *minIntensitiesIndex = intensity3;
                *minRawIndex = 2;
                *minPixelsIndex = pixel3;
            }else if(intensity2 <= intensity1){
                if(intensity3 <= intensity2){
                    *maxIntensitiesIndex = intensity1;
                    *maxRawIndex = 0;
                    *maxPixelsIndex = pixel1;

                    *minIntensitiesIndex = intensity3;
                    *minRawIndex = 2;
                    *minPixelsIndex = pixel3;
                } else {
                    if(intensity3 <= intensity1){
                        *maxIntensitiesIndex = intensity1;
                        *maxRawIndex = 0;
                        *maxPixelsIndex = pixel1;
                    } else {
                        *maxIntensitiesIndex = intensity3;
                        *maxRawIndex = 2;
                        *maxPixelsIndex = pixel3;
                    }

                    *minIntensitiesIndex = intensity2;
                    *minRawIndex = 1;
                    *minPixelsIndex = pixel2;
                }
            }else if(intensity3 > intensity1){
                if (intensity3 > intensity2) {

                    *maxIntensitiesIndex = intensity3;
                    *maxRawIndex = 2;
                    *maxPixelsIndex = pixel3;
                } else {
                    *maxIntensitiesIndex = intensity2;
                    *maxRawIndex = 1;
                    *maxPixelsIndex = pixel2;
                }

                *minIntensitiesIndex = intensity1;
                *minRawIndex = 0;
                *minPixelsIndex = pixel1;
            }else{

                *maxIntensitiesIndex = intensity2;
                *maxRawIndex = 1;
                *maxPixelsIndex = pixel2;
                
                *minIntensitiesIndex = intensity3;
                *minRawIndex = 2;
                *minPixelsIndex = pixel3;
            }

            int colToFillIndex = minRawIndex - minRaw;
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
            pixel_sum min = minPixels[minIndex];
            pixel_sum max = maxPixels[maxIndex];

            pixel_sum sum = {colSum->red + colSumIterNext->red + colSumIterNextNext->red - min.red- max.red,
            colSum->green + colSumIterNext->green + colSumIterNextNext->green - min.green- max.green,
            colSum->blue + colSumIterNext->blue + colSumIterNextNext->blue - min.blue- max.blue};

            if (sum.blue < 0){
                *indexImage = imagePixelsIndex->blue = 0;
            } else if(sum.blue > maxkernelScale) {
                *indexImage = imagePixelsIndex->blue = 255;
            } else {
                *indexImage = imagePixelsIndex->blue = (unsigned char) (sum.blue / kernelScale);
            }

            --indexImage;

            if (sum.green < 0){
                *indexImage = imagePixelsIndex->green = 0;
            } else if(sum.green > maxkernelScale) {
                *indexImage = imagePixelsIndex->green = 255;
            } else {
                *indexImage = imagePixelsIndex->green = (unsigned char) (sum.green / kernelScale);
            }

            --indexImage;

            if (sum.red < 0){
                *indexImage = imagePixelsIndex->red = 0;
            } else if(sum.red > maxkernelScale) {
                *indexImage = imagePixelsIndex->red = 255;
            } else {
                *indexImage = imagePixelsIndex->red = (unsigned char) (sum.red / kernelScale);
            }

            --indexImage;

            --imagePixelsIndex;

            ++minIntensitiesIndex;
            ++minRawIndex;
            ++minPixelsIndex;
            ++maxIntensitiesIndex;
            ++maxRawIndex;
            ++maxPixelsIndex;

            ++colSumIndex;
            if(colSumIndex == endColToFillIndex) {
                colSumIndex = colSum;
                minIntensitiesIndex = minIntensities;
                minRawIndex = minRaw;
                minPixelsIndex = minPixels;
                maxIntensitiesIndex = maxIntensities;
                maxRawIndex = maxRaw;
                maxPixelsIndex = maxPixels;
            }

            --srcKernelIndex;
            --srcKernelIndexIterNext;
            --srcKernelIndexIterNextNext;
        }
        srcKernelIndex -= 2;
        srcKernelIndexIterNext -= 2;
        srcKernelIndexIterNextNext -= 2;

        imagePixelsIndex->blue = *indexImage;
        --indexImage;
        imagePixelsIndex->green = *indexImage;
        --indexImage;
        imagePixelsIndex->red = *indexImage;
        --indexImage;

        --imagePixelsIndex;
    }
}
void charsToPixels(Image *charsImg, pixel* pixels, pixel* frame) {
    //calculating index efficiently
    unsigned char* indexImage = charsImg->data;
    pixel* indexPixels = pixels;
    pixel* indexFrame = frame;
    int row, col;
    int startColRow = m - 2;

    for (col = m; col > 0; --col) {
        indexPixels->red = indexFrame->red = *indexImage;
        ++indexImage;
        indexPixels->green = indexFrame->green = *indexImage;
        ++indexImage;
        indexPixels->blue = indexFrame->blue = *indexImage;
        ++indexImage;

        ++indexPixels;
        ++indexFrame;
    }

    for (row = startColRow; row > 0; --row) {
        for (col = m; col > 0; --col) {

            indexPixels->red = *indexImage;
            ++indexImage;
            indexPixels->green = *indexImage;
            ++indexImage;
            indexPixels->blue = *indexImage;
            ++indexImage;

            ++indexPixels;
        }
    }

    indexFrame = frame + m * (m - 1);
    for (col = m; col > 0; --col) {
        indexPixels->red = indexFrame->red = *indexImage;
        ++indexImage;
        indexPixels->green = indexFrame->green = *indexImage;
        ++indexImage;
        indexPixels->blue = indexFrame->blue = *indexImage;
        ++indexImage;

        ++indexPixels;
        ++indexFrame;
    }
}

void myfunction(Image *image, char* srcImgpName, char* blurRsltImgName, char* sharpRsltImgName, char* filteredBlurRsltImgName, char* filteredSharpRsltImgName, char flag) {
    if (m < 3) {
        if (flag == '1') {

        // write result image to file
        writeBMP(image, srcImgpName, blurRsltImgName);
        
        // write result image to file
        writeBMP(image, srcImgpName, sharpRsltImgName);
        return;
        }

        // write result image to file
        writeBMP(image, srcImgpName, filteredBlurRsltImgName);

        // write result image to file
        writeBMP(image, srcImgpName, filteredSharpRsltImgName);
        return;
    }

    size_t size = m*n*sizeof(pixel);
    pixel* pixelsImg = malloc(size);
    pixel* backupOrg = malloc(size);
    //caculating length once:
    charsToPixels(image, backupOrg, pixelsImg);
    if (flag == '1') {
        
        getBlurNoFilter(m, backupOrg, pixelsImg, image, 9);

        // write result image to file
        writeBMP(image, srcImgpName, blurRsltImgName);

        // sharpen the resulting image
        getSharppenNoFilter(m, pixelsImg, image);
        
        // write result image to file
        writeBMP(image, srcImgpName, sharpRsltImgName); 

        return;
    }

    getBlurWithFilter(m, backupOrg, pixelsImg, image, 7);

    // write result image to file
    writeBMP(image, srcImgpName, filteredBlurRsltImgName);

    // sharpen the resulting image
    getSharppenNoFilter(m, pixelsImg, image);

    // write result image to file
    writeBMP(image, srcImgpName, filteredSharpRsltImgName); 

    free(pixelsImg);
    free(backupOrg);
}
