rm Blur.bmp
rm Sharpen.bmp
rm Filtered_Blur.bmp
rm Filtered_Sharpen.bmp
rm report0.txt
rm report1.txt
rm report2.txt
make
# echo "*************************./showBMP 1x1.bmp *************************"
# ./showBMP 1x1.bmp
# echo "*************************./showBMP 1x1.bmp 1 *************************"
# ./showBMP 1x1.bmp 1
# echo "*************************./showBMP 1x1.bmp 2 *************************"
# ./showBMP 1x1.bmp 2
# echo "*************************cmp 1x1.bmp*************************"
# cmp Blur.bmp 1x1/Blur.bmp
# cmp Sharpen.bmp 1x1/Sharpen.bmp
# cmp Filtered_Blur.bmp 1x1/Filtered_Blur.bmp
# cmp Filtered_Sharpen.bmp 1x1/Filtered_Sharpen.bmp
# echo "*************************end cmp 1x1.bmp*************************"
# echo "*************************./showBMP 2x2.bmp *************************"
# ./showBMP 2x2.bmp
# echo "*************************./showBMP 2x2.bmp 1 *************************"
# ./showBMP 2x2.bmp 1
# echo "*************************./showBMP 2x2.bmp 2 *************************"
# ./showBMP 2x2.bmp 2
# echo "*************************cmp 2x2.bmp*************************"
# cmp Blur.bmp 2x2/Blur.bmp
# cmp Sharpen.bmp 2x2/Sharpen.bmp
# cmp Filtered_Blur.bmp 2x2/Filtered_Blur.bmp
# cmp Filtered_Sharpen.bmp 2x2/Filtered_Sharpen.bmp
# echo "*************************end cmp 2x2.bmp*************************"
echo "*************************./showBMP gibson_500.bmp *************************"
./showBMP gibson_500.bmp
gprof showBMP gmon.out > report0.txt
echo "*************************./showBMP gibson_500.bmp 1 *************************"
./showBMP gibson_500.bmp 1
gprof showBMP gmon.out > report1.txt
echo "*************************./showBMP gibson_500.bmp 2 *************************"
./showBMP gibson_500.bmp 2
gprof showBMP gmon.out > report2.txt
echo "*************************cmp gibson_500.bmp*************************"
cmp Blur.bmp expected/Blur_correct.bmp
cmp Sharpen.bmp expected/Sharpen_correct.bmp
cmp Filtered_Blur.bmp expected/Filtered_Blur_correct.bmp
cmp Filtered_Sharpen.bmp expected/Filtered_Sharpen_correct.bmp
echo "*************************end cmp gibson_500.bmp*************************"
rm gmon.out
make clean