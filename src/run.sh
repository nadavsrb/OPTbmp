rm Blur.bmp
rm Sharpen.bmp
rm Filtered_Blur.bmp
rm Filtered_Sharpen.bmp
rm report0.txt
rm report1.txt
rm report2.txt
make
echo "*************************./showBMP gibson_500.bmp *************************"
./showBMP gibson_500.bmp
gprof showBMP gmon.out > report0.txt
echo "*************************./showBMP gibson_500.bmp 1 *************************"
./showBMP gibson_500.bmp 1
gprof showBMP gmon.out > report1.txt
echo "*************************./showBMP gibson_500.bmp 2 *************************"
./showBMP gibson_500.bmp 2
gprof showBMP gmon.out > report2.txt
echo "*************************cmp files*************************"
cmp Blur.bmp expected/Blur_correct.bmp
cmp Sharpen.bmp expected/Sharpen_correct.bmp
cmp Filtered_Blur.bmp expected/Filtered_Blur_correct.bmp
cmp Filtered_Sharpen.bmp expected/Filtered_Sharpen_correct.bmp
echo "*************************end cmp files*************************"
rm gmon.out
make clean