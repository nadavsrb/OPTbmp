clear
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
# # echo "*************************./showBMP wild.bmp *************************"
# # ./showBMP wild.bmp
# # echo "*************************./showBMP wild.bmp 1 *************************"
# # ./showBMP wild.bmp 1
# # echo "*************************./showBMP wild.bmp 2 *************************"
# # ./showBMP wild.bmp 2

# # # ./showBMP wild.bmp  wild/Filtered_Blur.bmp Filtered_Blur.bmp
# # echo "*************************cmp wild.bmp*************************"
# # cmp Blur.bmp wild/Blur.bmp
# # cmp Sharpen.bmp wild/Sharpen.bmp
# # cmp Filtered_Blur.bmp wild/Filtered_Blur.bmp
# # cmp Filtered_Sharpen.bmp wild/Filtered_Sharpen.bmp
# # echo "*************************end cmp wild.bmp*************************"
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
echo "*************************./showBMP Blur1.bmp *************************"
./showBMP Blur1.bmp
gprof showBMP gmon.out > report0.txt
echo "*************************./showBMP Blur1.bmp 1 *************************"
./showBMP Blur1.bmp 1
gprof showBMP gmon.out > report1.txt
echo "*************************./showBMP Blur1.bmp 2 *************************"
./showBMP Blur1.bmp 2
gprof showBMP gmon.out > report2.txt
echo "*************************cmp Blur1.bmp*************************"
cmp Blur.bmp b/Blur.bmp
cmp Sharpen.bmp b/Sharpen.bmp
cmp Filtered_Blur.bmp b/Filtered_Blur.bmp
cmp Filtered_Sharpen.bmp b/Filtered_Sharpen.bmp
echo "*************************end cmp Blur1.bmp*************************"
echo "*************************./showBMP Sharpen1.bmp *************************"
./showBMP Sharpen1.bmp
gprof showBMP gmon.out > report0.txt
echo "*************************./showBMP Sharpen1.bmp 1 *************************"
./showBMP Sharpen1.bmp 1
gprof showBMP gmon.out > report1.txt
echo "*************************./showBMP Sharpen1.bmp 2 *************************"
./showBMP Sharpen1.bmp 2
gprof showBMP gmon.out > report2.txt
echo "*************************cmp Sharpen1.bmp*************************"
cmp Blur.bmp s/Blur.bmp
cmp Sharpen.bmp s/Sharpen.bmp
cmp Filtered_Blur.bmp s/Filtered_Blur.bmp
cmp Filtered_Sharpen.bmp s/Filtered_Sharpen.bmp
echo "*************************end cmp Sharpen1.bmp*************************"
echo "*************************./showBMP Filtered_Blur1.bmp *************************"
./showBMP Filtered_Blur1.bmp
gprof showBMP gmon.out > report0.txt
echo "*************************./showBMP Filtered_Blur1.bmp 1 *************************"
./showBMP Filtered_Blur1.bmp 1
gprof showBMP gmon.out > report1.txt
echo "*************************./showBMP Filtered_Blur1.bmp 2 *************************"
./showBMP Filtered_Blur1.bmp 2
gprof showBMP gmon.out > report2.txt
echo "*************************cmp Filtered_Blur1.bmp*************************"
cmp Blur.bmp bf/Blur.bmp
cmp Sharpen.bmp bf/Sharpen.bmp
cmp Filtered_Blur.bmp bf/Filtered_Blur.bmp
cmp Filtered_Sharpen.bmp bf/Filtered_Sharpen.bmp
echo "*************************end cmp Filtered_Blur1.bmp*************************"
echo "*************************./showBMP Filtered_Sharpen1.bmp *************************"
./showBMP Filtered_BluFiltered_Sharpen1r1.bmp
gprof showBMP gmon.out > report0.txt
echo "*************************./showBMP Filtered_Sharpen1.bmp 1 *************************"
./showBMP Filtered_Sharpen1.bmp 1
gprof showBMP gmon.out > report1.txt
echo "*************************./showBMP Filtered_Sharpen1.bmp 2 *************************"
./showBMP Filtered_Sharpen1.bmp 2
gprof showBMP gmon.out > report2.txt
echo "*************************cmp Filtered_Sharpen1.bmp*************************"
cmp Blur.bmp sf/Blur.bmp
cmp Sharpen.bmp sf/Sharpen.bmp
cmp Filtered_Blur.bmp sf/Filtered_Blur.bmp
cmp Filtered_Sharpen.bmp sf/Filtered_Sharpen.bmp
echo "*************************end cmp Filtered_Sharpen1.bmp*************************"
rm gmon.out
make clean