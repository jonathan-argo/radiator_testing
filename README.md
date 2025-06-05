# radiator_deploy_sim

pdfcrop latex_equations.pdf latex_equations_cropped.pdf

convert -density 300   -background transparent   latex_equations_cropped.pdf   -trim +repage   -bordercolor transparent -border 5x5   latex_equations.png


pdfcrop state_system_matrix.pdf state_system_matrix_cropped.pdf

convert -density 300   -background transparent   state_system_matrix_cropped.pdf   -trim +repage   -bordercolor transparent -border 5x5   state_system_matrix.png



pdfcrop hard_stops.pdf hard_stops_cropped.pdf

convert -density 300   -background transparent   hard_stops_cropped.pdf   -trim +repage   -bordercolor transparent -borde
r 5x5   hard_stops.png




git clone https://github.com/opencv/opencv.git
git -C opencv checkout 4.x

git clone https://github.com/opencv/opencv_contrib.git
git -C opencv checkout 4.x

cmake .. \
  -D CMAKE_BUILD_TYPE=Release \
  -D OPENCV_EXTRA_MODULES_PATH=../../opencv_contrib/modules \
  -D WITH_FFMPEG=ON \
  -D WITH_GTK=ON \
  -D BUILD_opencv_python=OFF \
  -D BUILD_opencv_python3=OFF \
  -D CMAKE_EXPORT_COMPILE_COMMANDS=ON

make -j$(nproc)