# radiator_deploy_sim

pdfcrop latex_equations.pdf latex_equations_cropped.pdf

convert -density 300   -background transparent   latex_equations_cropped.pdf   -trim +repage   -bordercolor transparent -border 5x5   latex_equations.png


pdfcrop state_system_matrix.pdf state_system_matrix_cropped.pdf

convert -density 300   -background transparent   state_system_matrix_cropped.pdf   -trim +repage   -bordercolor transparent -border 5x5   state_system_matrix.png



pdfcrop hard_stops.pdf hard_stops_cropped.pdf

convert -density 300   -background transparent   hard_stops_cropped.pdf   -trim +repage   -bordercolor transparent -borde
r 5x5   hard_stops.png