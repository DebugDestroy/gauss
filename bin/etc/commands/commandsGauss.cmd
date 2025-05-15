init

g 150 180 5 6 10
g 60 50 9 12 18
g 50 130 7 6 16
g 10 120 10 20 -10
g 80 75 10 15 -8
g 200 250 7 6 16
g 100 100 5 5 16
g 120 130 17 7 12
g 50 200 10 10 20
g 250 20 10 10 -29
g 170 50 20 10 20
g 10 10 10 10 7
g 290 290 10 10 7
g 290 10 10 10 7
g 10 290 10 30 7
g 230 90 20 10 7
g 180 130 10 20 10
g 200 190 20 10 7
g 230 90 10 10 7
g 230 90 10 10 7
g 30 200 10 10 7
g 100 250 8 8 10
g 200 100 10 10 -10
g 130 220 10 10 10
g 170 170 8 8 12
g 250 250 8 8 -20
g 270 150 6 6 25
g 220 30 10 10 10
g 140 40 6 10 -15
g 90 180 8 8 10
g 160 220 10 10 -12
g 70 250 6 6 10
g 270 70 10 10 -15
g 130 70 8 8 10
g 200 60 6 6 10
g 60 160 8 8 -10
g 110 30 6 6 10
g 180 30 8 8 8
g 250 120 6 6 10
g 40 90 6 6 -10
g 220 220 10 10 15


generate
gnuplot results/visualizations/gnuplot.png
bmp_write results/visualizations/Pole.bmp Full
bin 132 All
bmp_write results/visualizations/Slice.bmp Binary
wave 10
PlotMetedata results/visualizations/Metadata.png
k_means 10
bmp_write results/visualizations/kmeans.bmp Binary
triangulate
PlotVoronoi results/visualizations/Diagramma_Voronova.png
PlotDelaunay results/visualizations/Triangulation_Delone.png
find_path 20 27 100 298
PlotPath results/visualizations/Path.png
Plot3DPath results/visualizations/Plot3DPath.png
plotInteractive3DPath
end
