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
find_path 30 30 250 240
PlotPath results/visualizations/Path.png
Plot3DPath results/visualizations/Plot3DPath.png
plotInteractive3DPath
end
