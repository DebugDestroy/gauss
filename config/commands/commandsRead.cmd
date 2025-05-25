init
help
bmp_read results/visualizations/Read.bmp
gnuplot results/visualizations/gnuplot.png
bmp_write results/visualizations/Polew.bmp Full
bin 147 All
bmp_write results/visualizations/Slice.bmp Binary
wave 10
PlotMetedata results/visualizations/Metadata.png
k_means 5
bmp_write results/visualizations/kmeans.bmp Binary
triangulate
PlotVoronoi results/visualizations/Diagramma_Voronova.png
PlotDelaunay results/visualizations/Triangulation_Delone.png
find_path
PlotPath results/visualizations/Path.png
end
