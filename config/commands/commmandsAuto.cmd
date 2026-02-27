init

g_auto 

generate

save_g

gnuplot results/visualizations/gnuplot.png
bmp_write results/visualizations/Pole.bmp Full

bin 132 All
bmp_write results/visualizations/Slice.bmp Binary

wave 10
PlotMetedata results/visualizations/Metadata.png

k_means 100
bmp_write results/visualizations/kmeans.bmp Binary
k_means_kern 10
bmp_write results/visualizations/kmeans_with_kernel.bmp Binary

triangulate
PlotVoronoi results/visualizations/Diagramma_Voronova.png
PlotDelaunay results/visualizations/Triangulation_Delone.png

find_path_astar
PlotPath results/visualizations/AstarPath.png
Plot3DPath results/visualizations/AstarPlot3DPath.png
plotInteractive3DPath
find_path_dekstra
PlotPath results/visualizations/DekstraPath.png
Plot3DPath results/visualizations/DekstraPlot3DPath.png
plotInteractive3DPath
find_path_greedy
PlotPath results/visualizations/GreedyPath.png
Plot3DPath results/visualizations/GreedyPlot3DPath.png
plotInteractive3DPath
end
