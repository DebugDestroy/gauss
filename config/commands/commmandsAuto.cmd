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

k_means 10
PlotKmeans results/visualizations/kmeans.png
k_means_kern 10
PlotKmeans results/visualizations/kmeans_with_kernel.png

triangulate
PlotDelaunay results/visualizations/Triangulation_Delone.png
voronoi
PlotVoronoi results/visualizations/Diagramma_Voronova.png
build_nav_graph
PlotGraph results/visualizations/Graph.png
connect_to_graph
PlotGraph results/visualizations/ConnectedGraph.png

find_path_astar
save_metrics
PlotPath results/visualizations/AstarPath.png
Plot3DPath results/visualizations/AstarPlot3DPath.png
plotInteractive3DPath
find_path_dekstra
save_metrics
PlotPath results/visualizations/DekstraPath.png
Plot3DPath results/visualizations/DekstraPlot3DPath.png
plotInteractive3DPath
find_path_greedy
save_metrics
PlotPath results/visualizations/GreedyPath.png
Plot3DPath results/visualizations/GreedyPlot3DPath.png
plotInteractive3DPath
end
