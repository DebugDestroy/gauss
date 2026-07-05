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

astar_graph
PlotPath results/visualizations/AstarPathGraph.png
Plot3DPath results/visualizations/AstarPlot3DPathGraph.png
dekstra_graph
PlotPath results/visualizations/DekstraPathGraph.png
Plot3DPath results/visualizations/DekstraPlot3DPathGraph.png
greedy_graph
PlotPath results/visualizations/GreedyPathGraph.png
Plot3DPath results/visualizations/GreedyPlot3DPathGraph.png

grid
PlotGrid results/visualizations/Grid.png
build_nav_grid
PlotNavGrid results/visualizations/NavigationGrid.png
connect_to_grid
PlotNavGrid results/visualizations/ConnectedGrid.png

astar_grid
PlotGridPath results/visualizations/AstarPathGrid.png
dekstra_grid
PlotGridPath results/visualizations/DekstraPathGrid.png
greedy_grid
PlotGridPath results/visualizations/GreedyPathGrid.png
end
