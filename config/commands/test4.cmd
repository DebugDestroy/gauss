init

g 1 1 2.84496 2.03156 64.953
g 15 15 3.7345 2.74394 -29.815
g 1 15 3.4613 4.53558 18.4512

generate

bin 132 All
bmp_write results/visualizations/Slice.bmp Binary

wave 10
PlotMetedata results/visualizations/Metadata.png

triangulate
PlotDelaunay results/visualizations/Triangulation_Delone.png
voronoi
PlotVoronoi results/visualizations/Diagramma_Voronova.png
build_nav_graph
PlotGraph results/visualizations/Graph.png

end
