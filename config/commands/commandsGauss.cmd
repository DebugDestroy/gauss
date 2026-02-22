init

g 4 100 5 5 100
g 50 50 1.5 1.5 40 
g 30 40 4.5 4.5 140 
g 50 150 4.0 4.0 105 
g 75 275 3.0 3.0 105 
g 75 28 3.0 3.0 105 
g 75 15 3.0 3.0 105 
g 125 90 3.0 3.0 105 
g 275 250 3.0 3.0 105 
g 75 285 3.0 3.0 105 
g 90 50 5 5 100
g 150 250 5 5 100
g 150 110 4.5 4.5 140
g 150 5 3.0 3.0 105 
g 100 150 4.5 4.5 140
g 250 50 4.0 4.0 -110  
g 50 250 3.0 3.0 105  
g 280 200 3.0 3.0 105  
g 225 150 3.0 3.0 105 
g 240 200 3.0 3.0 105 
g 200 250 5 5 100
g 250 250 7.8 10.8 -105  
g 150 150 5.0 5.0 100  
g 100 200 4.5 4.5 -95  
g 200 100 4.2 4.2 90  
g 70 180 3.7 3.7 -85  
g 180 70 3.9 3.9 80  
g 120 120 4.8 4.8 -75  
g 220 220 4.3 4.3 70  
g 30 120 3.2 3.2 -65  
g 120 30 3.4 3.4 60  
g 270 180 3.6 3.6 -55  
g 180 270 3.3 3.3 50  
g 80 80 4.1 4.1 -45  
g 220 80 4.4 4.4 40  
g 80 220 3.8 3.8 -35  
g 170 170 5.5 5.5 30  
g 130 130 4.6 4.6 -25  
g 200 200 4.9 4.9 20  
g 60 60 3.5 3.5 -15  
g 240 240 4.2 4.2 10  
g 90 90 4.0 4.0 -5  
g 210 210 4.7 4.7 5  

generate
gnuplot results/visualizations/gnuplot.png
bmp_write results/visualizations/Pole.bmp Full
bin 132 All
bmp_write results/visualizations/Slice.bmp Binary
wave 10
PlotMetedata results/visualizations/Metadata.png
k_means 10
bmp_write results/visualizations/kmeans.bmp Binary
k_means_kern 5
bmp_write results/visualizations/kmeans_with_kernel.bmp Binary
triangulate
PlotVoronoi results/visualizations/Diagramma_Voronova.png
PlotDelaunay results/visualizations/Triangulation_Delone.png
find_path_astar 5 5 295 295
PlotPath results/visualizations/AstarPath.png
Plot3DPath results/visualizations/AstarPlot3DPath.png
plotInteractive3DPath
find_path_dekstra 5 5 295 295
PlotPath results/visualizations/DekstraPath.png
Plot3DPath results/visualizations/DekstraPlot3DPath.png
plotInteractive3DPath
find_path_greedy 5 5 295 295
PlotPath results/visualizations/GreedyPath.png
Plot3DPath results/visualizations/GreedyPlot3DPath.png
plotInteractive3DPath
end
