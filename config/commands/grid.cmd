init 1000 1000

g_auto 900 901 10 990 10 990 3 6 3 6 -20 20

generate

bin 5

wave 10

grid 5
build_nav_grid 5
connect_to_grid 3 997 997 3

astar_grid
save_metrics var/metrics/grid.csv

dekstra_grid
save_metrics var/metrics/grid.csv

greedy_grid
save_metrics var/metrics/grid.csv

end
