init 1000 1000

g_auto 900 901 10 990 10 990 3 6 3 6 -20 20

rrt 20000 3 997 997 3 1 5 30 30 2 2 2 20 10 0.2
save_metrics var/metrics/rrt.csv

rrt_star 20000 3 997 997 3 1 5 30 30 2 2 2 20 40 700 10 0.02
save_metrics var/metrics/rrt.csv

end
