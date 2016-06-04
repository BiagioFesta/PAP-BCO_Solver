# easy_parallel branch
#    generate 100 without compressione
#    O2
clc;
clear all;
close all;

samples_time_ms = [
1358 
1271 
1356 
1331 
1291 
1316 
1303 
1244 
1306 
1306 ];

samples_solution = [
86
85
86
85
86
85
85
85
86
86];

samples_num_edges = [
2406
2460
2438
2543
2503
2464
2437
2480
2483
2564 ];


mean(samples_num_edges)
mean(samples_time_ms)
mean(samples_solution)
max(samples_solution)
min(samples_solution)

