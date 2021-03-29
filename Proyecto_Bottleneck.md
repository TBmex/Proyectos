### Proyecto Bottleneck



~~~r
[carlos@Koch zrun_bottleneck]$ cat *.meancov
G0076BottleNeck_mean_depth	264.944949056
G0076BottleNeck_median_depth	274.0
G0076BottleNeck_genome_coverage	0.979959569601
G1180BottleNeck_mean_depth	362.549403019
G1180BottleNeck_median_depth	373.0
G1180BottleNeck_genome_coverage	0.98158961558
G1181BottleNeck_mean_depth	208.757102975
G1181BottleNeck_median_depth	213.0
G1181BottleNeck_genome_coverage	0.979491478244
[carlos@Koch zrun_bottleneck]$ pwd
/home/carlos/zrun_bottleneck
~~~
Ganmos 15 posiciones con el deep secuenci
~~~r
> Cluster04_deep %>% anti_join(Cluster04_no_deep, by = "Position")
   Position Freq_D Freq_R SNPs
1     21819 100.00 100.00    1
2    584438 100.00 100.00    1
3   1160770 100.00 100.00    1
4   2631968   9.81  10.20    1
5   2631971  10.31  10.61    1
6   2631977  13.59  14.15    1
7   3336528  11.52  16.04    1
8   3663889  35.06  26.04    1
9   3691005 100.00 100.00    1
10  3820562  14.62  17.26    1
11  3820565  15.71  17.47    1
12  4060201 100.00 100.00    1
13  4060230 100.00 100.00    1
14  4060334 100.00 100.00    1
15  4361162 100.00 100.00    1
~~~
### Como obtenemos las variables identificadas
