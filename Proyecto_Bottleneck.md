### Proyecto Bottleneck

En: /home/carlos/zrun_bottleneck

En yersin:
scp carlos@koch.ibv.csic.es:/home/carlos/zrun_bottleneck/*.DR.snp.final .

Es una forma que nosotsotros pensamos

y es esto y esto

Nosotros elegimos 4 baps inferiores, y los elegimos por estos criterios.

Presentacion pequena de 3 a 4 diapositivas
Grafisco, cuales elegi, explicar variables
 y por que los elegi
Hacer obbs ratio de variables
Explicar numeors de clusters

Muchos cluster s de pocos y pocos clusters de muchos

Fran descargo asi: Francisco José Martínez Martínez  6:31 PM
/data/fmartinez/MALAWI_DATASET_2020/descarga_dataset.py
Problemita de coverage
export LD_LIBRARY_PATH='/data/Software/miniconda3/lib/':"$LD_LIBRARY_PATH"


Revisar el municipio de residencia, codigos postales

Mycroreact XD, comentar

Pipeline: mirar en el tablet posiciones i es que hay posicone sconflictivas

Usar:
$FILE.DR.snp.final $FILE.snp.indel

~~~ Bash

[carlos@Koch zrun_bottleneck]$ cat *.meancov

G0076BottleNeck_mean_depth	264.944949056
G0076BottleNeck_median_depth	274.0
G0076BottleNeck_genome_coverage	0.979959569601
G76_mean_depth	56.5214263435
G76_median_depth	58.0
G76_genome_coverage	0.967151547354


G1180BottleNeck_mean_depth	362.549403019
G1180BottleNeck_median_depth	373.0
G1180BottleNeck_genome_coverage	0.98158961558
G1180_mean_depth	158.66146341
G1180_median_depth	165.0
G1180_genome_coverage	0.977134927277


G1181BottleNeck_mean_depth	208.757102975
G1181BottleNeck_median_depth	213.0
G1181BottleNeck_genome_coverage	0.979491478244
G1181_mean_depth	118.291111342
G1181_median_depth	122.0
G1181_genome_coverage	0.975182317617

~~~
> filter(Cluster04_Deep, SNPs==2)
  Position Freq_D Freq_R SNPs
1  2630158   9.00   0.00    2
2  2630161   8.97   0.00    2
3  3691003 100.00   0.00    2
4  2631968   0.00  10.20    2
5  2631971   0.00  10.61    2
6  4120926   0.00  42.51    2
> filter(Cluster10_Deep, SNPs==2)
  Position Freq_D Freq_R SNPs
1   711850  14.63   0.00    2
2  4120926  30.93   0.00    2
3  1404375   0.00  99.56    2
4  4364688   0.00 100.00    2
