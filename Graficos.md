### Grafico Sp_inclusters vs Sp_nocluster
>Descartar genotipos 11, 12, 16, 13

![](assets/Proyecto_BAPS-12da5bf8.jpeg)
- **Sp_incluster** =  Numero de españoles en clusters de transmicion
- **Sp_nocluster** =  Numero de españoles que no se encuentran en clusters de transmición

|Genotipo|N  |Sp |No_sp|Sp_x  |Sp_incluster|Sp_nocluster|Sp_incluster_x|N_incluster|
|--------|---|---|-----|------|------------|------------|--------------|-----------|
|1       |109|67 |42   |0.6147|39          |28          |0.5821        |59         |
|2       |218|120|98   |0.5505|69          |51          |0.575         |118        |
|3       |66 |33 |33   |0.5   |25          |8           |0.7576        |43         |
|4       |89 |58 |31   |0.6517|19          |39          |0.3276        |29         |
|5       |66 |51 |15   |0.7727|26          |25          |0.5098        |31         |
|6       |98 |44 |54   |0.449 |34          |10          |0.7727        |61         |
|7       |79 |55 |24   |0.6962|38          |17          |0.6909        |49         |
|8       |151|110|41   |0.7285|51          |59          |0.4636        |69         |
|9       |75 |55 |20   |0.7333|26          |29          |0.4727        |33         |
|10      |96 |46 |50   |0.4792|31          |15          |0.6739        |54         |
|11      |21 |4  |17   |0.1905|1           |3           |0.25          |10         |
|12      |2  |1  |1    |0.5   |1           |0           |1             |2          |
|13      |15 |5  |10   |0.3333|4           |1           |0.8           |7          |
|14      |35 |17 |18   |0.4857|8           |9           |0.4706        |16         |
|15      |30 |19 |11   |0.6333|6           |13          |0.3158        |8          |
|16      |27 |3  |24   |0.1111|1           |2           |0.3333        |11         |

- **Genotipo** = Genotipo
- **N** = Total de aislados
- **Sp** = Total de aislados españoles
- **No_sp** = Total de aislados NO españoles
- **Sp_x** = Proporcion de españoles en el genotipo
- **Sp_incluster** =  Numero de españoles en clusters de transmicion
- **Sp_nocluster** =  Numero de españoles que no se encuentran en clusters de transmición
- **Sp_incluster_x** = Proporcion de **"Sp_incluster"** en relacion a **"Sp"**
- **N_incluster** = Total de aislados en clusters de transmición

### Grafico Spanish_cluster_cases vs Foreign_cluster_cases
![](assets/Proyecto_BAPS-8ab2c7fe)
- **Spanish_cluster_cases** = Numero de aislados españoles en clusters de solo españoles
- **Foreign_cluster_cases** = Numero de aislados NO españoles en clusters de solo NO españoles

> Se calcula el **odds_ratio** basado en **Spanish_cluster_cases** vs **Foreign_cluster_cases**

> Genotipos 5 y 7 no poseen Foreign_cluster_cases, pero si Spanish_cluster_cases 15 y 18 casos.

> Descartar genotipo 15

|Genotipo|N_incluster|N  |Sp_incluster|Spanish_cluster_cases|Foreign_cluster_cases|Mixtos_cluster_cases|odds_ratio|conf.low |conf.high|pvalue   |
|--------|-----------|---|------------|---------------------|---------------------|--------------------|----------|---------|---------|---------|
|1       |59         |109|39          |7                    |2                    |50                  |2.7516    |0.4543   |30.0101  |0.2825   |
|2       |118        |218|69          |25                   |20                   |73                  |reference |reference|reference|reference|
|3       |43         |66 |25          |15                   |10                   |18                  |1.1969    |0.3998   |3.6808   |0.8036   |
|4       |29         |89 |19          |10                   |4                    |15                  |1.9774    |0.4773   |9.9632   |0.3611   |
|5       |31         |66 |26          |15                   |0                    |16                  |Inf       |2.4274   |Inf      |0.0011   |
|6       |61         |98 |34          |4                    |11                   |46                  |0.2969    |0.0597   |1.1978   |0.0747   |
|7       |49         |79 |38          |18                   |0                    |31                  |Inf       |2.9652   |Inf      |0.0003   |
|8       |69         |151|51          |35                   |8                    |26                  |3.4491    |1.2185   |10.6003  |0.0120   |
|9       |33         |75 |26          |12                   |2                    |19                  |4.6893    |0.8850   |47.9528  |0.0585   |
|10      |54         |96 |31          |15                   |9                    |30                  |1.3278    |0.4343   |4.2210   |0.6178   |
|11      |10         |21 |1           |0                    |8                    |2                   |0.0000    |0.0000   |0.5478   |0.0047   |
|12      |2          |2  |1           |0                    |0                    |2                   |0.0000    |0.0000   |Inf      |1.0000   |
|13      |7          |15 |4           |0                    |0                    |7                   |0.0000    |0.0000   |Inf      |1.0000   |
|14      |16         |35 |8           |3                    |3                    |10                  |0.8035    |0.0968   |6.6676   |1.0000   |
|15      |8          |30 |6           |2                    |0                    |6                   |Inf       |0.1392   |Inf      |0.5005   |
|16      |11         |27 |1           |0                    |5                    |6                   |0.0000    |0.0000   |1.0012   |0.0502   |

- **Genotipo** = Genotipo
- **N_incluster** = Total de aislados en clusters de transmición
- **N** = Total de aislados
- **Sp_incluster** =  Numero de españoles en clusters de transmicion
- **Spanish_cluster_cases** = Numero de aislados españoles en clusters de solo españoles
- **Foreign_cluster_cases** = Numero de aislados NO españoles en clusters de solo NO españoles
- **Mixtos_cluster_cases** = Numero de aislados que pertenecen a clusters mixtos

### OR y Grafico de "genotipos españoles" vs "genotipos no españoles" en nuestra poblacion
> basado en HIV Infection Disrupts the Sympatric Host–Pathogen Relationship in Human Tuberculosis

![](assets/Proyecto_BAPS-db75662f.png)

Dividir por cuadrantes agrupar
- Genotipos Españoles 5 9 7 8 15 4 1
- Genotipos mixtos 2 3 10 6 14 12
- Genotipos no Españoles 13 11 16

### Media y mediana de clusters españoles puros
- Para cada uno de los Genotipos
- Cual es el tamaño medio del cluster media y mediana
- Tamaño medio con percentiles plot

### Buscamos responder la pregunta si la relación es biología o azar
