## **New_Frecuencias**
### Tabla de New_Frecuencias de los 16 genotipos
|Genotipo|N  |Sp |For|NA_|Sp_x |For_x|NA_x |N_incluster|Sp_incluster|For_incluster|NA_incluster|Sp_incluster_x|For_incluster_x|NA_incluster_x|
|--------|---|---|---|---|-----|-----|-----|-----------|------------|-------------|------------|--------------|---------------|--------------|
|1       |109|69 |36 |4  |0.633|0.330|0.037|59         |41          |15           |3           |0.695         |0.254          |0.051         |
|2       |218|132|70 |16 |0.606|0.321|0.073|118        |73          |32           |13          |0.619         |0.271          |0.110         |
|3       |66 |33 |31 |2  |0.500|0.470|0.030|43         |25          |17           |1           |0.581         |0.395          |0.023         |
|4       |89 |61 |23 |5  |0.685|0.258|0.056|29         |20          |6            |3           |0.690         |0.207          |0.103         |
|5       |66 |52 |8  |6  |0.788|0.121|0.091|31         |27          |1            |3           |0.871         |0.032          |0.097         |
|6       |98 |47 |44 |7  |0.480|0.449|0.071|61         |37          |17           |7           |0.607         |0.279          |0.115         |
|7       |79 |56 |20 |3  |0.709|0.253|0.038|49         |38          |11           |0           |0.776         |0.224          |0.000         |
|8       |151|110|38 |3  |0.728|0.252|0.020|69         |51          |15           |3           |0.739         |0.217          |0.043         |
|9       |75 |55 |15 |5  |0.733|0.200|0.067|33         |27          |4            |2           |0.818         |0.121          |0.061         |
|10      |96 |46 |41 |9  |0.479|0.427|0.094|54         |32          |16           |6           |0.593         |0.296          |0.111         |
|11      |21 |5  |14 |2  |0.238|0.667|0.095|10         |2           |7            |1           |0.200         |0.700          |0.100         |
|12      |2  |1  |0  |1  |0.500|0.000|0.500|2          |1           |0            |1           |0.500         |0.000          |0.500         |
|13      |15 |5  |10 |0  |0.333|0.667|0.000|7          |4           |3            |0           |0.571         |0.429          |0.000         |
|14      |35 |17 |14 |4  |0.486|0.400|0.114|16         |7           |5            |4           |0.438         |0.312          |0.250         |
|15      |30 |20 |7  |3  |0.667|0.233|0.100|8          |6           |0            |2           |0.750         |0.000          |0.250         |
|16      |27 |3  |21 |3  |0.111|0.778|0.111|11         |1           |8            |2           |0.091         |0.727          |0.182         |

- Genotipo = Genotipo
- N = Total de aislados
- Sp = Total de aislados españoles
- For = Total de aislados extranjeros
- NA_ = Total de aislados NA
  - Sp_x = Proporcion de españoles en el genotipo
  - For_x = Proporcion de extranjeros en el genotipo
  - NA_x = Proporcion de NA en el genotipo
- N_incluster = Total de aislados en clusters de transmisión
  - Sp_incluster = Numero de españoles en clusters de transmisión
  - For_incluster = Numero de extranjeros en clusters de transmisión
  - NA_incluster = Numero de NA en clusters de transmisión
  - Sp_incluster_x = Proporción de "Sp_incluster" en relacion a "N_incluster"
  - For_incluster_x = Proporción de "For_incluster" en relacion a "N_incluster"
  - NA_incluster_x = Proporción de "NA_incluster" en relacion a "N_incluster"

### Genero un subset al descartar genotipos por diferentes razones
> - N < 25 (12, 13, 11)
- Na_x > 0.1 (12, 14, 16, 15)
- NA_incluster_x > 0.15 (12,14,15,16)
- Numero de clusters < 6 (11, 12, 13, 14, 15, 16)

|Genotipo|N  |Sp |For|NA_|Sp_x |For_x|NA_x |N_incluster|Sp_incluster|For_incluster|NA_incluster|Sp_incluster_x|For_incluster_x|NA_incluster_x|
|--------|---|---|---|---|-----|-----|-----|-----------|------------|-------------|------------|--------------|---------------|--------------|
|1       |109|69 |36 |4  |0.633|0.330|0.037|59         |41          |15           |3           |0.695         |0.254          |0.051         |
|2       |218|132|70 |16 |0.606|0.321|0.073|118        |73          |32           |13          |0.619         |0.271          |0.110         |
|3       |66 |33 |31 |2  |0.500|0.470|0.030|43         |25          |17           |1           |0.581         |0.395          |0.023         |
|4       |89 |61 |23 |5  |0.685|0.258|0.056|29         |20          |6            |3           |0.690         |0.207          |0.103         |
|5       |66 |52 |8  |6  |0.788|0.121|0.091|31         |27          |1            |3           |0.871         |0.032          |0.097         |
|6       |98 |47 |44 |7  |0.480|0.449|0.071|61         |37          |17           |7           |0.607         |0.279          |0.115         |
|7       |79 |56 |20 |3  |0.709|0.253|0.038|49         |38          |11           |0           |0.776         |0.224          |0.000         |
|8       |151|110|38 |3  |0.728|0.252|0.020|69         |51          |15           |3           |0.739         |0.217          |0.043         |
|9       |75 |55 |15 |5  |0.733|0.200|0.067|33         |27          |4            |2           |0.818         |0.121          |0.061         |
|10      |96 |46 |41 |9  |0.479|0.427|0.094|54         |32          |16           |6           |0.593         |0.296          |0.111         |


### Grafico Spanish in cluster (%) VS Foreing in cluster (%)

![](assets/New_Baps-951cd8c5.png)

- Spanish in cluster (%) = Sp_incluster_x = Proporción de "Sp_incluster" en relacion a "N_incluster"
- Foreing in cluster (%) = For_incluster_x = Proporción de "For_incluster" en relacion a "N_incluster"

> Agrupo genotipos > 70% Spanish in cluster

|Genotipos   |Sp_incluster|For_incluster|Odds |pvalue|
|------------|------------|-------------|-----|------|
|5,9,7,8     |143         |31           |2.081|0.0014|
|1,4,2,6,10,3|228         |103          |ref  |ref   |

### Grafico Spanish cluster cases (%) VS Foreing cluster cases (%)

![](assets/New_Baps-e8dabb1f.png)

### Tabla de clusters
|Genotipo|Mix|Sp |For|Na |Total|
|--------|---|---|---|---|-----|
|1       |9  |4  |1  |1  |15   |
|2       |13 |8  |3  |11 |35   |
|3       |4  |7  |3  |1  |15   |
|4       |2  |4  |2  |3  |11   |
|5       |0  |5  |0  |2  |7    |
|6       |3  |2  |3  |4  |12   |
|7       |7  |4  |0  |0  |11   |
|8       |5  |14 |1  |3  |23   |
|9       |2  |5  |1  |2  |10   |
|10      |4  |5  |2  |4  |15   |

- Genotipo = Genotipo
- Mix = Numero de clusters mixtos
- Sp = Numero de clusters de solo españoles
- For = Numero de clusters de solo extranjeros
- Na = Numero de clusters NA
- Total = Total de clusters

#### Eventos de transmisión en only spanish cluster
|Genotipo|Sp |Eventos|
|--------|---|-------|
|1       |4  |15     |
|2       |8  |22     |
|3       |7  |8      |
|4       |4  |9      |
|5       |5  |18     |
|6       |2  |2      |
|7       |4  |14     |
|8       |14 |21     |
|9       |5  |12     |
|10      |5  |10     |

- Genotipo = Genotipo
- Sp = Numero de clusters de solo españoles
- Eventos = Numero de eventos de transmisión

#### Tabla_only_cases_cluster

|Genotipo|N_incluster|Spanish_cluster_cases|Foreing_cluster_cases|Mixes_cluster_cases|Na_cluster_cases|
|--------|-----------|---------------------|---------------------|-------------------|----------------|
|1       |59         |19                   |2                    |32                 |6               |
|2       |118        |30                   |7                    |37                 |44              |
|3       |43         |15                   |8                    |18                 |2               |
|4       |29         |13                   |4                    |5                  |7               |
|5       |31         |23                   |0                    |0                  |8               |
|6       |61         |4                    |9                    |12                 |36              |
|7       |49         |18                   |0                    |31                 |0               |
|8       |69         |35                   |3                    |24                 |7               |
|9       |33         |17                   |2                    |6                  |8               |
|10      |54         |15                   |9                    |9                  |21              |

- Genotipo = Genotipo
- N_incluster = Total de casos en clusters de transmición
- Spanish_cluster_cases = Numero de casos españoles en clusters de solo españoles
- Foreing_cluster_cases = Numero de casos extranjeros en clusters de solo extranjeros
- Mixes_cluster_cases = Numero de casos en clusters mixtos
- Na_cluster_cases = Numero de casos en clusters NA
