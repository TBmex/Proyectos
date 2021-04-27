##  New Frecuencies
### Tabla de New_Frecuencias de los 16 genotipos
|Genotipo|N  |Sp |For|NA_|Sp_x |For_x|N_incluster|Sp_incluster|For_incluster|NA_incluster|Sp_incluster_x|For_incluster_x|NA_incluster_x|
|--------|---|---|---|---|-----|-----|-----------|------------|-------------|------------|--------------|---------------|--------------|
|1       |109|69 |36 |4  |0.633|0.330|59         |41          |15           |3           |0.695         |0.254          |0.051         |
|2       |218|132|70 |16 |0.606|0.321|118        |73          |32           |13          |0.619         |0.271          |0.110         |
|3       |66 |33 |31 |2  |0.500|0.470|43         |25          |17           |1           |0.581         |0.395          |0.023         |
|4       |89 |61 |23 |5  |0.685|0.258|29         |20          |6            |3           |0.690         |0.207          |0.103         |
|5       |66 |52 |8  |6  |0.788|0.121|31         |27          |1            |3           |0.871         |0.032          |0.097         |
|6       |98 |47 |44 |7  |0.480|0.449|61         |37          |17           |7           |0.607         |0.279          |0.115         |
|7       |79 |56 |20 |3  |0.709|0.253|49         |38          |11           |0           |0.776         |0.224          |0.000         |
|8       |151|110|38 |3  |0.728|0.252|69         |51          |15           |3           |0.739         |0.217          |0.043         |
|9       |75 |55 |15 |5  |0.733|0.200|33         |27          |4            |2           |0.818         |0.121          |0.061         |
|10      |96 |46 |41 |9  |0.479|0.427|54         |32          |16           |6           |0.593         |0.296          |0.111         |
|11      |21 |5  |14 |2  |0.238|0.667|10         |2           |7            |1           |0.200         |0.700          |0.100         |
|12      |2  |1  |0  |1  |0.500|0.000|2          |1           |0            |1           |0.500         |0.000          |0.500         |
|13      |15 |5  |10 |0  |0.333|0.667|7          |4           |3            |0           |0.571         |0.429          |0.000         |
|14      |35 |17 |14 |4  |0.486|0.400|16         |7           |5            |4           |0.438         |0.312          |0.250         |
|15      |30 |20 |7  |3  |0.667|0.233|8          |6           |0            |2           |0.750         |0.000          |0.250         |
|16      |27 |3  |21 |3  |0.111|0.778|11         |1           |8            |2           |0.091         |0.727          |0.182         |


- Genotipo = Genotipo
- N = Total de aislados
- Sp = Total de aislados españoles
- For = Total de aislados extranjeros
- NA_ = Total de aislados NA
  - Sp_x = Proporcion de españoles en el genotipo
  - For_x = Proporcion de extranjeros en el genotipo
- N_incluster = Total de aislados en clusters de transmisión
  - Sp_incluster = Numero de españoles en clusters de transmisión
  - For_incluster = Numero de extranjeros en clusters de transmisión
  - NA_incluster = Numero de NA en clusters de transmisión
  - Sp_incluster_x = Proporción de "Sp_incluster" en relacion a "N_incluster"
  - For_incluster_x = Proporción de "For_incluster" en relacion a "N_incluster"
  - NA_incluster_x = Proporción de "NA_incluster" en relacion a "N_incluster"

#### Genero un subset al descartar genotipos por diferentes razones
|Genotipo|N  |Sp |For|NA_|Sp_x |For_x|N_incluster|Sp_incluster|For_incluster|NA_incluster|Sp_incluster_x|For_incluster_x|NA_incluster_x|
|--------|---|---|---|---|-----|-----|-----------|------------|-------------|------------|--------------|---------------|--------------|
|1       |109|69 |36 |4  |0.633|0.330|59         |41          |15           |3           |0.695         |0.254          |0.051         |
|2       |218|132|70 |16 |0.606|0.321|118        |73          |32           |13          |0.619         |0.271          |0.110         |
|3       |66 |33 |31 |2  |0.500|0.470|43         |25          |17           |1           |0.581         |0.395          |0.023         |
|4       |89 |61 |23 |5  |0.685|0.258|29         |20          |6            |3           |0.690         |0.207          |0.103         |
|5       |66 |52 |8  |6  |0.788|0.121|31         |27          |1            |3           |0.871         |0.032          |0.097         |
|6       |98 |47 |44 |7  |0.480|0.449|61         |37          |17           |7           |0.607         |0.279          |0.115         |
|7       |79 |56 |20 |3  |0.709|0.253|49         |38          |11           |0           |0.776         |0.224          |0.000         |
|8       |151|110|38 |3  |0.728|0.252|69         |51          |15           |3           |0.739         |0.217          |0.043         |
|9       |75 |55 |15 |5  |0.733|0.200|33         |27          |4            |2           |0.818         |0.121          |0.061         |
|10      |96 |46 |41 |9  |0.479|0.427|54         |32          |16           |6           |0.593         |0.296          |0.111         |

#### Grafico las porporciones en cluster de españoles vs extranjeros

![](assets/New_Baps-d3d646aa)

- Spanish in cluster (%) = Proporción de "Sp_incluster" en relacion a "N_incluster"
- Foreing in cluster (%) = Proporción de "For_incluster" en relacion a "N_incluster"

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
|11      |2  |0  |0  |1  |3    |
|12      |0  |0  |0  |1  |1    |
|13      |3  |0  |0  |0  |3    |
|14      |1  |0  |0  |4  |5    |
|15      |0  |1  |0  |2  |3    |
|16      |0  |0  |1  |2  |3    |

#### Grafico los Only Spanish Clusters y cuento su numero de eventos de transmisión
![](assets/New_Baps-aa8a4ba6)
