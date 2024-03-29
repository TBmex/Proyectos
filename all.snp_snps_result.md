### Antibiograma

![](assets/all.snp_snps_result-0f6af92c.png)

### Resistencias WGS

![](assets/all.snp_snps_result-2875bfab.png)

### Buscamos si hay posiciones de resisetencia en los archivos con terminación ".snp_snps_result", estos archivos contienen snps asociados a resistencia de cada muestra.

Nota: No se encontró SNP que explique las discrepancias en la tabla.
Nota: Los SNPs estan fijados.

|Sample     |Position|Gene      |Antibiotic                                           |Freq |N_reads|WT |MUT|Codon_change|AA_change|PhyResSE       |ReSeqTB            |Comments                                    |
|-----------|--------|----------|-----------------------------------------------------|-----|-------|---|---|------------|---------|---------------|-------------------|--------------------------------------------|
|VH13894.snp|761139  |rpoB      |rifampicin (RMP)                                     |100  |63     |C  |G  |cac/gac     |His445Asp|High_confidence|High_confidence    |Check for double mutations in the same codon|
|VH13974.snp|2155168 |katG      |isoniazid (INH)                                      |100  |97     |C  |G  |agc/acc     |Ser315Thr|High_confidence|High_confidence    |                                            |
|VH14235.snp|1673425 |inhA-fabG1|isoniazid (INH) ethionamide (ETH) prothionamide (PTO)|100  |85     |C  |T  |-           |c-15t    |High_confidence|Moderate_confidence|                                            |
|VH14235.snp|2726136 |ahpC      |isoniazid (INH)                                      |100  |90     |C  |T  |-           |---      |Low_confidence |-                  |                                            |
|VH14529.snp|1673425 |inhA-fabG1|isoniazid (INH) ethionamide (ETH) prothionamide (PTO)|100  |62     |C  |T  |-           |c-15t    |High_confidence|Moderate_confidence|                                            |
|VH14822.snp|761155  |rpoB      |rifampicin (RMP)                                     |100  |89     |C  |T  |tcg/ttg     |Ser450Leu|High_confidence|High_confidence    |Check for double mutations in the same codon|
|VH14822.snp|765461  |rpoC      |rifampicin (RMP)                                     |100  |105    |A  |C  |aac/cac     |Asn698His|-              |-                  |Rifampicin compensatory mutation            |
|VH14822.snp|2155168 |katG      |isoniazid (INH)                                      |98,44|63     |C  |G  |agc/acc     |Ser315Thr|High_confidence|High_confidence    |                                            |
|VH14822.snp|2288859 |pncA      |pyrazinamide (PZA)                                   |99,07|106    |A  |C  |gtc/ggc     |Val128Gly|Low_confidence |High_confidence    |                                            |
|VH14825.snp|7572    |gyrA      |fluoroquinolones (FQ)                                |100  |65     |T  |C  |tcg/ccg     |Ser91Pro |High_confidence|High_confidence    |Resistance to Moxifloxacin                  |
|VH15161.snp|1673425 |inhA-fabG1|isoniazid (INH) ethionamide (ETH) prothionamide (PTO)|100  |73     |C  |T  |-           |c-15t    |High_confidence|Moderate_confidence|                                            |
|VH15161.snp|4247429 |embB      |ethambutol (EMB)                                     |100  |74     |A  |G  |atg/gtg     |Met306Val|High_confidence|-                  |                                            |

### Eligo las muestras que tienen resistencia fenotípica pero no se encontro resistencia genotípica.
- **VH13894 Buscar resistencias a S, H, E.**
- VH13974 OK (Coincide)
- VH14235 OK (Se detecta resistencia a Et por WGS)
- VH14346 OK (Sensible)
- **VH14529 Buscar resistencia a R**
- VH14558 OK (Sensible)
- **VH14822 Buscar resistencia a Et**
- VH14825 OK (Se detecta resistencia a Fq por WGS)
- VH14832 OK (Sensible)
- VH15161 OK (Se detecta resistencia a H, E, Et y prothionamide por WGS)

### Realice una búsqueda en el archivo ".snp" de las muestras seleccionadas.

- **VH13894 Buscar resistencias a S, H, E.**
  - gidB
    - 4407588	| synonymous_variant | 615G>A | A205A
  - rpsL, rrs: NA
  - inhA, fabG1, ahpC, ndh: NA
  - katG posicion NO asociada a resistencia.
    - 2154724 | missense_variant | 1388T>G | L463R
  - embA, embB, embC: NO asociadas a resistencia
    - 4240671	| embC | missense_variant | 809T>C  | I270T
    - 4242803	| embC | missense_variant | 2941G>C | V981L
    - 4247646	| embB | missense_variant | 1133C>A | A378E

- **VH14529 Buscar resistencia a R**
  - rpoB, posición filogenetica
    - 763031 | synonymous_variant | 3225C>T | A1075A
  - rpoC: NA

- **VH14822 Buscar resistencia a Et**
  - etoA, mshA, inhA: NA
  - ndhA
    - 471666 | missense_variant | 974T>C | M325T

### Conclusión: No se encontró explicación a las incongruencias entre la resistencia fenotípica y las genotípicas vistas en la tabla.
