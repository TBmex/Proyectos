### Tabla de ".snp_snps_result"

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

- **VH13894 Buscar resistencias a S, H, E.**
- VH13974 OK (Coincide)
- VH14235 OK (Se detecta resitencia a Et por WGS)
- VH14346 OK (Sencible)
- **VH14529 Buscar resistencia a R**
- VH14558 OK (Sencible)
- **VH14822 Buscar resistencia a Et**
- VH14825 OK (Se detecta resitencia a Fq por WGS)
- VH14832 OK (Sencible)
- VH15161 OK (Se detecta resistencia a H, E, Et y prothionamide por WGS)

### Busqueda en .snp por muestra

- **VH13894 Buscar resistencias a S, H, E.**
Nada en inhA y fabG1. En katG posicion no asociada a resistencia.
![](assets/all.snp_snps_result-806fbf91.png)

- **VH14529 Buscar resistencia a R**, Posicion filogenetica en rpoB
![](assets/all.snp_snps_result-8c974a6f.png)

- **VH14822 Buscar resistencia a Et**, Cambio cercano a fabG1
![](assets/all.snp_snps_result-a77aeffa.png)
ethA y ethR: Nada
