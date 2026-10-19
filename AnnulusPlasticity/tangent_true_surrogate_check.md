# Confronto diretto D_alg true / surrogate

Configurazione corrente SBM mixed, p2, proiezione analytical, lambda=0.5,
H=10000, dt=0.001. Due mesh: 10×10 e 20×20, entrambe convergenti.

Le matrici sono state lette direttamente dalle due istanze costitutive della
`SbmLoadSolidCondition`, alla soluzione convergente di ciascuno step,
**prima di finalizzare entrambe le storie materiali**. Questa è la tangente
algoritmica dello step, non una risposta ricalcolata dopo il commit.
La strumentazione è stata compilata in una copia della libreria in `/tmp`;
non sono stati modificati i sorgenti o la libreria installata di Kratos.

Differenza relativa: `||D_true-D_surrogate||_F / ||D_surrogate||_F`.

| Mesh | Tempo | Punti Neumann diversi / totali | Differenza media | Massima |
|---|---:|---:|---:|---:|
| 10×10 | 0.001 | 80 / 80 | 16.56% | 58.98% |
| 10×10 | 0.002 | 80 / 80 | 1.87% | 5.11% |
| 20×20 | 0.001 | 200 / 200 | 8.67% | 22.32% |
| 20×20 | 0.002 | 200 / 200 | 1.17% | 2.71% |

Le matrici sono considerate diverse se la differenza relativa supera 1e-6.
L'output numerico della strumentazione usa sei cifre significative: la
precisione è sufficiente per le differenze percentuali riportate.

Al primo step della mesh 10, in 16 punti su 80 una matrice è elastica e
l'altra è plastica. Negli altri confronti entrambe sono plastiche, ma le
tangenti non coincidono. Quindi, in questo studio p2, sostituire D_true con
D_surrogate non è un'identità numerica.

Controllo separato sotto snervamento sulle stesse due mesh, con due step a
t=0.0001 e t=0.0002: tutte le coppie coincidono, come atteso per la matrice
elastica costante. Riepilogo in
`diagnostics/tangent_true_surrogate_elastic/comparison_summary.json`.

Tutte le coppie di matrici sono in
`diagnostics/tangent_true_surrogate/comparison.json`; riepilogo e coppie con
massima differenza in `comparison_summary.json` nella stessa cartella.
