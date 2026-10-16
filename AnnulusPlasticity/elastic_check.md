# Controllo SBM in campo elastico

Verifica del 16 ottobre 2026 con la libreria Kratos installata aggiornata.

**Newton converge anche per p2 e la risposta materiale è elastica corretta.
La convergenza spaziale dell'annulus misto, però, non raggiunge gli ordini L2
attesi sulle mesh verificate. Il problema di accuratezza non è esclusivo
della plasticità.**

## Configurazione e prove

Sono stati eseguiti 19 solve, tutti completati con Newton–Raphson puro,
senza line search. I file di produzione e i sorgenti Kratos non sono stati
modificati. Le copie dei parametri e dei materiali sono nelle cartelle dei risultati.

- SBM mixed: Dirichlet esterna, Neumann interna, p1/p2, mesh 10/20/40/80.
- Materiale J2 corrente: E=67000, nu=0, yield=201, H=10000,
  `PENALTY_FACTOR=-1`.
- Parametri geometrici mixed correnti: lambda inner/outer=1,
  proiezione `linealized`, quadratura bulk 4×4.
- Annulus: arresto a t=0.0002, due step da 0.0001, forcing nullo.
  La tensione equivalente analitica massima sul dominio fisico è 45.225,
  ampiamente inferiore a 201; la soluzione analitica non contiene plasticità.
- Patch affine: u=(t*x,0), arresto a t=0.002, due step da 0.001,
  forcing nullo, trazione specifica del campo affine. La tensione equivalente
  analitica è 134, ancora inferiore a 201.

Comando riproducibile per i 16 controlli mixed:

```bash
python3 AnnulusPlasticity/diagnose_mixed.py \
  --cases annulus_elastic patch_elastic --degrees 1 2 \
  --meshes 10 20 40 80 --strategy newton_raphson \
  --output AnnulusPlasticity/diagnostics/elastic_check
```

## Risposta materiale e patch affine

In tutti i 19 solve la deformazione plastica equivalente campionata nei
punti di integrazione degli elementi è **esattamente zero**. Anche le tensioni
numeriche equivalenti massime campionate rimangono inferiori a 201.
La tensione restituita dalla legge coincide con la tensione elastica ricavata
dal gradiente numerico, entro errore relativo inferiore a 1e-11.

Il patch affine mixed dà errori L2 relativi compresi fra 4.06e-7 e 5.65e-7;
gli errori relativi di stress sono inferiori a 1.7e-6. Quindi p2 non fallisce
in modo generale nel caso elastico e riproduce praticamente il campo affine.
La piccola discrepanza residua non è una dimostrazione di esattezza geometrica
della discretizzazione della skin e delle sue normali.

## Annulus non affine, condizioni miste

Errori relativi L2 dello spostamento:

| Spans per direzione | p1 | p2 |
|---:|---:|---:|
| 10 | 0.244733 | 0.0817539 |
| 20 | 0.176098 | 0.0252796 |
| 40 | 0.101033 | 0.0143229 |
| 80 | 0.0523124 | 0.00500414 |

Gli ordini L2 misurati fra mesh 40 e 80 sono **0.95 per p1** e **1.52 per p2**,
inferiori rispettivamente a p+1=2 e p+1=3. Gli ordini H1 bulk sullo stesso
raffinamento sono 1.05 e 1.63. Gli errori relativi di stress alla mesh 80 sono
5.95% per p1 e 0.664% per p2.

La serie mostra riduzione dell'errore, ma non giustifica affermare che la
configurazione mixed abbia già convergenza spaziale ottimale. Non è stata
verificata una regione asintotica oltre la mesh 80.

## Controllo con soli Dirichlet e stessa geometria

Per separare il ruolo della Neumann ho copiato i parametri **mixed**, mantenuto
lambda inner/outer=1 e tutti gli altri parametri, e cambiato soltanto la
condition importata sulla skin interna da `SbmLoadSolidCondition` a
`SbmSolidCondition`. Anche sul cerchio interno viene così imposto lo spostamento
analitico, anziché la trazione. Non è stato usato direttamente il JSON Dirichlet
di produzione, che ha lambda=0.5 e avrebbe cambiato la geometria surrogata.

| Spans | p2 mixed, L2 relativo | p2 soli Dirichlet, L2 relativo |
|---:|---:|---:|
| 20 | 0.0252796 | 0.00206594 |
| 40 | 0.0143229 | 0.000457107 |
| 80 | 0.00500414 | 0.000108148 |

Alla mesh 80 l'errore mixed è circa **46 volte maggiore**. I soli Dirichlet
mostrano ordine L2 2.08 fra 40 e 80: migliorano molto l'accuratezza, ma neppure
questa serie dimostra l'ordine 3 atteso per p2. Il confronto indica che il
trattamento Neumann aumenta sensibilmente l'errore; non isola da solo quale
termine della formulazione o dell'implementazione ne sia responsabile.

## Artefatti

- `diagnostics/elastic_check/summary.json`: esiti dei 16 solve mixed.
- `diagnostics/elastic_check/accuracy.json`: metriche aggiuntive e ordini per
  tutti i 19 solve; errore di stress con norma tensoriale, peso 2 sul taglio.
- `diagnostics/elastic_check/elastic_convergence.png`: curve L2 e riferimenti
  h^p, h^(p+1).
- `diagnostics/elastic_dirichlet_control/summary.json`: tre controlli Dirichlet;
  i rispettivi `parameters.json` conservano la configurazione completa.

Il controllo conferma il funzionamento del solver e della risposta elastica
locale. Per capire il difetto di accuratezza globale conviene proseguire
dal caso elastico non affine, verificando consistenza dei termini Neumann,
estensione delle funzioni di base e trattamento geometrico/Dirichlet.
