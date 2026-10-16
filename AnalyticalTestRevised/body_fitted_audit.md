# Audit body-fitted plastico

Build esaminata: KratosPR. Le macro ANALYTICAL_TANGENT risultano nuovamente attive in SolidElement e SupportSolidCondition; la libreria installata e successiva alle modifiche sorgenti.

## Controlli verificati

- Processi e modelers di tutti i quattro JSON body-fitted plastici coincidono con build_parameters.build. I tempi sono personalizzati: polynomial termina a 0.002 con passo 0.001; cosine_middle e cosine_boundary a 0.0025 con passo 0.00125. Non sono la stessa storia di carico.
- Tangente AnalyticalConstitutiveMatrix verificata rispetto a differenze finite di AnalyticalStress su tre stati plastici (normali e taglio). Massima discrepanza assoluta circa 3.6e-6 su moduli dell'ordine di 67000. Le tensioni coincidono con la proiezione radiale usata in analytic_fields, per nu=0.
- Jacobiano locale di SupportSolidCondition verificato rispetto a -dRHS/du su un punto del bordo inferiore, campo polynomial, t=0.002, mesh 4x4. Discrepanza relativa p=1: 2.73e-11; p=2: 6.76e-12. Questo test locale non individua un errore di linearizzazione di mDBOld; non dimostra stabilita globale.
- Prove p=2, mesh 8x8, 4 thread, risposta analitica attiva, senza leggere la deformazione plastica equivalente della legge con memoria: polynomial L2 standard 0.113310, quadratura 6 punti per span 0.113324; cosine_boundary L2 standard 0.107448, quadratura 6 punti 0.107438. Nessun avviso globale di non convergenza in queste quattro prove. L'errore elevato resta anche aumentando la quadratura.

## Problemi reali da separare

1. Risposta ibrida: assemblaggio/stress esportati usano la proiezione analitica, ma ApplyConstitutiveLaw e FinalizeSolutionStep continuano a calcolare e aggiornare la legge J2 con memoria. La equivalent_plastic_strain letta dalla legge appartiene quindi a una risposta che non governa il residuo analitico. Scegliere esplicitamente il modello prima di confrontare deformazioni plastiche o interpretare i warning costitutivi. Non disabilitare indiscriminatamente la macro se si vuole studiare la risposta analitica.
2. In analysis_stage.py RunSolutionLoop legge is_converged ma procede comunque a FinalizeSolutionStep. Il runner puo produrre CSV e aggiornare lo stato materiale dopo un solve fallito. Serve interrompere prima della finalizzazione quando SolveSolutionStep restituisce False.
3. Formula analitica attiva anche per materiali elastici: AnalyticalStress/AnalyticalConstitutiveMatrix leggono YIELD_STRESS con GetValue non const. Se assente, questo inserisce una variabile nelle proprieta condivise durante l'assemblaggio parallelo. Il rischio di race identificato resta quando si riattiva globalmente la macro. Serve una selezione esplicita e accessi const/prevalidati, non un valore di snervamento artificiale aggiunto all'elasticita.
4. I due passi di carico attraversano presto lo snervamento. mDBOld e congelato nel passo, e viene ricostruito nel successivo: confrontare storie temporali identiche e piu risolte per isolare l'effetto. Ridurre il passo non ha risolto la divergenza nelle prove precedenti con la risposta costitutiva effettiva.
5. Plasticita perfetta al bordo: perdita di rigidezza tangente in alcune direzioni e stabilita della formulazione penalty-free sono ipotesi da verificare. PENALTY_FACTOR=-1 sceglie Nitsche penalty-free; non e un penalty numericamente piccolo. Un confronto controllato con Dirichlet forte o una penalizzazione dimensionata puo isolare il contributo del bordo. Il Jacobiano locale corretto non garantisce coercivita/stabilita globale.

## Limiti

Le prove sono locali e su mesh piccole. Non isolano ancora la causa dell'errore globale nei casi p=2. Non sono stati cambiati materiali, JSON, macro C++ o la formulazione Nitsche durante questo audit. parameters_used.json e scritto da run_and_post, non documenta automaticamente i dati di convergence.csv; non assumere che i due file appartengano alla stessa esecuzione.

## Rettifica sulla prova di quadratura

Il controllo successivo ha trovato che NurbsSurfaceGeometry::CreateIntegrationPoints ignorava IntegrationInfo e usava sempre degree+1 per direzione nella build KratosPR. Le precedenti prove nominalmente a 6 punti non variavano quindi la quadratura del bulk; i numeri riportati non dimostrano che la quadratura sia irrilevante. La chiave corretta del modeler e number_of_integration_points_per_span; la correzione della geometria deve essere ricompilata e verificata contando i punti effettivi.
