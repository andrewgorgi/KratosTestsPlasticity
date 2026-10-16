# Diagnosi della convergenza plastica

Caso: body_fitted, polynomial, plastic, p=1, 24 thread, KratosPR dopo la correzione della race.

Il JSON corrente usa end_time=0.0025 e time_step=0.00125. Riprodotto il messaggio `Maximum number of iterations in plasticity loop reached`: riguarda il return mapping costitutivo, non solo max_iteration del Newton globale. La simulazione produce poi una deformazione plastica equivalente non finita.

Prove diagnostiche senza modificare i materiali/configurazioni del progetto:

- time_step=0.000125, mesh 10x10: Newton non converge dal passo 11 (t=0.001375); proseguire fino al termine produce errori dell'ordine di 1e102, inutilizzabili.
- stesso passo con TANGENT_OPERATOR_ESTIMATION=0 in un materiale temporaneo: Newton continua a non convergere.

La precedente macro ANALYTICAL_TANGENT sovrascriveva sia tensioni sia tangente con formule J2 senza storia, anziche usare la risposta della legge costitutiva. Pertanto le vecchie esecuzioni e quelle attuali non risolvono la stessa risposta materiale oltre lo snervamento. Ripristinare la macro non costituisce una validazione della legge registrata.

Il riferimento analitico e la BODY_FORCE sono basati su una proiezione radiale 3D con sigma_zz. La legge registrata e ridotta a tre componenti: la corrispondenza oltre lo snervamento resta da verificare, come indicato nel README. Le prove non isolano ancora se la divergenza globale dipenda principalmente da questa discrepanza, dalla linearizzazione nonlineare delle condizioni Nitsche o dall'integrazione costitutiva.

Nessuna modifica al materiale plastico, ai passi temporali o alla formulazione plastica e stata applicata durante questa diagnosi. Il problema plastico non e ancora risolto.
