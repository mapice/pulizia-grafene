# Terzo giro: che cosa deve fare una molecola per staccare PMMA dal grafene

6 ottobre 2026. Nota teorica per selezionare o scartare interventi chimici. Nessuna energia libera di adesione del campione è stata misurata o calcolata atomisticamente in questa sessione. Le condizioni riferite includono numerosi solventi già falliti, agitazione a 500 rpm e scarso miglioramento della sequenza DMF–IPA–THF–IPA. La teoria non deve fingere che il collo di bottiglia sia già identificato.

## 1. Dissolvere PMMA libero e staccare PMMA aderente sono due problemi termodinamici distinti

Denotiamo polimero P, grafene G e liquido L. Il lavoro reversibile per unità di area necessario a sostituire un contatto P–G con due interfacce immerse è

    W_PG^L = γ_PL + γ_GL − γ_PG.

Le γ sono energie libere interfaciali a composizione e temperatura specificate, non tensioni superficiali del solvente puro intercambiabili fra loro. Se W>0, separare il contatto costa energia libera; se W<0, l'inserimento di liquido sarebbe favorevole per quell'area e per quegli stati, purché il liquido possa arrivarci. Un W positivo piccolo può comunque essere superato da energia elastica, lavoro esterno o guadagno di solvatazione delle catene. W negativo non implica velocità elevata: l'accesso può avere una barriera.

Per un'isola reale la variazione complessiva deve comprendere termini ulteriori:

    ΔG = A_separata W_PG^L + ΔG_conf + ΔG_elastica
         + ΔG_miscelazione + ΔG_bordo.

Non si deve contare due volte la solvatazione: nella formula a superfici macroscopiche essa entra già nelle γ delle fasi appropriate; aggiungere ΔG_miscelazione ha senso solo quando si definisce esplicitamente un ulteriore passaggio da polimero compatto a catene disperse. Per catene isolate adsorbite, una differenza di energia libera fra stato adsorbito e stato solvato è spesso meglio definita delle tre γ di un film nanometrico.

Un buon solvente per il polimero libero abbassa il costo dello stato disciolto. Non garantisce un basso costo di sostituzione del contatto, accesso al contatto, bassa barriera di distacco di una catena agganciata in più punti, o rimozione di legami covalenti. Questa distinzione spiega perché molti buoni solventi possano fallire senza postulare automaticamente reticolazione o saturazione del bagno.

### Una relazione utile, ma con ipotesi strette

Se gli stati delle superfici restano gli stessi, non vi sono rigonfiamento, dissoluzione, reazione o isteresi rilevante, le equazioni di Young permettono di scrivere

    W_PG^L = W_PG^v − γ_LV(cos θ_P + cos θ_G),

dove v indica il vapore, θ_P e θ_G sono angoli del medesimo liquido sui due solidi. Ciò mostra anche perché "tensione superficiale minore" non sia, da sola, il criterio giusto: entra un prodotto con i coseni di due angoli. Per PMMA in un buon solvente le ipotesi sono spesso violate; un angolo dinamico su un film che si sta dissolvendo non determina W. Il grafene va inoltre misurato sullo stesso supporto e nello stesso stato di contaminazione: non basta un valore universale di "energia superficiale del grafene".

## 2. Un dato primario quantitativo che NON va trasferito arbitrariamente al bagno

Pastore Carbone et al., *Wettability of graphene by molten polymers*, **Polymer 180**, 121708 (2019), [DOI 10.1016/j.polymer.2019.121708](https://doi.org/10.1016/j.polymer.2019.121708), [manoscritto accettato](https://arxiv.org/pdf/2108.10007), misurano PMMA circa 120 kg/mol fuso su grafene CVD fra 170 e 200 °C, sotto vuoto. Dal contatto di circa 31° ricavano un lavoro di adesione decrescente da circa 0,055 a 0,051 J/m². Non è W_PG^L in DMF, THF o miscele a temperatura ambiente. La stessa fonte precisa il ruolo di interazioni secondarie, assenza di incastro e tensioni residue. Il valore fornisce un riferimento interfaciale concreto, non un parametro utilizzabile direttamente per il residuo attuale.

In questa ricerca non è stata reperita una misura trasferibile di γ_PL, γ_GL e γ_PG per AR-P 672.045/grafene/SiO2 nel liquido effettivo e nella sua storia termica. Sarebbe scorretto produrre da questi dati una classifica numerica di solventi per energia di distacco.

## 3. Selezione molecolare: due proprietà da migliorare insieme

Per una miscela costituita da solvente portante S e componente D capace di competere all'interfaccia, occorrono due verifiche separate:

1. Il polimero liberato deve rimanere solvato lungo tutto il percorso di composizione; la miscela non deve precipitarlo nella fase di estrazione.
2. D deve diminuire il costo o la barriera del distacco in liquido, oppure ridurre il riaggancio dei contatti che si aprono spontaneamente. Essere un buon solvente del PMMA libero non prova questa seconda proprietà.

Un possibile criterio operativo è trovare una composizione φ per cui, rispetto al portante S,

    ΔΔG_ads(φ) = ΔG_ads(φ) − ΔG_ads(S) > 0,

con ΔG_ads definita negativa per adsorbimento favorevole, mantenendo una fase liquida che solva il polimero. Un criterio cinetico diverso è ΔΔG_distacco^‡(φ)<0 oppure un rapporto di riaggancio/distacco più piccolo. Le due cose non sono equivalenti: la stabilità di equilibrio non fissa da sola una barriera. Queste quantità possono essere misurate indirettamente con un confronto di desorbimento e riadsorbimento, oppure calcolate con un modello atomistico validato.

La composizione alla superficie può essere diversa da quella nel volume. L'attività del componente D all'interfaccia, non soltanto la percentuale introdotta, controlla la competizione. Una forte affinità di D per il grafene è però anche un motivo per cui D può diventare il nuovo contaminante. La verifica deve includere la sua successiva esportazione e una misura chimica sensibile a D; uno spostamento di drogaggio non basta.

### Hansen: filtro di solubilità, non energia di distacco

La distanza convenzionale

    R_A² = 4(δ_D,P−δ_D,L)² + (δ_P,P−δ_P,L)² + (δ_H,P−δ_H,L)²

e una sfera empirica R_A<R0 possono selezionare composizioni da provare per la solvatazione. Non danno W_PG^L, l'energia di un legame al grafene o una barriera di accesso. Le medie volumetriche δ_i(φ)≈(1−φ)δ_i,S+φδ_i,D sono approssimazioni; associazione specifica e composizione preferenziale possono renderle inaccurate. Nel modello lineare, se entrambi gli estremi sono dentro la stessa sfera convessa, il segmento lo è; non segue che tutte le miscele reali siano buone per il residuo modificato.

### Un controllo sperimentale recente vicino alla massa 950K

Nakamoto, Fujiie e Yamamoto, *Industrial & Engineering Chemistry Research* **64** (2025), 81–86, [DOI 10.1021/acs.iecr.4c03728](https://pubs.acs.org/doi/10.1021/acs.iecr.4c03728), Tabella S1: PMMA Mw=1.155.000, 50 mg in 0,5 mL, agitazione e lettura visiva a 24 h. THF, acetone, DMF, NMP, MEK e DMAc sono classificati solubili; DMSO, GBL e ciclopentanone insolubili nelle condizioni della prova. È una prova concentrata di dissoluzione, non pulizia del grafene; non stabilisce insolubilità universale a ogni temperatura e concentrazione. Aiuta tuttavia a scartare proposte fondate soltanto sull'etichetta "solvente forte". La fonte studia anche la dipendenza dei parametri dalla massa molecolare; non ho accesso ai valori completi del testo a pagamento e non ne invento una terna per 950K.

La classificazione di alcune sostanze, fra cui cloroformio, dipende molto dalle condizioni della prova: non usarla per negare i trattamenti effettivamente effettuati o tutta la letteratura su altri PMMA. La conclusione utile è che la prova di solubilità sul lotto e alle concentrazioni pertinenti rimane necessaria.

**Conseguenza per la scelta:** DMAc può essere un controllo aggiuntivo di solvente portante, avendo evidenza di dissoluzione di PMMA ad alta massa e non comparendo nell'elenco comunicato; non c'è qui evidenza che stacchi il contatto meglio del DMF già fallito. NMP ha evidenza di solvatazione ma nessuna vittoria deducibile sulla superficie. Non promuoverei DMSO, GBL o ciclopentanone come rimedio semplicemente più forte. Una proposta strategicamente diversa deve dimostrare competizione interfaciale, non soltanto cambiare il nome del buon solvente.

## 4. Perché molte interazioni deboli producono un residuo persistente

Una catena di PMMA nominale 950 kg/mol contiene circa 9.500 unità MMA, ma solo una frazione ignota è a contatto. Non è giustificato moltiplicare l'energia di un monomero per tutte le unità e chiamare il risultato barriera di distacco: le catene possono aprire i contatti in successione, formare anse e staccarsi dal bordo.

Un modello esplicito, utile perché mostra quale proprietà cambiare, considera m ancoraggi reversibili. Se j sono chiusi, si aprono al tasso totale j k_off; i rimanenti si richiudono al tasso (m−j)k_reb. Lo stato j=0 è assorbente solo nell'ipotesi che la catena venga esportata prima di tornare. Il rapporto

    q = k_reb/k_off

misura il riaggancio, non la solubilità massiva. Siano T_j i tempi medi di prima uscita verso zero. Le equazioni all'indietro sono

    −1 = (m−j)k_reb(T_{j+1}−T_j)
          + j k_off(T_{j−1}−T_j),   T_0=0.

Con Δ_j=T_j−T_{j−1}, si ottiene la ricorrenza esatta

    Δ_m = 1/(m k_off),
    Δ_j = [1+(m−j)k_reb Δ_{j+1}]/(j k_off),
    T_m = Σ_j Δ_j.

In particolare

    T_1 = [(1+q)^m−1]/(m q k_off),

con limite T_1=1/k_off per q=0. Se il riaggancio è assente, T_m=H_m/k_off; se q è grande, per m fissato il termine dominante è q^(m−1)/(m k_off). Dunque un componente che impedisca la richiusura può ridurre drasticamente il tempo senza spezzare la catena né postulare un grande aumento di k_off.

Questo modello non è una stima dei tempi reali: m, k_off e k_reb non sono misurati, i contatti polimerici sono correlati e cambiano geometria. È una derivazione condizionale che identifica il bersaglio: impedire il riaggancio dei contatti aperti, invece di inseguire soltanto la dissoluzione del polimero libero.

### Competitore che occupa il sito appena liberato

Se D occupa rapidamente un sito libero con costante di adsorbimento K_D e attività a_D, un modello locale a siti indipendenti dà frazione libera f0=1/(1+K_D a_D). Se l'apertura del contatto polimerico resta invariata e la chiusura richiede sito libero,

    q_eff = q0/(1+K_D a_D).

Inserirlo nella formula precedente rende esplicita una previsione: a parità di solvatazione, viscosità, temperatura e trasporto, l'accelerazione del rilascio dovrebbe correlare con la copertura interfaciale del competitore. D non deve necessariamente espellere istantaneamente una catena: può intercettare i contatti che si aprono. L'effetto si satura e D va poi rimosso. L'equazione è una riduzione modellistica da verificare, non una legge già stabilita per DMF/THF/IPA.

## 5. Perché una classifica di energie di adsorbimento nel vuoto non risolve la selezione

Patil e Caffrey, *Adsorption of common solvent molecules on graphene and MoS2 from first-principles*, [arXiv:1806.07843](https://arxiv.org/pdf/1806.07843), calcolano molecole isolate di toluene, benzaldeide, ciclopentanone, IPA, NMP e cloroformio. Le energie elettroniche di legame sono nell'intervallo circa −0,4/−0,79 eV, usando DFT con dispersione. Non sono energie libere di adsorbimento dal liquido. Il costo di togliere una molecola dalla sua solvatazione, la competizione con le altre molecole, entropia, copertura e volume escluso sono assenti da quella quantità. Non si può quindi convertire il valore più negativo in "miglior solvente per pulire PMMA".

Il dato è utile per controllare la componente dispersiva di un modello atomistico. Anche l'adesione del PMMA fuso può essere un controllo separato. Nessuno dei due controlli, da solo, valida la previsione di scambio PMMA/solvente su grafene immerso.

## 6. Un intervento falsificabile: gradiente che conserva la solvatazione e modifica l'adesione

Una miscela S/D è promettente se la sua aggiunta al bagno reale non precipita PMMA/PPC e se diminuisce l'adesione o la ricattura sul grafene. La prova preliminare deve usare almeno il lotto PMMA libero, il residuo dopo la storia termica e il materiale già estratto; i tre non sono automaticamente identici.

L'intervento proposto sul campione consiste nel mantenere un portante già compatibile con il polimero, variare gradualmente l'attività del componente interfaciale evitando la regione di precipitazione, esportare il materiale liberato e tornare al portante senza passaggio attraverso un liquido che lo faccia riprecipitare. Non c'è un rapporto numerico ottimale deducibile dai dati raccolti. Un gradiente è giustificato quando evita una regione di precipitazione/transitorio identificata: non è vantaggioso per definizione.

**Esperimento che distingue il meccanismo:** confrontare S, S/D e un controllo a viscosità e qualità solvente comparabili, misurando contemporaneamente rilascio del PMMA e occupazione da D di un testimone grafenico. Se D sostituisce il contatto, il suo arrivo all'interfaccia deve precedere o accompagnare l'esportazione del polimero; il trattamento di solo polimero libero non deve spiegare tutto il beneficio. Dopo il ritorno a S, misurare anche la rimozione di D. Uno spostamento del Dirac point o una minore rugosità da soli non verificano il meccanismo.

Misurare la forza di separazione di un colloide rivestito con il medesimo PMMA contro grafene supportato, nei due liquidi e a più velocità, può confrontare direttamente adesione e isteresi, ma occorre controllare il rigonfiamento/dissoluzione del rivestimento. In condizioni non stazionarie si parla di forza apparente dipendente dalla storia, non di lavoro termodinamico W. Un'alternativa è un esperimento spettroscopico o gravimetrico di desorbimento/riadsorbimento su interfaccia riprodotta, con correzione per solvente trattenuto.

### Un candidato da non confondere con pulizia: formammide

Suk et al., *Enhancement of the Electrical Properties of Graphene Grown by Chemical Vapor Deposition via Controlling the Effects of Polymer Residue*, [DOI 10.1021/nl304420b](https://doi.org/10.1021/nl304420b), [PDF degli autori](https://utw10193.utweb.utexas.edu/Archive/RuoffsPDFs/340.pdf), riportano miglioramenti elettrici dopo formammide ma morfologia AFM simile e presenza di formammide rilevata mediante XPS. Interpretano il risultato come compensazione del drogaggio da molecole solvatate nel residuo polimerico. Questa fonte non dimostra che formammide rimuova il PMMA; è un controllo concreto contro la promozione di un solvente sulla sola mobilità.

## 7. Il calcolo atomistico utile e quello che sarebbe una falsa promessa

Una simulazione realistica non è il confronto fra una molecola MMA nel vuoto e una superficie perfetta. Il problema minimo informativo comprende:

- oligomeri PMMA con più lunghezze e stereochimica controllata, per esempio 10/20/50 unità come modelli locali dichiarati;
- grafene sul supporto pertinente o una giustificazione dell'assenza del supporto, includendo casi separati di piano basale e difetto rappresentativo;
- solvente esplicito, composizioni effettive e acqua residua se presente;
- stati iniziali sia adsorbiti sia parzialmente aperti, repliche indipendenti e verifica che il risultato non dipenda dalla configurazione iniziale;
- una superficie di energia libera in almeno numero di contatti e grado di apertura/accesso del solvente, anziché un'unica distanza che nasconda barriere conformazionali.

La differenza ΔΔG_ads fra due liquidi e la barriera per sostituire un piccolo tratto adsorbito sono obiettivi più realistici della promessa di un tempo di pulizia per PMMA950K. La catena 950K comprende circa 9.500 unità: non è rappresentata da un trimero, e le sue lentezze collettive non emergono da pochi nanosecondi di dinamica molecolare.

Un esempio primario di scala modellata è *Gradient of Segmental Dynamics in Stereoregular Poly(methyl methacrylate) Melts Confined between Pristine or Oxidized Graphene Sheets*, [testo aperto](https://pmc.ncbi.nlm.nih.gov/articles/PMC7962820/): catene di 20 unità e geometria confinata. È un precedente per descrivere contatti e dinamica locale, non una validazione di cinetiche di pulizia in solvente del prodotto commerciale.

Prima di fidarsi di una classifica occorre controllare proprietà del liquido e del polimero, solvatazione degli oligomeri, adsorbimento del solvente sul grafene e interazioni incrociate. Combinare parametri di forza generici con la regola di Lorentz–Berthelot non costituisce da solo una validazione interfaciale. Una parametrizzazione non reattiva non può verificare scissione, reticolazione o idrolisi. Se l'effetto previsto cambia segno al variare ragionevole della parametrizzazione o della lunghezza dell'oligomero, il calcolo non seleziona ancora il reagente.

I metodi con campionamento vincolato possono fornire energie libere solo dopo controllo di sovrapposizione delle finestre, convergenza, campionamento dei contatti e incertezza. Non trasformare l'assenza di distacco spontaneo in una breve traiettoria in "residuo insolubile". In questa sessione non sono stati eseguiti tali calcoli né creati numeri da presentare come loro risultato.

## 8. Verdetto della teoria fisica

La teoria seleziona il tipo di intervento: **mantenere il polimero solvato mentre si diminuisce l'adesione in liquido o si intercettano i contatti che si riaprono, poi esportare sia polimero sia competitore**. Non seleziona ancora una molecola vincente, perché mancano energie libere o misure cinetiche interfaciali per i liquidi candidati. Questo limite è informativo: i dati di solubilità già escludono alcune scorciatoie e la prova di competizione interfaciale può distinguere un avanzamento da un altro bagno in un buon solvente.

Per materiale covalentemente legato, reticolato, carbonizzato o sepolto, questa classe di interventi può fallire. Le misure di identità e posizione già proposte restano necessarie per sapere se il bersaglio fisico descritto esiste davvero nel campione.
