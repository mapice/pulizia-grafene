# Terzo giro — una via chimica distinta: solvatazione competitiva con HFIP

Ricerca svolta il 6 ottobre 2026. Responsabilità: candidato chimico non ancora nell'elenco, meccanismo fisico e possibilità di smentirlo. Non sono stati modificati gli elaborati finali né eseguiti esperimenti.

## Decisione

Il candidato che selezionerei per una **qualificazione preliminare sui materiali** è l'**1,1,1,3,3,3-esafluoroisopropanolo (HFIP)**. La ragione è cercare una solvatazione dei carbonili qualitativamente diversa da DMF/THF, combinata con la capacità di dissolvere PMMA libero. Non lo presento come reagente che depolimerizza il PMMA, né come pulizia già dimostrata per AR-P950K/PPC/grafene/SiO2. Non ho trovato una prova primaria diretta di quest'ultimo risultato.

La proposta nuova per questo caso è: **stabilizzare chimicamente la catena distaccata rispetto alla catena aderente usando un solvente fortemente donatore di legami a idrogeno, anziché tentare soltanto un altro accettore polare o l'idrolisi acquosa del PMMA**. HFIP è un candidato materiale per questa ipotesi; la differenza di stabilizzazione deve essere provata. Questa idea non è una rivendicazione di priorità della chimica della solvatazione.

Non propongo concentrazione, durata o temperatura di pulizia del grafene: gli articoli trovati non le giustificano. Il tempo di preparazione di un polimero massivo in un lavoro non va ricopiato come bagno sul grafene. HCl48h, acidoacetico, DMF, THF e gli altri tentativi riferiti restano risultati negativi acquisiti.

## 1. Cosa è davvero documentato per HFIP

### 1.1 PMMA libero in HFIP puro, a temperatura ambiente

Yongxin Wang et al., *Eco-Friendly High-Performance Poly(methyl methacrylate) Film Reinforced with Methylcellulose*, ACS Omega5(2020),24256–24261, [DOI10.1021/acsomega.0c02249](https://pmc.ncbi.nlm.nih.gov/articles/PMC7528171/).

Il lavoro usa PMMA di massa dichiarata99.200 e HFIP99,5%. PMMA e metilcellulosa vengono sciolti a temperatura ambiente, con carico totale indicato circa1,5g/100gHFIP. Il confronto riporta cheHFIP è l'unico solvente fra quelli provati che scioglie **entrambi** i componenti: non dice cheTHF/DMF non sciolganoPMMA da solo. Gli autori colano poi film in uno stampo di vetro e li asciugano anche sotto vuoto a50°C.

**Portata:** prova di dissoluzione diPMMA libero, di massa circaun decimo del950K. Non prova desorbimento, assenza di residuiHFIP, compatibilità elettronica dell'ossido, né solubilità delPPC. Il testo chiama impropriamentePMMA un «polyester»: non ripetere questa classificazione; la sua catena principale è carboniosa.

### 1.2 Indizio a masse maggiori, ma in eluente additivato

Agilent, *Analysis of Polyamides by GPC in HFIP*, nota applicativa con cromatogrammi primari, [PDF5990-7978EN](https://www.agilent.com/cs/library/applications/5990-7978EN.pdf).

Mostra standardPMMA fra1.020 e1.900.000g/mol, compreso790.000, inHFIP con0,02M trifluoroacetato di sodio. Lo standard più grande è escluso dai pori della colonna: non confondere esclusione cromatografica con insolubilità.

**Portata:** supporta l'accessibilità analitica diPMMA ad alta massa in quel mezzo. Non dimostra cheHFIP puro sciolga completamente AR-P950K né che rimuova residui. **Il sale non viene incluso nella proposta di pulizia**: non trasferire additivi della cromatografia al sistemaSiO2/grafene. Il dato elimina solo l'obiezione ingenua cheHFIP sia utilizzabile esclusivamente con oligomeri.

### 1.3 Legame a idrogeno: evidenza molecolare, non rimozione

W. Liu et al., *Tacticity control in the radical polymerization of 2,2,2-trifluoroethyl methacrylate with fluoroalcohol*, J.FluorineChem.123(2003),147–151, [DOI10.1016/S0022-1139(03)00114-3](https://researchportal.hkust.edu.hk/en/publications/tacticity-control-in-the-radical-polymerization-of-222-trifluoroe/).

Gli autori osservano conFTIR interazioni a legameH fra fluoroalcoli e monomeri/specie in crescita; studiano anche copolimerizzazione conMMA inHFIP. Questo supporta l'interazione con il sistema metacrilico, ma non misura ilPMMA adsorbito né una costante utile per il nostro residuo.

**Differenza daHCl/acidoacetico già falliti:** acidità diBrønsted e capacità di solvatare un'interfaccia polimerica non sono la stessa variabile. Un fluoroalcol organico può unire penetrazione/solvatazione e donazione dilegameH. Non concludere però che sia «più acido» dell'acidoacetico o che ogni forte donatore distacchi il polimero.

## 2. Il meccanismo proposto, con un criterio che può fallire

Un estere delPMMA può interagire conHFIP via carbonile; ilPPC contiene gruppi carbonato, ma la sua risposta resta da verificare. L'ipotesi è che la catena distaccata esponga aHFIP più siti favorevoli, o siti meglio accessibili, rispetto alla catena aderente. In tal casoHFIP abbassa selettivamente l'energia libera dello stato distaccato.

### Modello di legame preferenziale, non un altro modello a stati cinetici

Per una catena consideriamoA=adsorbita eD=distaccata. Nel modello di siti indipendenti, n_A,n_D sono numeri efficaci di siti accessibili eK_A,K_D le costanti di associazione conHFIP. A parità delle altre condizioni:

    ΔG_des(a)=ΔG_des(0)
      −kBT [ n_D ln(1+K_D a) − n_A ln(1+K_A a) ],

con a attivitàHFIP eΔG_des=G_D−G_A. Ne segue

    d ln K_des / d ln a
      = n_D K_D a/(1+K_D a) − n_A K_A a/(1+K_A a).

È un modello termodinamico di qualificazione, senza costanti assegnate al campione. La proprietà decisiva è la **differenza** di solvatazione fraD eA. Se tutti i carbonili rimangono ugualmente accessibili da adsorbiti, n_A≈n_D eK_A≈K_D, il vantaggio si cancella anche con un donatore forte. SeHFIP stabilizza ponti con siti ossigenati superficiali, potrebbe persino stabilizzareA.

Questo è il punto che rende la proposta più specifica di «prova un altro solvente»: il bersaglio è una differenza misurabile di associazione e desorbimento. Il modello non garantisce che la stabilizzazione esista sul piano basale del grafene, dove dispersione e contatti estere-metile possono dominare. Non attribuire automaticamente l'adesionePMMA/grafene a una rete di legamiH.

### Perché non partire da una miscela arbitraria

THF e soprattuttoDMF sono anche accettori di legameH e possono legareHFIP, riducendo o modificando la sua attività rispetto alla frazione volumetrica. DiluireHFIP in un solvente già provato può perciò attenuare proprio l'interazione cercata; in altri casi può migliorare la solvatazione del segmento apolare. Senza diagramma di solubilità, attività o spettroscopia, non c'è base per scegliere50:50 o un altro rapporto come ottimo.

La prova preliminare deve quindi stabilire la risposta dei **materiali reali**, separati e insieme, e del residuo effettivo; non supporre cheHFIP sia un solvente universale. Non imporre un risciacquo finaleIPA prima di averne verificato il rischio di precipitazione per ciò cheHFIP ha estratto.

## 3. Qualificazione e criteri di abbandono

Questi sono obiettivi di ricerca, non una ricetta operativa:

1. **Solubilità effettiva:** AR-P672.045 ePPC degli stessi lotti, con la storia termica pertinente, devono risultare effettivamente estraibili nel mezzo candidato. Misurare polimero disciolto o frazione residua; la sola scomparsa visiva non basta. SePPC precipita o resta insolubile, il candidato non è qualificato per la miscela.
2. **Competizione interfaciale:** su un riferimento adsorbito e su un riferimento libero, cercare cambiamenti carbonilici coerenti con associazione reversibile e quantificare materiale esportato. Se il carbonile mostra interazione conHFIP ma il materiale resta aderente, è falsificata l'idea che tale solvatazione sia sufficiente a pulire quell'interfaccia.
3. **Origine della frazione resistente:** l'eluato deve avere l'impronta del residuo rimosso. Se resta una contaminazione di identità diversa, non aumentare arbitrariamente dose otempoHFIP: la proposta ha raggiunto il limite del suo bersaglio.
4. **Compatibilità reale:** su testimoni diSiO2 e grafene senza residuo, verificare ossido, continuità, Raman, drogaggio e contaminazionefluorurata. CheHFIP sia stato colato su vetro o usato per disperdere fogli di grafene non certifica un ossido elettronico o una membrana monostrato trasferita.
5. **Nessuna sostituzione dello sporco:** quantificareF1s/altre firmeHFIP dopo rimozione del solvente. Se diminuisce ilPMMA ma rimane materialeHFIP o cambia persistentemente il grafene, non dichiarare pulizia risolutiva. La volatilità non garantisce l'assenza di specie associate.
6. **Non scissione dichiarata senza prova:** distribuzione di massa e chimica del polimero devono restare coerenti con dissoluzione se questo è il meccanismo rivendicato. Una reazione inattesa colPPC richiede un'analisi separata; non trasformarla retroattivamente in depolimerizzazione selettiva desiderata.

Il primo risultato utile potrebbe quindi essere negativo:HFIP dissolve i materiali liberi ma non la frazione adsorbita. Questo discriminerebbe la nostra ipotesi di solvatazione da un problema di accesso, adsorbimento già stabile inHFIP o legamechimico/rete diversa.

## 4. Perché ho scartato gli altri candidati per questo giro

### TFE

Esiste uso primario diPMMA in2,2,2-trifluoroetanolo per elettrofilatura: [*Incorporation and Deposition of Spin Crossover Materials into and onto Electrospun Nanofibers*, Polymers15(2023),2365, DOI10.3390/polym15102365](https://www.mdpi.com/2073-4360/15/10/2365). La preparazione aggiunge1,35gPMMA a10mLTFE contenenti un complesso e agita12h. È ancora solubilità di materiale libero, non desorbimento. HFIP offre una perturbazione della capacità donatrice più distinta; scegliere entrambi come se fossero due soluzioni provate diluirebbe l'argomento. TFE resta un eventuale confronto meccanicistico, non una raccomandazione parallela.

### TFA e acidi aromatici

Non ho trovato una dimostrazione primaria di rimozione selettiva del residuoAR-P950K/PPC su grafeneSiO2 che giustifichi TFA. La protonazione/solvatazione non è scissione della catenaPMMA; un effetto sulPPC può essere diverso. Non trasferire una miscela usata inGPC o nella sintesi di blocchi polimerici a un campione di grafene senza dati interfaciali.

Piccoli aromatici potrebbero competere per il grafene, ma devono poi essere rimossi più facilmente del materiale sostituito. Un controesempio concreto: A.C.deOliveiraPimenta eJ.E.Kilduff, *Oxidative coupling and the irreversible adsorption of phenol by graphite*, [DOI10.1016/j.jcis.2005.06.075](https://www.sciencedirect.com/science/article/pii/S0021979705007356), misurano adsorbimento e rigenerazione con metanolo e considerano una frazione irreversibile. Non sceglierò fenolo, cresoli o un acido aromatico soltanto perché «aderiscono più forte»: quello stesso fatto può peggiorare la contaminazione.

### Cromatografia dei polimeri: lezione utile, ricetta non trasferibile

Anna M. Caltabiano, Joe P. Foley e André M. Striegel, *Organic solvent modifier and temperature effects in non-aqueous size-exclusion chromatography on reversed-phase columns*, [DOI10.1016/j.chroma.2017.11.027, testo primarioPMC6604611](https://pmc.ncbi.nlm.nih.gov/articles/PMC6604611/), studianoPMMA inTHF con altri modificatori su fasiC18,C4,fenile e ciano. Le interazioni non ideali dipendono dalla superficie e dal modificatore; in certi casi gli alcoli non eliminano l'adsorbimento. Autori eDOI verificati nei metadati editoriali.

Conclusione trasferibile: **buon solvente in soluzione e buon agente di desorbimento da una data superficie non sono sinonimi**. Non trattare silice funzionalizzata, grafene basale e difettiossigenati come un'unica superficie. La cromatografia rende più preciso il criterio di qualificazione diHFIP, ma non ne fornisce l'efficacia sul campione.

## 5. Vera via reattivaC–H: nuovi lavori, ma selettività non dimostrata sul grafene

Ho verificato nei metadati editoriali che i due lavori sotto erano già pubblicati entro il6ottobre2026.

- Ferdinando De Luca Bossa et al., *Depolymerization of Commercial Polymethacrylates Triggered by Hydroperoxides*, JACS2026, pubblicazioneonline3agosto, [DOI10.1021/jacs.6c10099](https://pubs.acs.org/jacsat/article/doi/10.1021/jacs.6c10099/5242090/Depolymerization-of-Commercial-Polymethacrylates). L'abstract primario documenta idroperossidiTBHP/CHP, formazione termica di radicali, estrazione diH eβframmentazione,220°C;70%depolimerizzato ma55%monomero intatto e15%prodotti laterali. È vera rottura della catena, non idrolisi laterale.
- Chih-Yin Lin et al., *Low-Temperature Thermal Depolymerization of Polymethacrylates*, JACS2026, pubblicazioneonline29agosto, [DOI10.1021/jacs.6c13508](https://pubs.acs.org/jacsat/article/doi/10.1021/jacs.6c13508/5361403/Low-Temperature-Thermal-Depolymerization-of). L'abstract primario riporta depolimerizzazione diPMMA non modificato a170°C conperossidi e>95%monomero in condizioni ottimizzate, conradicalialcossile ecloro; l'assenza di solventiclorurati richiede temperature maggiori. Non sono stati verificati tutti i dettagli del supplemento: non trasformare queste informazioni in un protocollo.

Il piano basale ideale delgrafene non ha legamiC–H, ma questo **non basta per la selettività**: radicali eossidanti possono reagire con la reteπ, bordi e difetti; residui polimerici possono anche essere innestati alla superficie. Le condizioni sono inoltre lontane dalvincolo di una semplice pulizia liquida delcampione, e non ci sono prove di conservazione delgrafene con quella chimica.

È una direzione reattiva reale che sarebbe interessante qualificare separatamente, ma non la candidata scelta. La ricerca necessaria sarebbe trovare una finestra in cui la velocità di scissionePMMA prevalga nettamente su funzionalizzazione/etching delgrafene e nuovi depositi. Nessun dato trovato permette di affermare che tale finestra esista sul sistema attuale.

## Conclusione operativa

**Un solo candidato nuovo,HFIP, con un'ipotesi chimica distinta e una soglia di qualificazione esplicita.** Evidenza attuale: solvatazione/dissoluzionePMMA libero e interazioni fluoroalcolo-metacrilico. Parte congetturale: preferenza energetica per catene distaccate che acceleri il desorbimento del residuo reale. Mancanze:AR-P950K puro, PPC, interfaccia ecompatibilità elettronica. Per questo la conclusione è «candidato da qualificare», non «soluzione funzionante» e non «immergere per untempo scelto a caso».

## 6. Aggiornamento cronologico e obiettivo finale: prelievo a secco per lo stack

Le particelle sono già presenti subito dopo il trasferimento suSiO2. «90% èPMMA» è una stima soggettiva dell'esperto, non un risultato composizionale. La membranaPC e la cupolaPDMS appartengono al successivo prelievo a secco: **non possono essere la fonte di particelle già presenti prima**. Non confonderePC del futuroprelievo conPPC del supporto precedentemente descritto.

L'intervento chimico avviene prima che quelPC/PDMS tocchi il campione. Non occorre quindi qualificare l'immersione di una cupolaPC/PDMS già presente; occorre qualificare che il grafene trattato si possa poi prelevare con buona resa e produca interfacce pulite nello stack. Gli esiti pertinenti sono: area/cristalli recuperati, continuità, residui trasportati nello stack, bolle e prestazioni richieste. Una diminuzione controllata di adesione aSiO2 potrebbe migliorare il prelievo, mentre distacco prematuro nel bagno farebbe perdere materiale. La conservazione dell'adesione iniziale non è dunque un requisito assoluto, ma la fase e il momento del distacco contano.

## 7. Valutazione indipendente del binario HFIP/toluene: non promosso

L'idea «HFIP solvata i carbonili, toluene occupa il grafene e impedisce il riattacco» è fisicamente immaginabile, ma **non ha ancora una ragione specifica sufficiente per essere la proposta prioritaria**. Non ho trovato evidenza primaria di spiazzamento diPMMA950K adsorbito sul grafene mediante toluene in una soluzioneHFIP. La sua presenza su una superficie in vuoto non è una misura della competizione in liquido.

### Fonti sulla miscela e sulle associazioni

- D.G.Hutton, brevetto sperimentale [US3284348A](https://patents.google.com/patent/US3284348A/en), descrive misure di equilibrio liquido-vaporeHFIP/toluene; afferma esplicitamente cheHFIP **non forma azeotropo con toluene**, a differenza delbenzene. Documenta l'uso del sistema liquido, non un diagramma di miscibilità per tutte le temperature/composizioni e certamente non il sistema conPMMA/PPC/acqua. Nessuna composizione del brevetto è trasferita alla pulizia.
- Le Lu e Ruimao Hua, *Dual XH–π Interaction of Hexafluoroisopropanol with Arenes*, Molecules26(2021),4558, [DOI10.3390/molecules26154558](https://pmc.ncbi.nlm.nih.gov/articles/PMC8347120/), calcolano interazioniHFIP conbenzene, toluene, anisolo e altri aromatici. Le energie pertinenti sono ottenute **in vuoto**, non costanti di associazione nella nostra miscela. Il dato basta per non chiamare iltoluene un diluente necessariamente inerte rispetto aHFIP; non permette di graduare numericamente il suo spiazzamento sul grafene.
- Federico Caporaletti et al., *Fast Collective Hydrogen-Bond Dynamics in Hexafluoroisopropanol Related to its Chemical Activity*, Angew.Chem.Int.Ed.63(2024),e202416091, [DOI10.1002/anie.202416091](https://pmc.ncbi.nlm.nih.gov/articles/PMC11656151/), mostrano mediante spettroscopia dinamiche di aggregatiHFIP e una netta asimmetria fra capacità donatrice e accettrice. Le attività donatrici di miscele non si ricavano semplicemente dalla frazione volumetrica. Non studiano iltoluene come agente pulente.

### Il limite termodinamico emerso dall'altra revisione

Indicando conP polimero,G grafene,O ossido,L liquido, il costo interfaciale di separazione per unità di area è

    W_PG^L = γ_PL + γ_GL − γ_PG,
    W_GO^L = γ_GL + γ_OL − γ_GO.

Pertanto

    W_PG^L−W_GO^L = γ_PL−γ_OL−γ_PG+γ_GO.

Il termineγ_GL si cancella. A parità delle altre ipotesi, un composto che stabilizzi soprattutto il grafene esposto può facilitare **entrambe** le separazioni e non crea da solo una preferenza per staccarePMMA piuttosto che grafene daSiO2. Accessibilità e barriere di ingresso alle due interfacce possono differire, ma devono essere dimostrate. Per ilprelievo successivo una certa riduzioneW_GO può risultare utile, senza che questo autorizzi la perdita del cristallo nelbagno.

HFIP potrebbe cambiareγ_PL eγ_OL in modo diverso, ma il segno relativo non è noto. Non affermare automaticamente «HFIP lega fortissimoSiOH»:HFIP è un forte donatore e un debole accettore di legameH, e l'ossido presenta siti diversi. La selettività richiede dati di interfaccia, non il conteggio nominale dei gruppiOH.

**Verdetto sulbinario:** non lo promuovo né assegno un rapporto. Aggiunge diluizione, interazioniOH–π, riorganizzazione degli aggregatiHFIP, possibile cambiamento della solubilità e diversa traiettoria di evaporazione. L'assunto vantaggio sulgrafene non è dimostrato. Ilbinario potrebbe diventare un confronto utile soltanto dopo un effettoHFIP specifico e riproducibile; non è un miglioramento già ricavato dalla teoria.

## 8. Un controllo semplice per legareHFIP ai carbonili

Il controllo più diretto resta spettroscopico e non richiede simulare un intero cristallo: seguire carbonilePMMA/PPC eOH diHFIP in una cella liquida chiusa, con riferimenti reali e quantità confrontabili. Prima si cerca un cambiamento reversibile compatibile conassociazione, senza interpretare uno spostamento generico come degradazione.

THF è un accettore che può competere perHFIP e non aggiunge un proprio carbonile. Una prova di aggiunta controllata, con adeguato controllo della diluizione e degli spettri dei soli solventi, può verificare se attenui la firmaHFIP–polimero. Esiste un precedente metodologico: Kunio Oka et al., *Thermochromism and solvatochromism of non-ionic polar polysilanes*, J.Organomet.Chem.611(2000),45–51, [DOI10.1016/S0022-328X(00)00433-2](https://doi.org/10.1016/S0022-328X(00)00433-2), riportano cheTHF spegne una risposta indotta daHFIP inpolisilani contenenti ossigeni eterei. **Non è una prova suPMMA né un fattore di spegnimento già misurato per il nostro sistema.**

Se l'attenuazione della firma di associazione accompagna la perdita di un vantaggio di desorbimento, senza precipitazione né perdita di grafene, il nesso proposto guadagna sostegno. Non sarebbe ancora una dimostrazione esclusiva:THF cambia anche la qualità del solvente e le attività. Se invece si verifica un'interazione carbonilica maHFIP non rimuove più materiale delcontrollo, la semplice «solvatazione dei carbonili» non basta a spiegare o risolvere ilcampione.
