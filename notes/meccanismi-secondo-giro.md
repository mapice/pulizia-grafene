# Secondo giro: critica meccanicistica e ipotesi discriminanti

6 ottobre 2026. Nota di lavoro basata su fonti primarie e sui nuovi dati forniti da Armando. Non contiene risultati sul suo campione, simulazioni atomistiche eseguite né una ricetta con efficacia dimostrata. Le formule sotto sono bilanci o modelli limite: servono a scegliere misure che distinguano meccanismi.

## 0. Nuovi dati che cambiano la decisione

Sono già stati provati: HCl per 48 h; acido acetico glaciale; acetato di etile; cloroformio; DMF; THF; MPK (sigla mantenuta letteralmente, non interpretata come NMP); acetone; IPA75/H2O25; remover indicato AR600.7; AR300-76 sia a temperatura ambiente sia a80°C. Non completare arbitrariamente le sigle commerciali.

Unico effetto, descritto come scarso: DMF30min → IPA5min → IPA3min → THF30min → IPA5min → IPA3min → N2, con agitatore a500rpm. La somma dei bagni è76min. Non sappiamo se ogni bagno fosse fresco, né geometria e tipo di agitazione.

**Conseguenza:** non raccomandare un'altra volta HCl come prima soluzione. Il suo precedente bibliografico rimane vero, ma il dato del campione lo declassa a tentativo già effettuato. La concentrazione non ancora specificata limita il confronto quantitativo, senza autorizzare a proporre di nuovo48h come se fosse un'idea nuova.

La questione promettente ora è se il modesto risultato della sequenza DMF/THF sia limitato dall'adsorbimento persistente, oppure parzialmente annullato dagli scambi con IPA e dall'asciugatura, oppure se la frazione residua sia un materiale diverso dal PMMA originale.

## 1. HCl: evidenza empirica sì, idrolisi rapida presunta no

### Polimetacrilati: il monomero non predice il polimero

Schönemann, Laschewsky e Rosenhahn, *Exploring the Long-Term Hydrolytic Behavior of Zwitterionic Polymethacrylates and Polymethacrylamides*, Polymers10(2018),639, [DOI10.3390/polym10060639](https://pmc.ncbi.nlm.nih.gov/articles/PMC6403559/). Nei polimetacrilati idrofili studiati gli autori non rilevano idrolisi del carbossilestere dopo un anno in HCl1M a temperatura ambiente, mentre i monomeri reagiscono. Sono polimeri differenti dal PMMA: **non è un limite cinetico numerico sul nostro materiale**. È però una confutazione della scorciatoia «c'è un estere, quindi HCl lo idrolizza utilmente in48h». Schermatura sterica, accessibilità e solvatazione contano.

Pizarro et al., *Synthesis of functional poly(styrene)-block-(methyl methacrylate/methacrylic acid)...*, J.Appl.Polym.Sci.129(2013),2076–2085, [DOI10.1002/app.38923](https://onlinelibrary.wiley.com/doi/10.1002/app.38923), dimostrano idrolisi parziale del blocco PMMA, che resta polimerico. Il testo sperimentale consultato usa polimero dissolto inTHF, HCl acquoso,60°C/48h; la concentrazione dell'acido non è indicata nel passaggio letto. Il prodotto contiene circa20% unitàMAA nel blocco inizialmente PMMA. Non trasferire quel risultato a residui adsorbiti in sola acqua a20–25°C, né copiare una preparazione di polimeri come pulizia del grafene.

Nel PMMA l'idrolisi del gruppo laterale COOCH3→COOH non rompe la catena carboniosa. Può lasciare materiale meno compatibile col solvente organico seguente. Neppure l'eventuale perdita del segnale esterico dimostra rimozione.

### PPC: legami nella catena, ma acidolisi non automaticamente veloce

Jae Hwan Jung, Moonhor Ree e Heesoo Kim, *Acid- and base-catalyzed hydrolyses of aliphatic polycarbonates and polyesters*, Catal.Today115(2006),283–287, [DOI10.1016/j.cattod.2006.02.060](https://www.sciencedirect.com/science/article/abs/pii/S0920586106001465). PPC disciolto inTHF con10% in massa di soluzione acquosa acida/basica,30°C, analizzato mediante viscosimetria/GPC: degradazione molto maggiore in condizioni fortemente basiche; quella acida è inferiore, e in condizioni moderate è molto bassa. Sono stati verificati abstract e passaggio sui metodi, non ricavate costanti dalle figure. È sbagliato dedurre da questa fonte un'emivita del PPC nel bagno proposto. Le basi restano fuori dal quesito.

**Verdetto sul meccanismo HCl:** il precedente diretto Xiao2019 può riflettere distacco, interfaccia, impurità associate, cambiamenti chimici limitati o combinazioni. Non dimostra scissione rapida del PMMA950K o PPC. Il fallimento riferito da Armando è compatibile con questa critica; non lo spiega da solo.

### Perché non inventare una miscela nucleofila

La scissione della catena carboniosa PMMA non è un risultato ottenibile semplicemente sommando «acido forte» e «buon solvente». Alcoholisi/transesterificazione modifica principalmente gruppi laterali; ammine e catalizzatori basici cambiano il vincolo e introducono nuovi residui; fluoruri minaccianoSiO2. Non è stata trovata una prova primaria abbastanza specifica per proporre BBr3 o un'altra miscela reattiva nuova sul campione. THF/HCl della sintesi sopra non è una nuova soluzione già validata. Occorre prima dimostrare trasformazione utile sul polimero di riferimento e compatibilità del sistema, non saltare direttamente al grafene.

## 2. L'adesione multipunto può sembrare insolubilità senza legami covalenti

Huston, Rice e Larson, *Forward Flux Sampling of Polymer Desorption Paths from a Solid Surface into Dilute Solution*, Polymers12(2020),2275, [DOI10.3390/polym12102275](https://pmc.ncbi.nlm.nih.gov/articles/PMC7601496/), simulano modelli di polimero adsorbito e campionamento di eventi rari. Trovano regimi distinti, fra distacco molto lento e diffusione, dipendenti dalla forza di adsorbimento. Non simulano PMMA950K/graphene/DMF: la loro utilità è mostrare perché la solubilità del materiale libero non determina il tempo di desorbimento. Non tutti i segmenti devono staccarsi simultaneamente in qualunque meccanismo, e non è lecito moltiplicare automaticamente energia di un monomero per l'intera catena.

Il lavoro *Structure, Dynamics, and Apparent Glass Transition of Stereoregular PMMA/Graphene Interfaces through Atomistic Simulations*, Macromolecules51(2018),7518–7532, [DOI10.1021/acs.macromol.8b01160](https://pubs.acs.org/doi/abs/10.1021/acs.macromol.8b01160), osserva adsorbimento dei PMMA tramite gruppi estere-metile. Le simulazioni riguardano basse masse molecolari e490–580K: non forniscono una velocità di pulizia a temperatura ambiente.

Per AR-P950K, il rapporto massa caratteristica/unità ripetitiva è dell'ordine di9,5×10^3 unità; non è il numero di contatti col grafene né il grado medio numerico della distribuzione. Il numero di contatti, la distribuzione di loop e tratti adsorbiti, l'età dello strato e la qualità del solvente restano ignoti.

**Conseguenza:** «non esce conTHF/DMF» non prova reticolazione. Un film che si dissolve quando isolato dal substrato ma persiste adsorbito indica un problema interfaciale. Viceversa, una frazione insolubile anche in riferimenti senza grafene giustifica indagini su materiale reticolato, degradato o diverso. La GPC vede soltanto la frazione solubile: non può escludere una rete insolubile senza un bilancio delle frazioni.

## 3. I tre meccanismi prioritari dopo i dati nuovi

### H1 — Rideposizione durante scambio e asciugatura

**Ipotesi specifica:** DMF e/oTHF estraggono parte del materiale; l'IPA successivo cambia la qualità del solvente e può far riapparire una frazione sulla superficie. Il secondo bagnoTHF recupera una quota, ma il passaggio finaleIPA/N2 la deposita di nuovo. È un'ipotesi, non la diagnosi del caso.

Prova minimale a76min, tutti i passaggi in recipienti freschi e con volumi/temperatura/geometria fissati:

| Gruppo | Sequenza, minuti |
| --- | --- |
| B, riferimento reale | DMF30; IPA5; IPA3; THF30; IPA5; IPA3; N2 |
| I, elimina IPA intermedio | DMF30; DMF5; DMF3; THF30; IPA5; IPA3; N2 |
| F, elimina IPA finale | DMF30; IPA5; IPA3; THF30; THF5; THF3; N2 |
| IF, elimina entrambi | DMF30; DMF5; DMF3; THF30; THF5; THF3; N2 |

I5+3min sostitutivi sono due bagni freschi distinti, non una semplice attesa nella soluzione già carica. Campioni assegnati casualmente da ciascun trasferimento; almeno un riferimento di processo non trattato per il medesimo tempo.500rpm non dà da solo la sollecitazione sul campione: tipo di agitatore, dimensioni del recipiente, livello, posizione del campione e movimento del liquido devono essere conservati.

**Che cosa misura:** sequenze complete. I confrontiI−B eIF−F isolano la sostituzione dell'IPA intermedio;F−B eIF−I quella finale, ciascuna nel proprio contesto. Il secondo fattore cambia anche volatilità, bagnamento e forze capillari dell'ultimo liquido. Un miglioramentoF non prova da solo precipitazione chimica; potrebbe essere un effetto dell'asciugatura o della membrana.

**Previsioni discriminanti:** diminuzione dei residui eliminando IPA, accompagnata da polimero rilevato nei bagni e conservazione del grafene, sostiene l'ipotesi. La comparsa delle particelle soltanto dopo l'asciugatura la rafforza. Se superfici misurate nello stesso liquido sono già ugualmente contaminate prima di asciugare, oppure la massa di polimero trascinabile è insufficiente, l'ipotesi di rideposizione finale dominante perde forza. Nessun incremento di cicli senza diagnosi: ripetere la stessa alternanza può ripetere lo stesso deposito.

**Fonti:** Tanaka et al., Langmuir24(2008),296–301, [DOI10.1021/la702132t](https://repository.kulib.kyoto-u.ac.jp/handle/2433/67605), mostrano tramite riflettometria neutronica che anche non-solventi possono gonfiare PMMA e modificarne interfaccia/modulo; non dimostrano precipitatiIPA sul grafene. Konda et al., *Transient Swelling During Development of PMMA Resist*,2024, [DOI10.2494/photopolymer.37.81](https://www.jstage.jst.go.jp/article/photopolymer/37/1/37_81/_article), misurano regimi di dissoluzione/rigonfiamento dipendenti da massa e sviluppatore. Questi risultati motivano il controllo del percorso di scambio, senza provare il meccanismo diH1 nel campione.

### H2 — Distacco interfaciale lento, localizzato ai bordi

**Ipotesi:** la frazione resistente resta adsorbita per contatti multipli o confini di accesso. Gli acidi falliscono perché non cambiano abbastanza l'interfaccia; i buoni solventi sciolgono il materiale libero ma non liberano rapidamente la catena o l'isola aderente.

**Modello limite e previsione:** per isole circolari di altezza costante, attacco dal bordo a velocitàv porta a r(t)=max(r0−vt,0), M(t)/M0=[max(1−vt/r0,0)]² e t_scomparsa=r0/v. Un meccanismo volumetrico semplice dM/dt=−kM ha invece uguale tempo di rimozione frazionale per isole di diversa dimensione. Le formule richiedono isole isolate, nessuna coalescenza, altezza nota e dissoluzione che non ricresca: non sono una legge universale del PMMA.

**Prova:** campi registrati con distribuzione iniziale dei raggi, misure temporali in liquido compatibile o provini sacrificati a tempi diversi, a parità di altezza e trattamento. Cambiamenti al centro senza arretramento del bordo contraddicono il modello di bordo costante; la dipendenza da raggio può anche derivare da diffusione, perciò occorre verificare più di una firma. Raccogliere materiale rimosso e confrontarne distribuzione delle masse, senza presentare la solaGPC come bilancio totale. Se il materiale persiste anche in un riferimento identico privo di grafene, la sola adesione al grafene non spiega il fenomeno.

**Risultato teorico utile:** tempi di immersione più lunghi e un buon solvente non eliminano necessariamente la barriera interfaciale; nessun tempo di completamento può essere calcolato senzav,k o una distribuzione di barriere misurata. L'adsorbimento è una spiegazione alternativa alla reticolazione, non la sua confutazione.

### H3 — Composizione diversa o mista, oppure posizione non accessibile

**Ipotesi:** particelle che sembranoPMMA sonoPPC, oligomeri di altro supporto, adesivo, inclusioniCu/ossido/sali o contaminanti sotto la membrana. Le originiPDMS/nastro sono candidate soltanto se quei materiali sono realmente entrati nel processo; non vanno introdotte retroattivamente come fatti.

**Previsioni:** improntePMMA/PPC distinte su riferimenti degli stessi lotti; frammenti organosilossanici coincidono con un bianco del materiale ausiliario; oppureCu è localizzato nei residui e ritrovato negli eluati. IlSi2p globale sul substratoSiO2 non identificaPDMS. Un piccoC1s carbonilico non distingue i due polimeri. Serve accordo di più frammenti/spettri, con risoluzione spaziale adeguata alle particelle e limiti di rilevazione dichiarati.

**Prove:** riferimentiPMMA, PPC e loro sovrapposizione sottoposti alla vera sequenza76min; prova separata di nastro/supporto effettivamente usato e bianchi dei solventi; ToF-SIMS su gemelli sacrificabili con spettri completi e controlli di matrice; eventuale nano-IR se sensibilità e scala sono sufficienti. Un riferimentoPMMA marcato isotopicamente permette, in un nuovo trasferimento, di misurare il contributoPMMA senza identificarlo dal carbonio totale. L'assenza sotto il limite di un singolo picco non certifica assenza del materiale.

**Fonti:** [Wang et al., Chem.Mater., DOI10.1021/acs.chemmater.6b03875](https://doi.org/10.1021/acs.chemmater.6b03875) usanoPMMA marcato eToF-SIMS; [Lupina et al., ACSNano, DOI10.1021/acsnano.5b01261](https://arxiv.org/abs/1505.00889) mostrano contaminazione metallica anche dopo delaminazione elettrochimica. [Low-Temperature, Dry Transfer-Printing of a Patterned Graphene Monolayer](https://pmc.ncbi.nlm.nih.gov/articles/PMC4673461/) osserva particelle contenentiSi in un trasferimento che impiega effettivamentePDMS, attribuendole a oligomeri silossanici: precedente per una categoria, non identificazione del materiale diArmando.

**Implicazione perHCl:** un acido potrebbe rimuovere ossidi o sali associati senza degradare il polimero; il suo insuccesso non escludeCu metallico o contaminazione sepolta. Non assumere cheHCl senza ossidante dissolva efficacemente ogni forma diCu. Questo resta una diagnosi da misurare, non una ragione per ripetere il bagno.

## 4. Bilancio quantitativo che può falsificare la rideposizione

Sia c_f la concentrazione massica di materiale non volatile nell'ultimo liquido che resta sul campione, V_f il suo volume trascinato e A l'area considerata. Se quel liquido è la sola sorgente di nuovo deposito,

    M_ridep ≤ c_f V_f,
    Γ_ridep = M_ridep/A ≤ c_f h_f,   h_f=V_f/A,
    t_equiv ≤ c_f h_f/ρ.

ρ è la densità del residuo se ne si vuole inferire uno spessore equivalente. La prima disuguaglianza di massa è preferibile perché non richiedeρ. Il limite vale sull'area di bilancio, non sull'altezza locale di una particella: l'evaporazione può concentrare massa in una zona piccola.

Conversione dimensionale verificata:

    t_equiv[nm] ≤ 0.001 c_f[mg/L] h_f[µm] / ρ[g/cm³].

Esempio deliberatamente illustrativo: c_f=1mg/L, h_f=10µm, ρ=1.18g/cm³ produce t_equiv≤0.0085nm. Per spiegare0.5nm medi con quello spessore liquido e quella densità servirebbero almeno59mg/L. Nessuno di questi numeri è stato misurato sul campione. La concentrazione nel bagno raccolto può sottostimare quella locale trascinata: bisogna delimitarne l'incertezza e il volume residuo, non usare il valore medio come limite certo. Le masse adsorbite durante il bagno prima dell'estrazione non sono comprese in questa specifica ipotesi di deposito da film finale.

**Vantaggio:** una spiegazione che sembra plausibile a parole diventa quantitativamente smentibile. Se il polimero massimo trasportabile non basta, bisogna cercare materiale già adsorbito, una sorgente durante il bagno o un'altra contaminazione.

## 5. Bilancio cinetico: plateau non identifica una rete insolubile

Un modello minimo reversibile è

    J = k_d Γ − k_a c (Γ_max − Γ),
    dΓ/dt = −J,
    V dc/dt = A J − Q c.

J è flusso di massa dalla superficie; Γmassa/area; c massa/volume; k_d1/tempo; k_a volume/(massa·tempo); Q volume/tempo. Il modello è omogeneo e non descrive il dettaglio molecolare, il cambio solvente o un'isola specifica. In bagno chiusoQ=0 può generare unplateau con materiale ancora reversibilmente adsorbito. Rinnovare il bagno resetta c e può far ripartire la perdita; un plateau che non risponde al ricambio suggerisce invece desorbimento lento, accesso insufficiente, trasformazione o rete insolubile. Nessuna di queste cause è identificata univocamente dal solo andamento.

Questo schema giustifica misurare gli eluati e distinguere bagni freschi da un'immersione prolungata nella stessa soluzione. Non consente di prevedere un'efficacia in percentuale scegliendo arbitrariamentek_d,k_a.

## 6. Artefatti che possono cancellare un risultato reale

- AFM può spostare materiale o cambiare contrasto con punta contaminata. Usare zone non preanalizzate, registrare parametri/direzione della scansione, controllare punta su riferimento e ripetere con punta nuova se forme identiche si riproducono. Nemes-Incze et al., Carbon46(2008),1435–1442, [DOI10.1016/j.carbon.2008.06.022](https://arxiv.org/abs/0812.0690), mostrano errori di altezza fino a1nm da parametri di oscillazione sul grafene: non scambiare un cambiamento della risposta punta-superficie con film rimosso.
- Per osservazioni in liquido, confrontare i campioni nello stesso mezzo d'immagine: cambiare daIPA aTHF cambia anche le interazioniAFM. Un'immagine in liquido e una a secco non hanno automaticamente altezze geometriche confrontabili.
- ToF-SIMS è distruttiva; i provini misurati non tornano nella prova di recupero come intatti. Dose diXPS/elettroni e Raman va documentata.
- N2 deve avere gli stessi parametri di getto e geometria; fissare intervallo fino alle misure ed esposizione all'aria. Un controllo trasferito e conservato per lo stesso tempo separa l'invecchiamento ambientale dalla chimica.
- Riduzione residui solo sull'intersezione delle aree dove il grafene resta prima/dopo; regioni asportate sono danno. Il99% globale di integrità non esclude che il restante1% contenesse quasi tutto il residuo.
- Un biancoSiO2 verifica contaminazione dei bagni, ma non nucleazione selettiva sul grafene. Dove possibile aggiungere un testimone grafenico pulito, riconoscendo le differenze di superficie.

## 7. Decisione raccomandata al coordinatore

Il secondo giro dovrebbe sostituire la centralità diHCl con una diagnosi causale della sola sequenza che ha mostrato un effetto. Il4bracci76min è una proposta concreta con reagenti già provati: controlla l'ordine degli scambi, non ripete un elenco di solventi. Abbinarlo al bilancio di massa della rideposizione e a riferimenti chimici del residuo lo rende informativo anche se tutte le varianti falliscono. Le tre ipotesi sopra possono coesistere. Non chiamarle nuove scoperte o prova cheArmando abbia torto: la novità del lavoro qui è formulare e discriminare alternative sul caso, senza una rivendicazione di priorità.

# Addendum richiesto: teoria di riorganizzazione ciclica del residuo

## 8. Un ruolo benefico dell'IPA è possibile, ma non automatico

L'ipotesiH1 non è l'unica interpretazione della debole efficaciaDMF→IPA→THF. L'IPA potrebbe anche contrarre un residuo rigonfio e ridurne temporaneamente i contatti con il grafene; ilTHF successivo estrarrebbe la frazione appena resa mobile. Questa è una proposta meccanicistica ulteriore, specifica per il caso e da verificare. Non rivendichiamo che il principio generale dei cicli di solvente sia nuovo.

**Passaggi fisici ipotizzati:**

1. DMF/THF penetra e rende mobile la frazione accessibile del residuo.
2. Il cambio versoIPA modifica l'equilibrio fra contatti polimero-polimero e polimero-grafene. In una finestra favorevole alcune catene si raccolgono in configurazioni con meno punti di adesione.
3. Il ritorno rapido a un buon solvente permette l'estrazione prima che quelle catene ristabiliscano molti contatti.
4. Il liquido fresco allontana il materiale estratto. Il ciclo è aperto, alimentato dallo scambio dei solventi e dalla rimozione del materiale; non è un dispositivo che viola l'equilibrio.

L'IPA può invece produrre l'esito opposto: precipitazione, aumento dell'adsorbimento, o arresto della mobilità. Una catena compatta può restare saldamente appoggiata sulla superficie. «Collasso» non significa «desorbimento». Il segno dell'effetto non si conosce senza il bilancio delle energie interfaciali e della cinetica.

### Condizione energetica necessaria, non ancora misurata

Per confrontare due configurazioni dello stesso volume polimerico in un dato liquidoL, una descrizione grossolana è

    F = γ_PS A_PS + γ_PL A_PL + γ_SL A_SL + F_conf,
    ΔF = (γ_PS−γ_SL)ΔA_PS + γ_PL ΔA_PL + ΔF_conf,

conS grafene, P polimero, A_SL=A_tot−A_PS eF_conf energia libera conformazionale/elastica. Il distacco parziale richiede perdere area di contattoΔA_PS<0; compattazione e penetrazione del liquido cambiano ancheA_PL. Il cambio solvente cambiaγ_PL,γ_SL eF_conf. Solo se la configurazione a meno contatti diventa accessibile e la sua estrazione precede la riadesione può esserci beneficio. Non sono disponibili valori calibrati per il residuoPMMA/PPC reale: la formula non autorizza a prevedere il segno né a trattare il residuo come un film elastico macroscopico.

### Fonte primaria e limiti

- Martins, Plascak e Bachmann, *Adsorption of flexible polymer chains on a surface: Effects of different solvent conditions*, J.Chem.Phys.148(2018),204901, [DOI10.1063/1.5027270](https://arxiv.org/abs/1805.11459): simulazione di polimeri generici mostra che compattazione e adsorbimento dipendono insieme dalla qualità del solvente e dai contatti. Non è una simulazionePMMA/PPC sul grafene.
- Tenopala-Carmona et al., *Real-time observation of conformational switching in single conjugated polymer chains*, Sci.Adv.4(2018),eaao5786, [DOI10.1126/sciadv.aao5786](https://pmc.ncbi.nlm.nih.gov/articles/PMC5817931/): cambi reversibili di conformazione di cateneP3HT ancorate a un'estremità al vetro, al cambio di solvente. Le catene sono ancorate covalentemente e non vengono pulite via: evidenza del cambiamento di conformazione, non del nostro distacco.
- *A Closed-Loop Solvent Recycling Device for Polymer Removal in Graphene Transfer Process*, Separations12(2025),295, [DOI10.3390/separations12110295](https://www.mdpi.com/2297-8739/12/11/295): la pulizia ciclica descritta migliora l'aspetto macroscopico, maXPS resta comparabile all'immersione e mostra polimero residuo. Il dispositivo riguarda ricambio/riciclo solvente e non dimostra i cicli di qualità del solvente ipotizzati qui. È un motivo per non proclamare una soluzione da sole immagini migliori.

## 9. Modello matematico a stati, distinto da una percentuale di successo inventata

DefiniamoL come frazione aderente e poco accessibile, M come frazione temporaneamente mobile ma ancora sul campione, E come estratta nel liquido. Le masse sono normalizzate alla popolazione suscettibile; resta eventualmente una frazionef non rimossa nella finestra considerata.

Per il vettore x=(L,M)^T, ciascun solventej determina

    dx/dt = K_j x,
    K_j = [ −u_j       v_j       ]
          [  u_j  −(v_j+e_j)    ],
    dE/dt = e_j M.

u_j: sblocco; v_j: riadesione; e_j: estrazione. Sono tassi non negativi da misurare, non proprietà note dei bagni. Il solvente povero può avereu_P>0 nella finestra ipotizzata ma anchev_P elevato ede_P piccolo; il buon solvente può averee_G alto. Non si assume che questa sia la realtà diIPA eTHF.

PerN alternanze a esposizioni totaliT_P eT_G:

    x_N = [ exp(K_G T_G/N) exp(K_P T_P/N) ]^N x_0,
    r_N = f + (1−f) (1,1) x_N.

La sequenza raggruppata con le stesse esposizioni è

    x_blocco = exp(K_G T_G) exp(K_P T_P) x_0.

Le matrici non commutano in generale: modificare l'ordine può cambiare l'estrazione perché il sistema conserva memoria conformazionale. Questa è una ragione teorica precisa per provare il ciclo; il vantaggio dipende dai tassi, non discende dal nome del solvente.

**Caso semplificato:** se ogni ciclo termina riportando integralmente la frazione non estratta nello statoL, se i cicli sono identici e il bagno elimina la riadsorzione dal liquido, si ottiene r_N=f+(1−f)(1−η)^N. η è la resa per ciclo da stimare, non un valore assegnato al campione. La formula non è valida automaticamente seM sopravvive, se i siti cambiano con l'età, se i bagni si caricano o se si deposita nuovo materiale.

**Cautela su un falso optimum matematico:** prendendoη=α(τ_P)[e_G/(e_G+v_G)]{1−exp[−(e_G+v_G)τ_G]} eα=a[1−exp(−k_uτ_P)], conτ=T/N, il modello a reset può prevedereη~N^−2 e quindi nessuna rimozione quandoN→∞. Quel limite dipende dal reset arbitrariamente istantaneo, non è una previsione generale. Nel modello continuo con memoria, il limite corretto è

    lim_(N→∞) [exp(K_G T_G/N)exp(K_P T_P/N)]^N
      = exp(K_G T_G + K_P T_P).

Non promettere dunque un numero ottimale di cicli senza stimare i tempi di riorganizzazione/miscela. Un optimum può emergere da cinetiche o ritardi effettivi, ma non è imposto dal principio dell'alternanza.

### Un risultato falsificabile più robusto

Il modello di pura rimozione a un solo stato,

    dm/dt = −k_j m,

predice m_finale=m_iniziale exp(−Σ_j k_j T_j), indipendente dall'ordine dei solventi. Un vantaggio riproducibile dell'alternanza, a uguali tempi totali, ricambi, manipolazioni e stato finale, falsifica questa spiegazione semplice. Non dimostra da solo la nostra sequenzaL→M→E: cambiamenti di composizione locale, trasporto, precipitazione e meccanica sono spiegazioni concorrenti da verificare.

## 10. Esperimento specifico che mette alla prova la teoria

MantenereDMF30min comune,IPA totale16min,THF totale30min, quindi76min totali. Confronto iniziale aN=4:

- **Raggruppato:** DMF30; quattro bagniIPA da4min consecutivi; quattro bagniTHF da7,5min consecutivi; N2.
- **Alternato:** DMF30; quattro ripetizioni diIPA4min→THF7,5min; N2.

Sono gli stessi otto bagni successivi alDMF, con identici volumi, ricambi, manipolazioni, geometria, temperatura e solvente finaleTHF. I bagni consecutivi del gruppo raggruppato devono essere davvero freschi: non sostituirli con un'unica immersione, altrimenti si confondono alternanza e rinnovo. Il protocollo reale diArmando, che terminaIPA, resta un riferimento separato; non chiamare il gruppo raggruppato una replica letterale.

Misure prima/dopo su gemelli; alcuni sacrificati appena dopoIPA e appena dopoTHF. Predizione favorevole alla teoria: dopoIPA diminuisce l'impronta laterale e cresce l'altezza con massa/volume quasi costanti; dopoTHF diminuisce il residuo totale e il polimero compare nell'eluato, senza diminuzione della copertura del grafene. Nessuna depolimerizzazione è richiesta. Se cambia soltanto la forma, la pulizia non è avvenuta. Se eliminareIPA dà prestazioni superiori alla sequenza ciclica, il ruolo beneficoIPA perde sostegno. Se l'alternanza non supera il raggruppamento entro una sensibilità adeguata, questa teoria non ha prodotto un miglioramento per la finestra provata.

La selezioneN=4 è un numero progettuale per distinguere un impulso da un blocco, non un optimum previsto. Qualsiasi simulazione deve esplorare regimi favorevoli, neutrali e sfavorevoli e riportare chiaramente che i tassi sono non calibrati. Simulare soltanto tassi scelti per far vincere il ciclo non aggiungerebbe evidenza.
