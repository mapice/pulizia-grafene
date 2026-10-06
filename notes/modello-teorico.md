# Modello minimale: distacco, trasporto, riadsorbimento e asciugatura

Nota teorica del 6 ottobre 2026. Non contiene dati misurati sul campione di Armando, né dimostra una procedura risolutiva. Il modello è una proposta falsificabile. In particolare, il riadsorbimento del PMMA presente nel bagno è un'ipotesi da verificare, non la diagnosi del campione. La bassa concentrazione media nel bagno non esclude una partizione favorevole alla superficie, ma nemmeno la dimostra. Non si invoca la saturazione della solubilità del PMMA nel volume del solvente.

**Aggiornamento sperimentale comunicato dall'utente.** I trattamenti già tentati erano agitati a 500 rpm; anche HCl per 48 ore e numerosi solventi non hanno risolto. La sequenza DMF 30 min → IPA 5 min → IPA 3 min → THF 30 min → IPA 5 min → IPA 3 min → N2 ha dato solo un modesto miglioramento. Questi sono resoconti del caso, non misure indipendenti. Non si può attribuire il problema a bagni quiescenti. Il solo numero di giri, senza geometria, volume e viscosità, non determina km o δ. La priorità sperimentale è confrontare la sequenza nota variando separatamente IPA intermedio e IPA finale, e verificare se la massa nell'eluato può spiegare il residuo finale; il modello seguente serve soprattutto a porre vincoli alle spiegazioni.

## 1. Variabili, unità e conservazione

Si consideri una superficie di area A immersa in volume liquido V costante, con ricambio Q (volume/tempo). Il volume esterno è ben mescolato; un coefficiente di trasporto km riassume lo strato liquido vicino alla superficie. I parametri sono effettivi, riferiti a una specie o popolazione di residui e a una composizione/temperatura fissata.

| Simbolo | Significato | Unità SI |
|---|---|---|
| Γ | massa superficiale accessibile e reversibilmente scambiabile | kg m^-2 |
| I | massa superficiale non ancora accessibile al solvente | kg m^-2 |
| Γirr | massa immobile nella finestra sperimentale | kg m^-2 |
| c, cs | concentrazione nel volume e immediatamente alla superficie | kg m^-3 |
| kd | frequenza intrinseca di distacco della popolazione accessibile | s^-1 |
| ka | coefficiente lineare di adsorbimento | m s^-1 |
| km | coefficiente di trasporto attraverso il liquido | m s^-1 |
| ki | frequenza di conversione I → Γ | s^-1 |
| J | flusso netto dalla superficie al liquido | kg m^-2 s^-1 |

Il modello lineare di scambio è

    J = kd Γ − ka cs = km (cs − c).

Eliminando cs:

    J = λ Γ − b c,
    λ = kd /(1 + ka/km),
    b = ka /(1 + ka/km),
    cs = (kd Γ + km c)/(ka + km).

Il riadsorbimento locale può quindi rallentare l'uscita anche quando il volume esterno è pulito. La resistenza del liquido, qui, entra attraverso il ritorno delle catene alla superficie; se ka=0, il modello quasi statico dà J=kdΓ. Per accumulo locale transitorio importante occorre conservare anche la massa nello strato limite, anziché eliminarla algebricamente.

I bilanci sono

    dI/dt = −ki I,
    dΓ/dt = ki I − J,
    V dc/dt = A J + Q(cin − c),
    dΓirr/dt = 0.

Segue esattamente

    d[A(I + Γ + Γirr) + Vc]/dt = Q(cin − c).

Questa identità impedisce di chiamare pulizia la sola redistribuzione. Per cin=0, la massa esportata E(t)=∫ Qc dt completa il bilancio: A(I+Γ+Γirr)+Vc+E è costante. Un bilancio sperimentale richiede anche controlli sull'adsorbimento alle pareti e sul materiale trattenuto nelle tubazioni; in caso contrario la perdita di massa dal campione non coincide necessariamente con l'effluente misurato.

La legge ka cs assume bassa occupazione o un intervallo abbastanza ristretto da permettere una linearizzazione. Non è una descrizione microscopica universale di un polimero multisito. Una possibile estensione è ka cs(1−Γ/Γmax), che conserva i bilanci ma rende la cinetica non lineare. Se coesistono frammenti diversi, ogni specie richiede i propri parametri e il bilancio totale è la somma. Non si deve confondere un coefficiente ottenuto su catene sciolte con quello di un residuo modificato dal trattamento.

## 2. Tre confronti senza dimensioni

Poniamo h=V/A, q=Q/V, K=ka/kd. K ha le dimensioni di una lunghezza: all'equilibrio Γ=Kc. Con Γ0>0, x=Γ/Γ0, y=hc/Γ0, z=I/Γ0 e τ=λt:

    x' = −x + βy + εz,
    y' = x − (β+r)y + r y_in,
    z' = −εz,

dove

    Da = ka/km,
    β = b/(hλ) = AK/V,
    r = q/λ,
    ε = ki/λ.

Da confronta la ricattura alla superficie con il trasporto attraverso il liquido; non misura l'energia di legame del residuo. β confronta la capacità lineare di partizione superficiale con il volume del bagno; non è un rapporto di saturazione del solvente. r confronta il ricambio del volume con il tempo di rilascio effettivo. ε confronta l'accesso alla popolazione schermata con lo scambio della popolazione già accessibile.

Un eventuale rapporto c/csat, necessario per un modello di dissoluzione prossimo alla solubilità, sarebbe un parametro distinto. Qui non è noto e non viene posto vicino a uno.

## 3. Risultati esatti della versione a una popolazione

In questa sezione I0=0, Γirr=0, cin=0. Il bagno parte senza polimero: c0=0.

### Bagno statico

Per Q=0:

    Γ(t)/Γ0 = [β + exp{−λ(1+β)t}]/(1+β).

Il limite è β/(1+β). Un tempo di permanenza molto lungo non supera questo equilibrio del modello. Il plateau esiste anche con concentrazioni molto minori della solubilità del polimero. Se c0 non è zero:

    Γ∞ = K(Γ0 + h c0)/(h+K).

Il bagno contenente il PMMA già disciolto dallo strato di sostegno può quindi aumentare, anziché diminuire, Γ, se la sua concentrazione supera Γ0/K. Questa è una previsione condizionale, non una conclusione sul campione.

### Ricambio continuo

I due tassi adimensionali positivi per r>0 sono

    μ± = [1+β+r ± sqrt((1+β+r)^2−4r)]/2.

La soluzione è

    x(τ) = [(μ+−1) exp(−μ−τ) + (1−μ−) exp(−μ+τ)]/(μ+−μ−).

Se β=0, semplicemente x=exp(−τ), indipendentemente dal ricambio. Per r piccolo rispetto a 1+β, μ−≈r/(1+β). Per r→∞, μ−→1: il tasso resta limitato da λ. Aumentare Q a parità di tutti gli altri parametri diminuisce x e y, per confronto fra sistemi cooperativi, ma il beneficio si esaurisce. Questo risultato non confronta protocolli a uguale consumo di solvente o con geometrie differenti.

### Limite che nessun ricambio supera

Poiché c,I≥0:

    dΓ/dt ≥ −λΓ, dunque Γ(t) ≥ Γ0 exp(−λt).

Inoltre λ≤kd. Il ricambio elimina il ritorno dal bagno; non spezza gratuitamente legami né accelera oltre kd il distacco intrinseco. Se ki è il processo lento, il ricambio non elimina quel ritardo. Con volume esterno idealmente pulito,

    I(t)=I0 exp(−ki t),
    Γ(t)=Γ0 exp(−λt)
          + ki I0 [exp(−ki t)−exp(−λt)]/(λ−ki),

con limite Γ=(Γ0+ki I0 t)exp(−λt) quando ki=λ. Una frazione Γirr resta come limite inferiore; "immobile" significa immobile nelle condizioni e nei tempi specificati, non impossibile da rimuovere in assoluto.

### Bagni successivi: volume e tempo contano entrambi

N bagni puliti di uguale volume V, lasciati ciascuno raggiungere l'equilibrio, danno una frazione [β/(1+β)]^N. Il solvente totale impiegato è però NV.

Se il volume totale è Vtot e ogni bagno usa Vtot/N, definiamo B=AK/Vtot. Nell'ipotesi che ogni stadio raggiunga l'equilibrio:

    xN = [NB/(1+NB)]^N → exp(−1/B).

La versione corretta a tempo totale fissato ttot, con λ costante e durata ttot/N per ogni stadio, è invece

    xN = {[NB + exp{−(1+NB)τtot/N}]/(1+NB)}^N,
    τtot = λ ttot.

Per N→∞:

    x∞ = exp{−[1−exp(−Bτtot)]/B}.

Dimostrazione: la base è 1−[1−exp(−Bτtot−τtot/N)]/(1+NB); moltiplicando il suo logaritmo per N si ottiene il limite indicato. Per B→0 si recupera exp(−τtot); per τtot→∞ si recupera exp(−1/B). Quest'ultimo limite non vale a tempo finito. La costanza di λ e l'istantaneità dei ricambi con un volume per stadio arbitrariamente piccolo sono idealizzazioni: non forniscono un protocollo fisicamente realizzabile con infiniti bagni.

## 4. Un parametro di trasporto effettivamente misurato

Arai, Sawatari, Yoshizaki, Einaga e Yamakawa, *Macromolecules* **29** (1996), 2309–2314, DOI [10.1021/ma951274n](https://pubs.acs.org/doi/10.1021/ma951274n), tabelle 1 e 3, misurano mediante diffusione dinamica della luce PMMA atattico in acetone a 25 °C. Per Mw=9,52×10^5 g/mol riportano D=2,93×10^-7 cm²/s=2,93×10^-11 m²/s; per Mw=4,82×10^5, D=4,29×10^-11 m²/s. Sono campioni a distribuzione stretta e catene già disciolte. Il primo valore può fissare un ordine di grandezza illustrativo per catene di massa simile al PMMA 950K, senza identificarlo con il prodotto commerciale né con i residui dopo lavorazione.

La pagina ufficiale verifica autore, titolo e DOI; i valori completi delle tabelle sono stati leggibili attraverso l'indicizzazione web della pagina editoriale riprodotta da un servizio universitario. Non sono stati ricavati da una curva o da una stima autonoma.

Con D=2,93×10^-11 m²/s, un modello di strato liquido di spessore δ suggerisce km≈D/δ e tempo diffusivo caratteristico tD≈δ²/D. Le seguenti δ sono geometrie ipotetiche dichiarate, non misure del bagno:

| δ | km | δ²/D |
|---|---:|---:|
| 10 µm | 2,93×10^-6 m/s | 3,41 s |
| 100 µm | 2,93×10^-7 m/s | 341 s, circa 5,7 min |
| 1 mm | 2,93×10^-8 m/s | 34.130 s, circa 9,5 h |

Il fattore esatto del tempo diffusivo dipende da geometria e condizione al bordo; δ²/D è una scala, non un tempo di rimozione al 99%. D non determina kd, ka o ki; non misura l'ingresso del solvente in PMMA secco, né la diffusione sotto il grafene. Applicare il valore a questi processi sarebbe ingiustificato. Né questo calcolo né un bagno prolungato provano che il trasporto sia il limite reale: lo strato limite può essere molto minore in presenza di convezione.

## 5. Perché l'asciugatura è un processo separato

Una superficie estratta può trattenere un film liquido di volume Vf, spessore medio equivalente ℓf=Vf/A. Se questo contiene concentrazione cexit di polimero non volatile, la massa trasportata dal film è Vf cexit. Una frazione η, con 0≤η≤1, può depositarsi nell'area osservata:

    Γfinale = Γbagnata + η ℓf cexit.

È un bilancio terminale, non una cinetica universale. Presuppone che non arrivino altre sorgenti e che η riassuma la distribuzione fra superficie, bordo e altrove. Se esiste eterogeneità, ηℓfcexit va interpretato localmente solo dopo aver misurato anche il trasporto laterale. L'evaporazione può cambiare composizione e conformazione, non soltanto aumentare c.

La diminuzione del residuo dopo un ultimo bagno fresco può quindi derivare da minore massa nel film di uscita, anche senza aumentare il distacco in immersione. Per distinguerlo occorre mantenere separati un confronto durante l'immersione e un confronto dopo identica asciugatura. AFM, contatto elettrico o drogaggio da soli non chiudono il bilancio di massa.

## 6. Esperimento discriminante: cambiare c senza cambiare la forza del flusso

Il confronto principale usa la stessa superficie e lo stesso solvente, temperatura e velocità locale del liquido. Un circuito di ricircolo fissa il flusso sul campione e quindi, per quanto possibile, km. Alimentazione e scarico separati del serbatoio cambiano Q senza cambiare la velocità sul campione. V va mantenuto costante e il tempo di transito va misurato con un tracciante compatibile.

1. Si caratterizza il residuo iniziale e si segue una fase a solvente fresco. Si raccolgono frazioni dell'effluente e si includono un circuito privo di campione e un campione senza PMMA come controlli.
2. A parità di idrodinamica si impone un gradino noto della concentrazione di PMMA nel solvente. Servono almeno un'aggiunta del PMMA commerciale nello stesso solvente e, separatamente, il bagno raccolto durante la rimozione; le due soluzioni non sono chimicamente equivalenti a priori. Un tracciante isotopico del PMMA, se disponibile e controllato, distingue nuova deposizione da semplice redistribuzione del residuo originale.
3. Si torna a solvente fresco. Non si asciuga il medesimo campione fra i tre tratti. Repliche sacrificate terminano ogni tratto per analisi chimica della superficie; le tecniche in liquido devono essere controllate per rigonfiamento e perturbazione indotta dalla misura.

Predizione locale, dopo il tempo di transito e prima che Γ cambi apprezzabilmente:

    ΔJ = −b Δc.

Questo cambio di segno o di velocità non è un adattamento arbitrario di una curva di decadimento: è la risposta a una variabile imposta. Se ka>0, aumentando c il rilascio rallenta e, oltre c*=λΓ/b=Γ/K, diventa deposito. La comparsa del tracciante sulla superficie dimostra un percorso liquido → superficie nelle condizioni del test. Non dimostra da sola che esso sia il principale responsabile del residuo originale: va confrontata la massa deposta con quella del plateau.

Una concentrazione sufficientemente bassa o una durata troppo breve possono produrre un esito nullo pur con ka>0. Un esito nullo va riportato come limite superiore alla risposta, specificando sensibilità, Δc e tempi; non come ka=0. La concentrazione deve restare nello stesso regime di solvatazione e viscosità, altrimenti cambiano anche km e gli altri coefficienti.

Misurati b dal gradino e λ dall'andamento iniziale in bagno pulito, si predicono senza riadattarli i plateau al variare di V/A e i tassi al variare di Q/V. Una sola curva monotona non separa desorbimento lento, accesso lento e distribuzione di legami; questa variazione controllata rende il modello confutabile. Le misure di b e λ devono usare un indicatore calibrato di massa/flussso e una popolazione coerente, non il solo valore quadratico medio della rugosità.

Esiti operativi:

- Dipendenza da c e deposito tracciato: il ricambio ha un meccanismo verificato; quantificarne la rilevanza sul residuo reale.
- Dipendenza dalla velocità locale, a c uguale: trasporto nello strato liquido plausibilmente rilevante; stimare km con geometria e D.
- Nessuna dipendenza da c o idrodinamica entro limiti documentati: priorità a distacco, accesso e composizione, senza dichiararli già identificati.
- Cambiamento soltanto dopo asciugatura: controllare inventario del film di uscita e deposito da evaporazione.

## 7. Geometria delle isole: test aggiuntivo, non diagnosi automatica

Un'isola cilindrica di altezza quasi costante h0 e raggio R ha M=ρh0πR². Se il bordo arretra a velocità costante v, R=R0−vt, M/M0=(1−vt/R0)² e il tempo di scomparsa scala come R0. Se la rimozione avviene da tutta la faccia superiore a flusso per area j costante, il tempo è ρh0/j, indipendente da R0 a parità di altezza. Per una legge omogenea di primo ordine, M/M0=exp(−kt), anche il tempo di perdita di una frazione prefissata è indipendente da R0. Un accesso dominato da diffusione laterale fornisce invece la scala R0²/Dint.

Nessuna di queste leggi è universale: altezza correlata al raggio, composizione non omogenea, forma, frammentazione, rigonfiamento e coalescenza possono confondere il confronto. Serve seguire le stesse isole, confrontare gruppi di altezza simile e definire prima una frazione di massa residua o un criterio geometrico osservabile. Il contatto della punta AFM può spostare i residui; finestre di controllo non scandite e repliche sono essenziali per interpretare il test. Un esponente vicino a 1 o 2 non dimostra, da solo, un meccanismo.

## 8. Portata del contributo

Il contributo teorico è un insieme di bilanci, limiti e risposte a perturbazioni controllate. Non identifica il solvente migliore, non fornisce energie di adesione, non prova pulizia atomica e non garantisce che il residuo sia PMMA intatto. Le simulazioni eventualmente costruite da queste equazioni devono chiamarsi esempi del modello e riportare parametri scelti, ipotesi e osservabili. Nessun parametro cinetico di superficie qui è stato misurato sul campione.

## 9. Un limite quantitativo utile per respingere una spiegazione da asciugatura

Se tutto il PMMA presente nel film di uscita ricade sulla medesima area, η=1 fornisce un limite superiore:

    Γrideposizione ≤ cexit ℓf,
    hequivalente ≤ cexit ℓf/ρ.

Per esempio, scegliendo dichiaratamente cexit=1 mg/L, ℓf=10 µm e ρ=1,18 g/cm³, il limite è 0,00847 nm di spessore polimerico equivalente. Questi sono valori di esempio; concentrazione e film liquido vanno misurati, mentre la densità scelta è un riferimento di PMMA denso e non un dato del residuo poroso. Per produrre 1 nm uniforme con quel medesimo film di 10 µm occorrerebbero almeno 118 mg/L, prima ancora di perdite laterali. Una piccola quantità può formare isole alte su una frazione minuscola dell'area: il confronto va fatto con la massa o il volume equivalente per area, non con l'altezza della singola particella.

Se la massa superficiale osservata supera anche il limite ottenuto con l'estremo superiore credibile di cexit e ℓf, la sola evaporazione del liquido di uscita non può esserne la causa principale. Resta possibile riadsorbimento durante tutta l'immersione, alimentato da un volume maggiore: non va indebitamente escluso dal limite sul solo film finale.

## 10. Teoria di una sequenza di solventi: l'ordine può aiutare oppure peggiorare

Questa estensione rappresenta due ipotesi opposte sull'IPA: può mobilizzare una popolazione bloccata, oppure può far tornare bloccata una popolazione che il solvente precedente aveva reso mobile. Sono ipotesi dinamiche; i parametri non sono noti e la loro esistenza non dimostra il meccanismo chimico proposto.

Siano L, M e R masse normalizzate rispettivamente bloccata, mobile ed esportata. R è uno stato assorbente: si assume rimozione effettiva dal sistema durante il contatto con il buon solvente. Il precedente modello con concentrazione del bagno è necessario quando il ritorno dal liquido non è trascurabile. Non si rinomina tale ritorno con il parametro c seguente.

Durante il contatto con IPA (P), transizioni L→M a frequenza a e M→L a frequenza b; durante il buon solvente (G), M→R a frequenza k e M→L a frequenza c:

    P = [[−a, b, 0], [a, −b, 0], [0, 0, 0]],
    G = [[0, c, 0], [0, −(k+c), 0], [0, k, 0]].

Ogni colonna somma a zero e ogni elemento fuori diagonale è non negativo. Ne seguono conservazione L+M+R e positività. Qui c descrive riadesione della quota mobile alla superficie; una quota già nel bagno richiederebbe almeno un quarto stato per distinguere liquido ed esportazione.

Le mappe esatte, senza integrazione numerica, sono semplici. Per un tratto IPA di durata t, s=L+M resta costante e

    Mnuovo = a s/(a+b) + [M−a s/(a+b)] exp[−(a+b)t],
    Lnuovo = s−Mnuovo.

Per un tratto G di durata u, una frazione 1−exp[−(k+c)u] di M lascia quello stato: la frazione k/(k+c) arriva a R e c/(k+c) a L. I casi a+b=0 o k+c=0 sono l'identità.

### 10.1 Predizione di ordine con segno verificabile

Con v=(L,M,R), fare IPA e poi G dà exp(uG)exp(tP)v, mentre l'ordine inverso dà exp(tP)exp(uG)v. La differenza nella massa esportata è esattamente

    ΔR = [k/(k+c)] [1−exp{−(k+c)u}]
          × [1−exp{−(a+b)t}] (aL−bM)/(a+b).

Per piccoli tempi ΔR=k(aL−bM)tu+termini di ordine superiore. È il componente R del commutatore (GP−PG)v. IPA prima del solvente estraente aiuta se aL>bM; danneggia se aL<bM. Dire soltanto che i generatori non commutano non dimostra un vantaggio: il segno dipende dallo stato e dai tassi.

Il confronto matematico termina in solventi diversi. Per un esperimento occorre applicare un passaggio finale identico o misurare la massa esportata prima dell'asciugatura, altrimenti si confonde cinetica con composizione del film finale. Un passaggio comune può attenuare o cancellare il vantaggio d'ordine; va incluso nella stessa composizione di mappe, non assunto irrilevante.

### 10.2 Teorema elementare: alternare più spesso può peggiorare

Poniamo b=c=0, L0=1, M0=R0=0, e fissiamo i tempi totali TI e TG con aTI=kTG=w>0. N alternanze di uguale durata IPA→G danno esattamente

    RN = 1−exp(−w)[1+N(1−exp(−w/N))].

Dimostrazione: a ogni tratto IPA la quota L si moltiplica per p=exp(−w/N). Nel tratto G successivo anche M si moltiplica per p; per induzione LN=p^N e MN=N(1−p)p^N. La conservazione dà la formula. La funzione n[1−exp(−w/n)] cresce strettamente per n>0, perché la sua derivata è 1−exp(−w/n)(1+w/n)>0. Quindi RN diminuisce strettamente con N.

Un solo ciclo rimuove (1−exp(−w))²; infiniti cicli a identici tempi totali rimuovono 1−(1+w)exp(−w). Per w=1 sono 0,399576 e 0,264241. Sono quantità adimensionali di un controesempio teorico, non previsioni sperimentali. Quando la mobilizzazione è irreversibile, conviene in questo modello mobilizzare prima: ogni porzione ha più tempo successivo per essere estratta.

### 10.3 Regime opposto: sottrarre l'intermedio a una reazione inversa

Con b>0 il tratto IPA può ricondurre M a L. L'estrazione fra tratti IPA può allora essere vantaggiosa perché sottrae M a questa reazione inversa. Un esempio adimensionale aTI=1, bTI=9, kTG=100, c=0 dà circa 0,100 rimosso con un ciclo e 0,480 con dieci cicli, mantenendo uguali TI e TG. I valori esatti sono calcolabili con le mappe precedenti. Nessuno di questi numeri è una misura del processo.

Nel limite in cui ogni tratto G estrae praticamente tutto M, ponendo fN=[a/(a+b)]{1−exp[−(a+b)TI/N]}, la quota rimasta dopo N cicli è (1−fN)^N. Questo limite richiede kTG/N grande; non resta valido per N→∞ a kTG fissato. A tempi totali e frequenze fissati, il limite corretto delle mappe è la formula di Trotter:

    limN→∞ [exp(TG G/N) exp(TI P/N)]^N = exp(TG G+TI P).

L'equivalenza non garantisce pulizia migliore. I due regimi precedenti mostrano proprio che non esiste un principio generale "più alternanze = meno residuo".

### 10.4 Esperimento collegato alla sequenza riportata

Il primo confronto tiene uguali DMF, THF, temperatura, agitazione e metodo finale, variando separatamente presenza dei lavaggi IPA intermedi e presenza di quelli finali. I gruppi senza IPA devono avere tempi e maneggiamento controllati. Se serve separare l'effetto della durata totale, si aggiungono controlli con la stessa durata totale nei solventi buoni. Le frazioni di eluato raccolte a ogni transizione permettono di confrontare massa esportata, specie presenti e residuo finale.

Il secondo confronto, solo se il primo dà un effetto riproducibile, fissa i tempi totali di ogni solvente e cambia la segmentazione. Per il ramo mobilizzazione reversibile, l'estrazione subito dopo il tratto IPA deve mostrare una quota esportabile che svanisce o torna bloccata se si attende in IPA; per il ramo precipitazione, introdurre IPA su una quota già mobilizzata deve ridurre l'estrazione successiva. Le predizioni vanno controllate direttamente sull'eluato e sulla superficie, non solo adattando la rugosità finale.

Questo è il contributo teorico più forte disponibile senza dati: un modello conservativo con predizioni di segno e un controesempio analitico contro la superiorità universale dei cicli. La storia sperimentale seleziona le prove da fare, ma non identifica ancora il generatore effettivo di IPA, DMF o THF.
