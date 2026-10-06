from pathlib import Path
import json
HERE=Path(__file__).resolve().parent
paper=HERE.parent.parent/'outputs/pulizia-grafene-pmma-ppc.tex'
text=paper.read_text()
text=text.replace('Revisione 4: traiettorie interfaciali e controlli dei modelli polarizzabili','Revisione 5: orientazione interfaciale e modello molecolare periodico sul Mac')
text=text.replace('revisione 4}', 'revisione 5}')
assert 'Una variabile nascosta osservata' not in text
summary=json.loads((HERE/'results/neural-liquid-numerics/summary.json').read_text())
assert len(summary['nve_branches'])==2
old=r'''\`E avviata una campagna pi\`u fitta di undici finestre,
500\,ps ciascuna, con $K=750$ e centri da 0,45 a 1,45\,nm,
per HFIP e DMF. Le finestre completate migliorano la sovrapposizione
nella regione aderente; il percorso inverso, altre conformazioni e
la convergenza delle variabili nascoste sono ancora necessari.'''
new=r'''La campagna pi\`u fitta usa undici finestre di 500\,ps,
con $K=750$ e centri da 0,45 a 1,45\,nm. Le undici finestre
HFIP sono completate: rimane una disconnessione osservata fra
i centri 0,75 e 0,85\,nm. La campagna DMF e il percorso inverso
HFIP proseguono. Un controllo delle conformazioni rivela inoltre
il limite del solo centro di massa, descritto qui sotto.'''
assert text.count(old)==1;text=text.replace(old,new)
marker=r'\subsection{Un generatore polarizzabile di parametri e il limite della verifica nel vuoto}'
addition=r'''\subsection{Una variabile nascosta osservata: l'orientazione della catena}

Abbiamo confrontato due campioni HFIP di 500\,ps con lo stesso
vincolo $z_0=1{,}05$\,nm e $K=750$, avviati dalla catena aderente
e dalla catena libera. Negli ultimi 250\,ps la distribuzione della
distanza ha sovrapposizione empirica circa 0,355; la distribuzione
congiunta di distanza ed estensione normale non ha sovrapposizione
nelle celle usate, di 0,01 e 0,005\,nm. La catena conserva
orientazioni diverse. Questo \`e un risultato del modello e
delle traiettorie effettive, non una barriera del residuo reale.

Definiamo l'estensione normale con le masse $m_i$:
\begin{equation}
 q=R_{g,z}=\left[\frac{\sum_i m_i(z_i-z_{\rm CM})^2}
 {\sum_i m_i}\right]^{1/2}.
\end{equation}
Per isolare l'effetto a distanze effettivamente uguali consideriamo
le configurazioni con $1{,}00\le z<1{,}04$\,nm.
Le rispettive altezze medie coincidono entro $8\times10^{-5}$\,nm,
mentre cambiano estensione e prossimit\`a al foglio
(Tabella~\ref{tab:orientation}). I punti sono correlati;
le medie non sono stime su repliche indipendenti.

\begin{table}[ht]
\centering\small
\caption{Confronto a pari intervallo di altezza, nel modello
PMMA4/grafene rigido/HFIP. L'energia \`e l'interazione diretta
Lennard--Jones, non un'energia libera di adsorbimento.}
\label{tab:orientation}
\begin{tabular}{lrr}
\toprule
 & Dallo stato aderente & Dallo stato libero\\
\midrule
Configurazioni salvate &52&63\\
$\langle z\rangle$ (nm)&1,02075&1,02067\\
$\langle R_{g,z}\rangle$ (nm)&0,27694&0,15384\\
Distanza minima atomi pesanti--G (nm)&0,51385&0,73865\\
Interazione diretta P--G (kJ/mol)&$-5,97$&$-0,85$\\
Legami H HFIP--carbonile&3,25&3,54\\
\bottomrule
\end{tabular}
\end{table}

La catena proveniente dallo stato aderente resta quasi verticale;
un'estremit\`a \`e pi\`u vicina al foglio. La seconda resta
quasi parallela. La solvatazione del carbonile non separa da sola
i due stati. Il risultato fornisce una variabile concreta per il
campionamento, anzich\'e attribuire il mancato distacco a una
generica scarsa solubilit\`a.

Se $F(z,q)$ comprende la misura statistica delle due coordinate,
il profilo lungo $z$ soddisfa, a meno di una costante,
\begin{equation}
 F(z)=-RT\log\!\int e^{-F(z,q)/(RT)}\,\mathrm dq+C.
\end{equation}
Una traiettoria che mantiene una sola popolazione in $q$ non
esegue questa integrazione, anche se le distribuzioni di $z$
sembrano collegabili. Abbiamo implementato $R_{g,z}$ nel motore
di simulazione e verificato i valori su 1.002 configurazioni:
errore massimo inferiore a $4\times10^{-10}$\,nm rispetto al
calcolo indipendente. Le forze del vincolo concordano con la
formula analitica entro $5\times10^{-7}$ e con le derivate
numeriche entro $3{,}7\times10^{-4}$\,kJ\,mol$^{-1}$\,nm$^{-1}$.
Due nuovi campioni di 100\,ps con vincoli simultanei su $z$ e
$R_{g,z}$ rimangono distinti. Il controllo numerico \`e riuscito;
il campionamento resta insufficiente e il profilo non viene accettato.

'''
assert text.count(marker)==1;text=text.replace(marker,addition+marker)
marker=r'\subsection{Vie indipendenti e lettura critica delle nuove proposte}'
addition=r'''\subsection{Esecuzione del modello polarizzabile nel liquido sul Mac}

Abbiamo reso eseguibile MACE-POLAR-1-M su una cella periodica
di 64 molecole HFIP, 768 atomi, sulla GPU del Mac disponibile.
I pesi del modello sono preservati. Il calcolo delle reti radiali,
dei prodotti e delle somme sui vicini procede a blocchi;
le derivate vengono ricalcolate dai rispettivi ingressi.
Le operazioni sono equivalenti per le prime derivate, con
parametri congelati. Non viene usato questo percorso per Hessiane
o addestramento.

Una contrazione onerosa pu\`o inoltre essere riordinata senza
costruire tutte le coppie di canali:
\begin{equation}
 \sum_{v,d}w_{uv}x_{bud}y_{bvd}
 =\sum_d x_{bud}\left(\sum_v w_{uv}y_{bvd}\right).
\end{equation}
La sottrazione dei riferimenti atomici additivi \`e una scelta
dello zero di energia che preserva forze e differenze a composizione
fissa. Una cella periodica diluita di 96 atomi controlla il percorso
contro il calcolo originale in doppia precisione: le modifiche
al riordino delle operazioni preservano le forze entro
$2\times10^{-15}$\,eV/\AA. Il confronto con la GPU in precisione singola
ha scarto relativo quadratico medio $1{,}5\times10^{-5}$.
Sulla cella completa il confronto processore--GPU, entrambi in
precisione singola, d\`a $1{,}1\times10^{-6}$; lo scarto energetico
\`e $1{,}3\times10^{-4}$\,eV. Il controllo completo in doppia
precisione sul processore supera il limite di memoria imposto
e non viene considerato concluso.

La GPU \`e limitata al 25\% della memoria assegnabile al processo;
il calcolo effettivo delle forze richiede circa 9,9 secondi
per configurazione. Le derivate numeriche di energia rispetto
a spostamento atomico, traslazione di una molecola e deformazione
della cella concordano con le derivate calcolate entro lo 0,34\%,
1,24\% e 0,84\% ai passi documentati. Due traiettorie a energia
totale costante, con identica partenza e durata di 8\,fs, hanno
escursioni energetiche massime di 0,03929 e 0,00979\,eV per
passi di 0,5 e 0,25\,fs. La riduzione \`e coerente con l'errore
quadratico del metodo di integrazione.

Queste verifiche rendono concretamente disponibile un modello
pi\`u avanzato per ulteriori controlli. La configurazione iniziale
proviene dal modello classico; otto femtosecondi non misurano
densit\`a o permittivit\`a di equilibrio. La qualificazione del
liquido e dell'interfaccia rimane un passaggio distinto.

'''
assert text.count(marker)==1;text=text.replace(marker,addition+marker)
paper.write_text(text)
print('Revision5 source updated from completed raw data; research remains active')
