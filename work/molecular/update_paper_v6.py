from pathlib import Path
HERE=Path(__file__).resolve().parent
p=HERE.parent.parent/'outputs/pulizia-grafene-pmma-ppc.tex'
text=p.read_text();assert 'Pressione coerente con il potenziale' not in text
text=text.replace('revisione 5}', 'revisione 6}')
text=text.replace('Revisione 5: orientazione interfaciale e modello molecolare periodico sul Mac',
                  'Revisione 6: pressione corretta, campionamento molecolare e autoassociazione dei solventi')
text=text.replace(r'''HFIP sono completate: rimane una disconnessione osservata fra
i centri 0,75 e 0,85\,nm. La campagna DMF e il percorso inverso
HFIP proseguono.''',r'''HFIP e le undici DMF sono completate: in HFIP rimane una
disconnessione osservata fra i centri 0,75 e 0,85\,nm. Sono
completate anche cinque finestre inverse HFIP, da 1,05 a 0,65\,nm;
le due storie conservano differenze e lacune di campionamento.''')
marker=r'\subsection{Vie indipendenti e lettura critica delle nuove proposte}'
addition=r'''\subsection{Pressione coerente con il potenziale: una correzione necessaria}

Nel percorso installato del modello polarizzabile la deformazione
virtuale della cella modificava i vettori locali, mentre posizioni,
cella reciproca e volume usati dalla parte elettrostatica restavano
non deformati. Questo non altera l'energia a deformazione nulla,
ma omette contributi alla sua derivata. Abbiamo quindi verificato
la pressione in doppia precisione su una configurazione di
16 molecole HFIP, confrontandola con
\begin{equation}
 \sigma_{ii}=\frac{1}{V}\frac{\partial U}{\partial\epsilon_{ii}},
 \qquad P_{\rm pot}=-\frac{\operatorname{tr}\sigma}{3}.
\end{equation}
Lo scarto isotropo \`e circa 183,56\,bar e persiste dimezzando
il passo di deformazione: non \`e un errore della precisione
singola. La prima traiettoria a cella variabile viene quindi
conservata come tentativo rifiutato, senza ricavarne propriet\`a
del liquido.

Abbiamo reso differenziabili anche le geometrie elettrostatiche
deformate e i fattori reciproci gaussiani. Le operazioni a
deformazione nulla mantengono esattamente l'energia e cambiano
le forze di meno di $1{,}5\times10^{-15}$\,eV/\AA\ nel confronto
in doppia precisione. Lo scarto residuo dalla derivata numerica
\`e 0,0106\,bar al passo $5\times10^{-5}$, con riduzione quadratica
al dimezzamento del passo. I pesi appresi sono preservati.

Ripetendo la stessa partenza e lo stesso integratore a cella
variabile, negli iniziali 50\,fs l'escursione massima dell'energia
conservata estesa passa da 0,04393 a 0,001056\,eV: una riduzione
di circa 42 volte. Questo verifica la correzione fisica del
calcolo. La nuova traiettoria di 100\,fs rimane transitoria:
la densit\`a scende dalla condizione iniziale 1,607 a circa
1,236\,g/cm$^3$. Non \`e una densit\`a di equilibrio accettata.

\subsection{Campionare le orientazioni senza attendere la dinamica lenta}

Il calcolo della sola energia del medesimo modello richiede
circa 0,77 secondi per una configurazione da 192 atomi;
il controllo contro il percorso con forze differisce al massimo
di $7{,}7\times10^{-6}$\,eV nei tre punti verificati.
Abbiamo quindi eseguito un campionamento a pressione imposta
che combina traslazioni, rotazioni molecolari, rotazioni del
legame OH e brevi proposte di dinamica hamiltoniana. Tutte le
coordinate interne restano aggiornabili; non si identifica
il liquido con molecole permanentemente rigide.

Per una proposta simmetrica in $\log V$, che scala solo i
centri molecolari conservando le coordinate relative,
l'accettazione \`e
\begin{equation}
 \min\left\{1,\exp\left[-\beta(\Delta U+P\Delta V)
 +(N_{\rm mol}+1)\log\frac{V'}{V}\right]\right\}.
\end{equation}
Il fattore $N_{\rm mol}$ viene dal determinante della
trasformazione dei centri; l'unit\`a aggiuntiva deriva dalla
proposta in volume logaritmico. Per 16 centri ideali, con
$\beta P=1$, il volume deve avere media e varianza 17.
Il controllo indipendente restituisce rispettivamente 16,993
e 17,008. Le proposte hamiltoniane usano quattro passi di
0,25\,fs e vengono accettate attraverso la variazione
dell'energia totale; il numero di mosse non misura un tempo
fisico di pulizia.

La prima serie effettiva di 1.000 mosse comprende
traslazioni, rotazioni, torsioni OH e 19 proposte hamiltoniane.
Le densit\`a campionate scendono fino a circa 1,268\,g/cm$^3$;
non si osservano rifiuti dovuti ai limiti di volume imposti.
L'equilibrio, le partenze indipendenti e l'effetto della
dimensione della cella restano da controllare. Il campionamento
continua preservando configurazione, contatori e stato del
generatore casuale. Non adattiamo il modello alla densit\`a
sperimentale per dichiarare riuscito il confronto.

\subsection{Un controllo chimico ulteriore: autoassociazione dei solventi}

Una simile attrazione locale per il carbonile non implica
una simile competizione nel solvente. Abbiamo aggiunto i
dimeri HFIP--HFIP e acido acetico--acido acetico, mantenendo
lo stesso livello $\omega$B97M-V/def2-TZVPD e la correzione
di base. Sulla griglia non locale iniziale otteniamo interazioni
congelate di $-29{,}744$ e $-88{,}060$\,kJ/mol. L'affinamento
HFIP porta a $-29{,}932$, con variazione di griglia 0,187\,kJ/mol;
MACE-POLAR-M restituisce $-29{,}984$, scarto circa 0,052\,kJ/mol.
Per il dimero dell'acido la griglia affinata d\`a $-88{,}074$,
variazione 0,014\,kJ/mol; il modello d\`a $-88{,}924$. Il dimero dell'acido
ha due legami a idrogeno; il confronto non li tratta come un
singolo legame equivalente. Il dimero non qualifica il liquido intero.

Per evitare di confondere interazione e deformazione,
abbiamo anche costruito, nel medesimo modello, la reazione
elettronica nel vuoto
\begin{equation}
 S_2+\mathrm{acetato}\longrightarrow S\!:\!\mathrm{acetato}+S,
\end{equation}
usando gli stessi monomeri raffinati come riferimento comune.
L'energia risultante \`e circa $-8,57$\,kJ/mol per HFIP e
$+32,75$\,kJ/mol per l'acido acetico. Questo individua una
competizione chimica diversa, assente dal criterio basato sul
solo legame solvente--estere. \`E una reazione fra piccoli
aggregati: entropia, molteplici partner nel liquido e interfaccia
sono ancora assenti. Inoltre la solvatazione del carbonile
pu\`o gi\`a essere presente sia nella catena aderente sia in
quella libera. Il risultato sostiene un'ipotesi da approfondire,
non seleziona da solo un bagno che pulisca il campione.

'''
# Keep the source syntactically standalone.
addition=addition.replace(r'\+non', 'non')
assert text.count(marker)==1;text=text.replace(marker,addition+marker)
p.write_text(text);print('Revision6 source updated from completed computation; scope remains open')
