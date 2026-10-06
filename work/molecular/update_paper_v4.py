"""Update the same standalone manuscript with executed interface/model checks."""
from pathlib import Path
import json
HERE=Path(__file__).resolve().parent
PROJECT=HERE.parent.parent
paper=PROJECT/'outputs/pulizia-grafene-pmma-ppc.tex'
old=(PROJECT/'work/versione-3/pulizia-grafene-pmma-ppc.tex').read_text()
def number(x):return f'{x:.2f}'.replace('.',',')
interface=[]
for name,label in [('HFIP','HFIP'),('DMF','DMF'),('THF','THF'),('AceticAcid','Acido acetico'),('Anisole','Anisolo')]:
    d=json.loads((HERE/'results/gaff-systems'/f'{name}_n4_seed20261006/interface-analysis.json').read_text())['statistics_second_half']
    interface.append(label+' & '+' & '.join(number(d[k]['mean']) for k in
        ['ester_methyl_G_LJ_kJ_mol','ester_core_G_LJ_kJ_mol','backbone_and_side_methyl_G_LJ_kJ_mol','polymer_graphene_LJ_energy_kJ_mol'])+r' \\')
bff=json.loads((HERE/'results/byteff-pol/checks.json').read_text())
force=next(x for x in bff['pairs'] if x['pair']=='HFIP__MethylAcetate')
new=old.replace('Revisione 3: calcoli quantistici, modello molecolare appreso e simulazione del liquido HFIP',
                'Revisione 4: traiettorie interfaciali e controlli dei modelli polarizzabili')
new=new.replace('6 ottobre 2026 --- revisione 3','6 ottobre 2026 --- revisione 4')
new=new.replace(r'\textbf{Stato: calcoli elettronici e controllo del liquido completati;'+
                '\n'+r'energia libera del distacco interfaciale non ancora qualificata.}',
                r'\textbf{Stato: traiettorie interfaciali eseguite su un oligomero;'+
                '\n'+r'nessuna energia libera di distacco o pulizia reale ancora qualificata.}')
new=new.replace('Aggiungiamo calcoli molecolari quantistici e con un potenziale appreso, oltre a una simulazione atomistica del liquido HFIP.',
                'Aggiungiamo calcoli quantistici, due modelli polarizzabili appresi e traiettorie di un oligomero su grafene in cinque liquidi.')
new=new.replace(r'HFIP \`e la nuova candidata chimica discussa nella Sezione~\ref{sec:molecular}.',
                r'HFIP rimane una candidata da verificare, discussa nella Sezione~\ref{sec:molecular}. I nuovi controlli non la promuovono sopra DMF.')
needle=r'\textbf{Restano aperti i controlli di sensibilit\`a 1--4 del liquido,'
start=new.index(needle)
finish=new.index(r'\section{Decisione scientifica',start)
addition=r'''
\subsection{Controllo della dimensione del modello e della rotazione di HFIP}

L'anomalia della versione piccola non viene generalizzata all'intera
famiglia. Ripetendo il controllo con MACE-POLAR-1-M, a nuclei fissi,
otteniamo diagonali positive per HFIP, DMF e THF. Per HFIP
$\alpha_{zz}=0{,}4947\,e\,\text{\AA}^2/\mathrm V$, da confrontare
con $0{,}4242$ di PBE0. La traslazione della molecola neutra non cambia
l'energia oltre $10^{-12}$\,eV. Abbiamo verificato nel codice la
convenzione: la molteplicit\`a del singoletto \`e 1, mentre il numero
di elettroni spaiati \`e zero. Anche con le geometrie centrate la
versione S d\`a diagonali negative per tutte e tre le molecole.
Il problema \`e quindi specifico a quella parametrizzazione o architettura,
non una conclusione contro tutti i potenziali appresi.

Abbiamo ruotato soltanto l'idrogeno OH in dodici posizioni, mantenendo
fissi tutti gli altri nuclei e confrontando le energie relative.
La barriera campionata \`e 15,075\,kJ/mol con
PBE0-D3(BJ)/def2-TZVP, 14,145 con MACE-POLAR-1-M e 15,320 con S.
GAFF2 d\`a 17,793\,kJ/mol, ma cambia anche la forma del profilo:
il massimo \`e a $0^\circ$, mentre nei riferimenti quantistici
\`e vicino a $\pm90^\circ$. Il controllo riguarda la molecola isolata;
non fornisce la popolazione dei conformeri nel liquido o la sua
permittivit\`a. Il modello M supera questi controlli molecolari,
senza che ne segua una qualificazione di grafene o \SiO.

\subsection{Traiettorie effettive del polimero sulla superficie}

Sono state eseguite cinque traiettorie a 298,15\,K con quattro unit\`a
PMMA, 65 atomi, e un foglio periodico rigido di 448 atomi di carbonio.
I parametri PMMA e grafene vengono dal supplemento di
Behbahani e Harmandaris \cite{Behbahani2021}; i solventi
da GAFF2/AM1-BCC, con neutralit\`a verificata. La cella misura
$3{,}4032\times3{,}4385\times3{,}8$\,nm. Seguiamo 50\,ps iniziali
a volume fisso e 100\,ps di osservazione per liquido.

\textbf{Il sistema contiene liquido su entrambe le facce del foglio:
non simula ancora grafene appoggiato sull'ossido.} L'assenza
di \SiO\ limita sia la selettivit\`a sia l'affinit\`a quantitativa:
il substrato reale pu\`o modificare il campo e le interazioni attraverso
il monostrato. Il grafene viene tenuto fisso; la sua immobilit\`a
verifica il vincolo numerico, non l'assenza di danni del trattamento.
La catena corta non rappresenta PMMA 950K, un deposito multistrato,
una frazione covalente o un residuo sepolto.

Nessuna di queste traiettorie produce distacco completo spontaneo.
La Tabella~\ref{tab:interface} scompone l'interazione diretta
polimero--grafene, mediata sulla seconda met\`a dei 100\,ps.
\textbf{Questa energia non \`e l'energia libera di adsorbimento.}
Non include l'intero bilancio polimero--liquido e grafene--liquido,
n\'e l'entropia o una dimostrazione di equilibrio.

\begin{table}[htbp]\centering\small
\caption{Energie dirette polimero--grafene nel modello con quattro
unit\`a PMMA, in kJ/mol di oligomero. Una sola realizzazione per
liquido; medie temporali, non una classifica di pulizia.}
\label{tab:interface}
\begin{tabular}{lrrrr}\toprule
Liquido & Metile estere & Nucleo estere & Scheletro/metili & Totale\\\midrule
INTERFACE_ROWS
\bottomrule\end{tabular}
\end{table}

In anisolo i gruppi estere perdono molti contatti, mentre lo scheletro
rimane aderente. In HFIP il residuo aderente forma gi\`a circa
2,9 legami idrogeno con i carbonili, contro circa 3,5 nello stato
distaccato ottenuto con un vincolo di separazione. Questa osservazione
del modello \`e coerente con una solvatazione gi\`a presente nello
stato aderente: contare legami idrogeno non basta a stabilire il vantaggio
del distacco. Le distribuzioni conformazionali e la loro convergenza
rimangono da controllare. L'anisolo \`e inoltre il solvente del
PMMA 950K commerciale; ci\`o motiva il controllo, ma la solubilit\`a
del materiale libero non dimostra il recupero del campione.

\subsection{Separazione controllata e motivo del mancato profilo di energia libera}

Per HFIP, DMF e anisolo abbiamo eseguito tre finestre di 125\,ps,
con centri a 0,60, 0,90 e 1,35\,nm dal foglio e molla
$K=2000$\,kJ\,mol$^{-1}$\,nm$^{-2}$. A 1,35\,nm la seconda
met\`a di ciascun campione non ha contatti stretti con il foglio,
e l'interazione diretta residua \`e trascurabile. \`E un distacco
forzato di un oligomero, non una pulizia ottenuta da un bagno.
La forza restituita dal programma coincide con quella ricostruita
da $K(z_0-\langle z\rangle)$.

Le tre distribuzioni di distanza non hanno sovrapposizione osservata
fra finestre adiacenti. Non si possono collegare i loro livelli
di energia libera senza ulteriori campioni. Inoltre, nella finestra
HFIP a 1,35\,nm la forza media passa da circa $+134$ a $-70$
kJ\,mol$^{-1}$\,nm$^{-1}$ fra blocchi iniziale e finale di 25\,ps:
\`e una deriva, non una differenza di equilibrio rispetto a DMF.
\textbf{Non viene quindi calcolato un costo di distacco dai tre punti.}

\`E avviata una campagna pi\`u fitta di undici finestre,
500\,ps ciascuna, con $K=750$ e centri da 0,45 a 1,45\,nm,
per HFIP e DMF. Le finestre completate migliorano la sovrapposizione
nella regione aderente; il percorso inverso, altre conformazioni e
la convergenza delle variabili nascoste sono ancora necessari.
Il codice che ricostruisce le distribuzioni \`e stato verificato
su campioni indipendenti di un potenziale armonico noto:
errore quadratico medio del profilo 0,068\,kJ/mol.
Questo qualifica il controllo numerico, non le traiettorie del materiale.
Nessun tempo di pulizia sperimentale viene letto dai percorsi accelerati.

\subsection{Un generatore polarizzabile di parametri e il limite della verifica nel vuoto}

Abbiamo eseguito ByteFF-Pol \cite{ByteFFPol2026}, che predice
parametri molecolari con una rete e poi usa un potenziale esplicito
con polarizzazione. Sono stati generati sette insiemi di parametri
e verificati tredici complessi sulle stesse geometrie degli altri
confronti. Le energie sono state scomposte e ricomposte;
la forza netta sul frammento coincide con la derivata numerica
dell'energia entro $2,3\times10^{-6}$\,kJ\,mol$^{-1}$\,\AA$^{-1}$
nel controllo HFIP--acetato di metile. A grande distanza
l'interazione tende a zero.

L'interazione calcolata \`e $-34,57$\,kJ/mol per HFIP--acetato,
$-28,81$ per acido acetico--acetato e $-20,86$ per IPA--acetato.
Sono valori diversi dai riferimenti PBE0 della Tabella~\ref{tab:quantum}
e da MACE-POLAR-1-M, che d\`a rispettivamente $-46,54$, $-45,10$
e $-27,15$. Il corretto controllo delle forze non risolve
questa incertezza fisica. \`E avviato un riferimento con il
funzionale e la base dell'addestramento, $\omega$B97M-V/def2-TZVPD,
includendo VV10 senza aggiungere D3, e controllando la griglia
d'integrazione. Finch\'e quel controllo non \`e completato,
la differenza non viene attribuita interamente al modello.

Il modello non include silicio. La sua esecuzione nel vuoto
non dimostra densit\`a, permittivit\`a o desorbimento corretti.
L'esportazione verso OpenMM richiede inoltre le scale elettrostatiche
intramolecolari specifiche della pubblicazione; la versione installata
non \`e stata modificata globalmente e non viene usata come
se fosse gi\`a quel motore. Questa \`e una via da qualificare,
non una soluzione dedotta dal nome del programma.

\subsection{Vie indipendenti e lettura critica delle nuove proposte}

Abbiamo tentato anche la teoria statistica del liquido DRISM/KH,
imponendo la permittivit\`a sperimentale. I modelli organici
non raggiungono il criterio di convergenza, neppure nei tentativi
documentati di variazione del solutore e prosecuzione in densit\`a
e carica. L'acqua di controllo distribuita con il programma
converge in 27 iterazioni. Non usiamo suscettivit\`a o energie
di solvatazione organiche non convergenti.

Il trattamento con nitrito del 2025 \cite{LeeNitrite2025} avviene
prima del deposito, lascia particelle e riduce il rapporto XPS
C--O--C da 20\% a 17\%; non prova il salvataggio richiesto.
Un trattamento elettrico senza bolle \cite{Park2018} \`e invece
applicato dopo il trasferimento. Richiede per\`o contatto elettrico
al grafene e un elettrolita: resta fuori dalla pulizia con soli bagni
e non viene presentato come immediatamente compatibile con
cristalli isolati da prelevare. Un miglioramento elettrico deve
comunque essere distinto dalla scomparsa del polimero.

\textbf{La ricerca computazionale \`e attiva. Restano aperti il
campionamento interfaciale, la qualificazione dei solventi,
PPC e lunghezze maggiori, il substrato reale e la verifica sul campione.}
I risultati disponibili non autorizzano una classifica dei pulenti
o una dichiarazione di soluzione.

'''
addition=addition.replace('INTERFACE_ROWS','\n'.join(interface))
new=new[:start]+addition+new[finish:]
new=new.replace(r'\appendix',r'''Il nuovo risultato da cercare \`e una differenza di energia libera
e di barriera che resista ai controlli di modello e campionamento,+insieme a una finestra di selettivit\`a sul substrato reale. La scelta
dei bagni deve poi superare le prove di massa, topografia, chimica
e resa del prelievo gi\`a specificate. La campagna molecolare non
ha ancora raggiunto questo risultato.

\appendix'''.replace('campionamento,\\+insieme','campionamento,\ninsieme'))
refs=r'''
\bibitem{ByteFFPol2026}
T. Zheng, X. Xu, Z. Wang et al.,
\emph{Bridging quantum mechanics to liquid properties via a universal
organic force field}, Nature Communications \textbf{17} (2026), 8276.
\href{https://doi.org/10.1038/s41467-026-73566-3}{DOI: 10.1038/s41467-026-73566-3}.

\bibitem{LeeNitrite2025}
K. Lee, J. Kil, J. Park, S. Yang e B. Park,
\emph{Rapid and Efficient Polymer/Contaminant Removal from Single-Layer
Graphene via Aqueous Sodium Nitrite Rinsing for Enhanced Electronic Applications},
Polymers \textbf{17} (2025), 689.
\href{https://doi.org/10.3390/polym17050689}{DOI: 10.3390/polym17050689}.

\bibitem{Park2018}
B. Park, J. N. Huh, W. S. Lee e I.-G. Bae,
\emph{Simple and rapid cleaning of graphenes with a `bubble-free'
electrochemical treatment}, Journal of Materials Chemistry C
\textbf{6} (2018), 2234--2244.
\href{https://doi.org/10.1039/C7TC05695H}{DOI: 10.1039/C7TC05695H}.

'''
new=new.replace(r'\end{thebibliography}',refs+r'\end{thebibliography}')
# Reject unresolved template markers and a known accidental string splice.
new=new.replace('campionamento,+insieme', 'campionamento,\ninsieme')
new=new.replace(r'campionamento,\+insieme','campionamento,\ninsieme')
assert 'INTERFACE_ROWS' not in new and r'\+insieme' not in new
assert force['component_sum_independent_error_kJ_mol']<1e-6
paper.write_text(new)
print('Revision4 source saved, same standalone editor path:',paper)
