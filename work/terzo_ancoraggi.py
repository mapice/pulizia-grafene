"""Controllo analitico del modello di ancoraggi multipli; nessuna calibrazione.

Tassi in unita' scelte e numero di contatti illustrativo. Nessuna simulazione
atomistica o misura di PMMA/grafene e' rappresentata da questo codice.
"""
import math


def increments(m, q, koff=1.0):
    d = [0.0]*(m+1)
    d[m] = 1.0/(m*koff)
    for j in range(m-1, 0, -1):
        d[j] = (1+(m-j)*q*koff*d[j+1])/(j*koff)
    return d[1:]


def from_one(m, q, koff=1.0):
    if q == 0:
        return 1/koff
    return math.expm1(m*math.log1p(q))/(m*q*koff)


def main():
    error = 0.0
    for m in range(1, 16):
        last=0
        for q in (0., .01, .1, 1., 2., 10.):
            d=increments(m,q)
            exact=from_one(m,q)
            error=max(error,abs(d[0]-exact)/exact)
            assert abs(d[0]-exact)/exact < 1e-13
            assert sum(d)>=last
            last=sum(d)
        harmonic=sum(1/j for j in range(1,m+1))
        assert abs(sum(increments(m,0))-harmonic)<1e-14
    print('Verificati T1 chiuso, limite senza riaggancio e monotonia in q.')
    print('Massimo errore relativo ricorrenza / formula:',error)
    print('Nessun tasso del campione e nessun tempo fisico sono stati stimati.')


if __name__=='__main__':
    main()
