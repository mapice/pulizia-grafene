"""Confronto indipendente dei metodi numerici e delle formule del documento."""
import math
import json
import numpy as np
import model_cleaning as m
import model_review as r

errors=[]
for beta, renewal, release, tau, initial in [
    (0.1,0.2,0.1,2,[1,0,0.6,0]),
    (2,3,0.1,10,[1,0,0.6,0]),
    (10,1,1,3,[0.4,0.1,0.5,0]),
]:
    x=m.exp_generator(m.transport_generator(beta,renewal,release),tau)@np.array(initial)
    y,_,_=r.integrate(tau,beta,renewal,release,tuple(initial))
    errors.append(float(np.max(np.abs(x-np.array(y)))))
assert max(errors)<1e-11
cyclic=[]
for w in [0.1,1,2,3]:
    for n in [1,2,4,8,20,100]:
        d=m.conformation_generator(w,0,0)
        t=m.conformation_generator(0,0,w)
        calculated=m.cycle_residue(d,t,n)
        exact=math.exp(-w)*(1+n*(-math.expm1(-w/n)))
        cyclic.append(abs(calculated-exact))
assert max(cyclic)<1e-11
print(json.dumps({'uniformization_RK4_max_error':max(errors),'cycles_closed_form_max_error':max(cyclic),'meaning':'Verifica delle equazioni, non del materiale reale'},indent=2))
