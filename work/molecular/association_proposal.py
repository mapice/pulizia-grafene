"""Rigid molecular relocation with an exactly normalized association proposal.

Independent mixture of uniform poses and an OO-shell/OH-cone proposal.
The Hastings ratio counts ALL eligible acceptors, including overlaps.
No energy bias is added. This geometry module evaluates no force field.
"""
import numpy as np


def unit_sphere(rng):
    v=rng.normal(size=3);return v/np.linalg.norm(v)


def rotation_axis(axis,angle):
    a=axis/np.linalg.norm(axis)
    K=np.array([[0.,-a[2],a[1]],[a[2],0.,-a[0]],[-a[1],a[0],0.]])
    return np.eye(3)+np.sin(angle)*K+(1-np.cos(angle))*(K@K)


def rotation_between(v,w):
    v=v/np.linalg.norm(v);w=w/np.linalg.norm(w);c=float(np.clip(np.dot(v,w),-1.,1.))
    axis=np.cross(v,w);norm=np.linalg.norm(axis)
    if norm<1e-12:
        if c>0:return np.eye(3)
        helper=np.eye(3)[int(np.argmin(np.abs(v)))];axis=np.cross(v,helper)
        return rotation_axis(axis,np.pi)
    return rotation_axis(axis/norm,np.arctan2(norm,c))


def cone_direction(axis,cosine_min,rng):
    axis=axis/np.linalg.norm(axis);helper=np.eye(3)[int(np.argmin(np.abs(axis)))]
    e=np.cross(axis,helper);e/=np.linalg.norm(e);f=np.cross(axis,e)
    c=float(rng.uniform(cosine_min,1.));phi=float(rng.uniform(0.,2*np.pi))
    return c*axis+np.sqrt(max(0.,1-c*c))*(np.cos(phi)*e+np.sin(phi)*f)


def eligible_acceptors(xyz,length,molecule,hydrogen=10,rlo=2.6,rhi=3.3,cosine_min=np.cos(np.pi/6)):
    assert 0<rlo<rhi<length/2
    oxygen=xyz[molecule,0];axis=xyz[molecule,hydrogen]-oxygen
    axis/=np.linalg.norm(axis)
    delta=xyz[:,0]-oxygen;delta-=length*np.rint(delta/length)
    r=np.linalg.norm(delta,axis=1);valid=(r>rlo)&(r<rhi)
    valid[molecule]=False
    ids=np.flatnonzero(valid)
    return ids[(delta[ids]@axis)/r[ids]>cosine_min]


def density_without_Haar_constant(count,n,length,bias,rlo,rhi,cosine_min):
    shell=4*np.pi*(rhi**3-rlo**3)/3
    cone=(1-cosine_min)/2
    return (1-bias)/length**3+bias*count/((n-1)*shell*cone)


def propose(xyz,length,molecule,rng,hydrogen=10,bias=.5,rlo=2.6,rhi=3.3,cosine_min=np.cos(np.pi/6)):
    assert 0<bias<1 and len(xyz)>1 and 0<rlo<rhi<length/2
    old_count=len(eligible_acceptors(xyz,length,molecule,hydrogen,rlo,rhi,cosine_min))
    old=xyz[molecule];old_axis=old[hydrogen]-old[0];old_axis/=np.linalg.norm(old_axis)
    biased=bool(rng.random()<bias)
    if biased:
        targets=np.delete(np.arange(len(xyz)),molecule);target=int(rng.choice(targets))
        radial=unit_sphere(rng);r=float(rng.uniform(rlo**3,rhi**3))**(1/3)
        oxygen=xyz[target,0]+r*radial
        axis=cone_direction(-radial,cosine_min,rng)
    else:
        oxygen=rng.uniform(0.,length,size=3);axis=unit_sphere(rng);target=None
    # Haar rotation: uniform direction (or cone) plus independent uniform
    # roll. The oxygen-anchored translation has the same SE(3) Jacobian as
    # center-of-mass translation; internal geometry remains unchanged.
    R=rotation_axis(axis,float(rng.uniform(0.,2*np.pi)))@rotation_between(old_axis,axis)
    new=xyz.copy();new[molecule]=(old-old[0])@R.T+oxygen
    new_count=len(eligible_acceptors(new,length,molecule,hydrogen,rlo,rhi,cosine_min))
    if biased:assert target in eligible_acceptors(new,length,molecule,hydrogen,rlo,rhi,cosine_min)
    qold=density_without_Haar_constant(old_count,len(xyz),length,bias,rlo,rhi,cosine_min)
    qnew=density_without_Haar_constant(new_count,len(xyz),length,bias,rlo,rhi,cosine_min)
    return new,float(np.log(qold)-np.log(qnew)),dict(biased_component=biased,old_acceptors=old_count,
         new_acceptors=new_count,proposal_density_old=qold,proposal_density_new=qnew)
