"""Independent SI check of molecular simulation unit conversions."""
import math,json
from pathlib import Path
NA=6.02214076e23
n=125;M=168.0404;rho=.8
volume_nm3=n*M/(602.214076*rho)
# Independent derivation in kg and m^3, then convert to nm^3.
mass_kg=(n*M/1000)/NA
volume_SI_nm3=(mass_kg/(rho*1000))*1e27
assert math.isclose(volume_nm3,volume_SI_nm3,rel_tol=2e-15)
assert 40<volume_nm3<50
water_density=(33.428*18.01528)/602.214076
assert .99<water_density<1.01
r=dict(HFIP_molecules=n,molar_mass_g_mol=M,density_g_cm3=rho,volume_nm3=volume_nm3,independent_SI_volume_nm3=volume_SI_nm3,box_side_nm=volume_nm3**(1/3),water_density_control_g_cm3=water_density,unit_constant='602.214076 = NA*1e-21',initial_invalid_pilot='Using 0.602214076 enlarged volume by 1000 and box by 10. Those pilots invalidated; no conclusions drawn.')
Path(__file__).with_name('results').joinpath('unit-checks.json').write_text(json.dumps(r,indent=2)+'\n')
print(r)
