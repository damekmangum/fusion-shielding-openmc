# %%
import os
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
import openmc

openmc.config['cross_sections'] = '/home/jovyan/nuclear_data/cross_sections.xml'

print("OpenMC version:", openmc.__version__)
print("Cross sections:", openmc.config.get('cross_sections'))

os.makedirs('runs', exist_ok=True)
os.makedirs('figures', exist_ok=True)

# %%
# Materials

steel = openmc.Material(name = 'stainless steel')
steel.set_density('g/cm3', 7.95)
steel.add_element('Fe', 0.655, 'wo')
steel.add_element('Cr', 0.170, 'wo')
steel.add_element('Ni', 0.120, 'wo')
steel.add_element('Mo', 0.025, 'wo')
steel.add_element('Mn', 0.020, 'wo')
steel.add_element('Si', 0.010, 'wo')

wc = openmc.Material(name = 'tungsten carbide')
wc.set_density('g/cm3', 15.63)
wc.add_element('W', 1, 'ao')
wc.add_element('C', 1, 'ao')

b4c = openmc.Material(name = 'boron carbide')
b4c.set_density('g/cm3', 2.52)
b4c.add_element('B', 4, 'ao')
b4c.add_element('C', 1, 'ao')

copper = openmc.Material(name = 'copper')
copper.set_density('g/cm3', 8.96)
copper.add_element('Cu', 1.0)

magnet = openmc.Material.mix_materials(
    [steel, copper], [0.6, 0.4], 'vo', name = 'magnet (60/40 Steel/Copper)'
)

shield_materials = {'Steel': steel, 'WC': wc, 'B4C': b4c}

# %%
PLASMA_R = 50.0         # 50 cm radius, void around point source
VV_CM = 5.0             # 5 cm thick vacuum vessel
MAGNET_FRONT_CM = 1.0   # 1 cm thick slice where results are tallied
MAGNET_BACK_CM = 9.0    # 9 cm, rest of magnet - backscatter

def build_model(shield_mat, shield_cm, particles=100000, batches=10, survival_biasing=True):
    # Layer from the center
    layers = [('plasma', None, PLASMA_R), ('vessel', steel, VV_CM)]
    if shield_cm > 0:
        layers.append(('shield', shield_mat, shield_cm))
    layers += [('magnet_front', magnet, MAGNET_FRONT_CM),
               ('magnet_back', magnet, MAGNET_BACK_CM)]

# Nested Spheres
    cells, inner, r = {}, None, 0.0
    radii = {}
    for i, (name, mat, t) in enumerate(layers):
        r_in, r = r, r + t
        radii[name] = (r_in, r)
        bt = 'vacuum' if i == len(layers) - 1 else 'transmission'
        surf = openmc.Sphere(r = r, boundary_type = bt)
        region = -surf if inner is None else (+inner & -surf)
        cells[name] = openmc.Cell(name = name, fill = mat, region = region)
        inner = surf

    geometry = openmc.Geometry(list(cells.values()))
    materials = openmc.Materials({c.fill for c in cells.values() if c.fill is not None})

# 14.1 MeV isotropic point source at the center
    source = openmc.IndependentSource()
    source.space = openmc.stats.Point((0.0, 0.0, 0.0))
    source.angle = openmc.stats.Isotropic()
    source.energy = openmc.stats.Discrete([14.1e6], [1.0])
    source.particle = 'neutron'

    settings = openmc.Settings()
    settings.run_mode = 'fixed source'
    settings.source = source
    settings.particles = particles
    settings.batches = batches
    settings.survival_biasing = survival_biasing

# Tallies on front slice of magnet
    cell_filter = openmc.CellFilter(cells['magnet_front'])

    fast_tally = openmc.Tally(name = 'fast_flux')
    fast_tally.filters = [cell_filter, openmc.EnergyFilter([0.0, 1.0e5, 2.0e7])]
    fast_tally.scores = ['flux']

    spec_tally = openmc.Tally(name = 'spectrum')
    spec_tally.filters = [cell_filter, openmc.EnergyFilter.from_group_structure('CCFE-709')]
    spec_tally.scores = ['flux']

    tallies = openmc.Tallies([fast_tally, spec_tally])

# Volume of tally slice -- for normalization
    r_in, r_out = radii['magnet_front']
    volume = 4.0 / 3.0 * np.pi * (r_out**3 - r_in**3)

    return openmc.Model(geometry, materials, settings, tallies), volume

# %%
def run_case(shield_name, shield_cm, **kwargs):
    model, volume = build_model(shield_materials[shield_name], shield_cm, **kwargs)
    workdir = Path('runs') / f'{shield_name}_{shield_cm}cm'
    workdir.mkdir(parents = True, exist_ok = True)

    sp_path = Path(model.run(cwd = workdir, output = False))
    if not sp_path.is_absolute():
        sp_path = workdir / sp_path

    with openmc.StatePoint(sp_path) as sp:
        fast = sp.get_tally(name = 'fast_flux')
        mean, std = fast.mean.ravel(), fast.std_dev.ravel()     # [slow bin, fast bin]

        spec = sp.get_tally(name = 'spectrum')
        edges = spec.find_filter(openmc.EnergyFilter).values
        s_mean = spec.mean.ravel() / volume
        s_std = spec.std_dev.ravel() / volume

    fast_flux = mean[1] / volume
    fast_rel = std[1] / mean[1] if mean[1] > 0 else np.nan
    return {'fast': fast_flux, 'fast_rel': fast_rel,
            'spectrum': s_mean, 'spectrum_std': s_std, 'edges': edges}
    

# %%
# Check
base = run_case('Steel', 0) # No shield, material doesn't matter
test = run_case('WC', 20)

print(f"No shield : {base['fast']:.3e} n/cm^2 per source neutron (rel. error {base['fast_rel']:.1%})")
print(f"20 cm WC  : {test['fast']:.3e} n/cm^2 per source neutron (rel. error {test['fast_rel']:.1%})")
print(f"Reduction : {base['fast'] / test['fast']:.1f}x")

# %%
thicknesses = [0, 10, 20, 30]

results = {name: {0: base} for name in shield_materials}
for name in shield_materials:
    for t in thicknesses[1:]:
        r = run_case(name, t)
        results[name][t] = r
        print(f"{name:5s} {t:3d} cm : fast flux {r['fast']:.3e} (rel. error {r['fast_rel']:.1%})")


# %%
# Empty gap of the same thickness: separates distance from shielding
shield_materials['Void'] = None
try:
    void = {t: run_case('Void', t) for t in thicknesses[1:]}
finally:
    del shield_materials['Void']

print("Reduction vs. an empty gap of the same thickness:")
for name in shield_materials:
    factors = [void[t]['fast'] / results[name][t]['fast'] for t in thicknesses[1:]]
    print(f"{name:5s}: " + ", ".join(f"{t} cm: {f:.1f}x" for t, f in zip(thicknesses[1:], factors)))

# %%
plt.figure(figsize=(7, 5))
for name in shield_materials:
    flux = np.array([results[name][t]['fast'] for t in thicknesses])
    err = flux * np.array([results[name][t]['fast_rel'] for t in thicknesses])
    plt.errorbar(thicknesses, flux, yerr = err, marker = 'o', capsize = 3, label = name)
plt.yscale('log')
plt.xlabel('Shield thickness (cm)')
plt.ylabel('Fast flux at magnet, >0.1 MeV\n(n/cm^2 per source neutron)')
plt.title('Fast neutron flux at the magnet vs shield thickness')
plt.grid(True, which = 'both', alpha = 0.3)
plt.legend()
plt.tight_layout()
plt.savefig('figures/fast_flux_vs_thickness.png', dpi = 200)
plt.show()

# %%
T_COMPARE = 20
names = list(shield_materials)
factors = [void[T_COMPARE]['fast'] / results[n][T_COMPARE]['fast'] for n in names]

plt.figure(figsize=(6, 4))
plt.bar(names, factors)
plt.ylabel('Fast-flux reduction factor\n(vs. empty gap of same thickness)')
plt.title(f'Shield performance at {T_COMPARE} cm')
plt.tight_layout()
plt.savefig(f'figures/reduction_at_{T_COMPARE}cm.png', dpi = 200)
plt.show()

for n, f in zip(names, factors):
    print(f"{n:5s}: {f:.1f}x reduction")

# %%
def plot_spectrum(T):
    plt.figure(figsize=(8, 5))
    for label, r in [('No shield', base)] + [(n, results[n][T]) for n in names]:
        e = r['edges']
        lethargy = np.log(e[1:] / e[:-1])
        plt.step(e[:-1], r['spectrum'] / lethargy, where='post',
                 label=label if label == 'No shield' else f'{label}, {T} cm')
    plt.xscale('log')
    plt.yscale('log')
    plt.xlim(1e-1, 2e7)
    plt.ylim(1e-9, 2e-3)
    plt.xlabel('Neutron energy (eV)')
    plt.ylabel('Flux per unit lethargy\n(n/cm² per source neutron)')
    plt.title(f'Neutron spectrum at the magnet, {T} cm shields')
    plt.axvline(1e5, color='gray', ls='--', lw=1)
    plt.grid(True, which='both', alpha=0.3)
    plt.legend()
    plt.tight_layout()
    plt.savefig(f'figures/spectrum_{T}cm.png', dpi=200)
    plt.show()

plot_spectrum(10)
plot_spectrum(20)

# %%
def flux_above(r, E_min):
    e = r['edges']
    return r['spectrum'][e[:-1] >= E_min].sum()

for E in [1e5, 1e6]:
    print(f"Flux above {E:.0e} eV")
    for name in shield_materials:
        print(f"  {name:5s}", [f"{flux_above(results[name][t], E):.2e}" for t in thicknesses])


