import csv, statistics as st
P, R = 9.8, 0.55
def hit(p, r): return p < P and r >= R

def f(x):
    try: return float(x)
    except (TypeError, ValueError): return None

rows = []
for x in csv.DictReader(open('data/results/screen/screen_metrics.csv')):
    pae, rep = f(x['PAE_IF']), f(x['epitope_reprod'])
    if pae is None or rep is None or not x['antigen_len'].isdigit():
        continue
    x['len'] = int(x['antigen_len']); x['pae'] = pae; x['rep'] = rep; x['H'] = hit(pae, rep)
    rows.append(x)

fl = [r for r in rows if r['job'] == 'flagship']
comm = {r['run']: (f(r['PAE_IF_mean']), f(r['epitope_reprod']))
        for r in csv.DictReader(open('data/results/cofold_metrics.csv'))
        if f(r['PAE_IF_mean']) is not None and f(r['epitope_reprod']) is not None}
fzd5, ulbp2 = comm['cofold-fzd5'], comm['cofold-ulbp2']
def band(d, lo, hi): return [r for r in d if lo <= r['len'] < hi]

print(f"=== CHECK A: Flagship FPR by antigen length (n={len(fl)}) ===")
for lo, hi in [(0,80),(80,150),(150,300),(300,600),(600,1001)]:
    b = band(fl, lo, hi)
    if b: print(f"  {lo:4d}-{hi:<4d} aa: {len(b):3d} proteins  FPR {100*sum(r['H'] for r in b)/len(b):5.1f}%")
valid = [r for r in fl if r['len'] >= 150]
print(f"  >>> VALID REGIME (>=150aa): FPR {100*sum(r['H'] for r in valid)/len(valid):.1f}%  ({sum(r['H'] for r in valid)}/{len(valid)})")

print("\n=== CHECK B: FZD5/ULBP2 enrichment WITHIN size band (150-300aa) ===")
vb = band(fl, 150, 300)
for name, (pae, rep) in [('FZD5', fzd5), ('ULBP2', ulbp2)]:
    br = sum(1 for r in vb if r['rep'] < rep); bp = sum(1 for r in vb if r['pae'] > pae)
    print(f"  {name:6s} reprod {rep:.3f} > {br}/{len(vb)} same-size decoys ({100*br/len(vb):.0f}%);  "
          f"PAE {pae:.2f} better than {100*bp/len(vb):.0f}%;  hit={hit(pae, rep)}")

print("\n=== CHECK C: Anchor-2 PF4 (70aa) — LENGTH-CONTROLLED ===")
a2 = [r for r in rows if r['job'] == 'anchor2']
pf4 = [r for r in a2 if r['gene'] == 'PF4'][0]; dec = [r for r in a2 if r['gene'] != 'PF4']
sd = band(dec, 40, 110)
print(f"  PF4: PAE {pf4['pae']:.2f}  reprod {pf4['rep']:.3f}  hit={pf4['H']}")
print(f"  PF4-sized (40-110aa) anchor-2 decoys: n={len(sd)}  FPR {100*sum(r['H'] for r in sd)/len(sd):.1f}%")
print(f"  PF4 reprod > {sum(1 for r in sd if r['rep'] < pf4['rep'])}/{len(sd)} same-size decoys "
      f"({100*sum(1 for r in sd if r['rep'] < pf4['rep'])/len(sd):.0f}%);  "
      f"PAE better than {100*sum(1 for r in sd if r['pae'] > pf4['pae'])/len(sd):.0f}%")

print("\n=== CHECK D: antibody stickiness in VALID regime (>=150aa) ===")
for job, ab in [('flagship', 'SHR-1210'), ('anchor2', 'ABT-736')]:
    v = [r for r in rows if r['job'] == job and r['len'] >= 150 and r['gene'] != 'PF4']
    print(f"  {ab:9s} FPR(>=150aa) {100*sum(r['H'] for r in v)/len(v):.1f}%  ({sum(r['H'] for r in v)}/{len(v)})")

print("\n=== CHECK E: non-physical over-docking (PAE < real off-target FZD5=5.69) ===")
npf = [r for r in fl if r['pae'] < 5.69]
print(f"  flagship decoys PAE<5.69: {len(npf)}/{len(fl)} ({100*len(npf)/len(fl):.0f}%);  median len {st.median([r['len'] for r in npf])}")
print(f"  of those, >=150aa (real concern): {sum(1 for r in npf if r['len'] >= 150)}")

print("\n=== TOP FLAGSHIP HITS in VALID regime (>=150aa) — candidate novel off-targets ===")
vh = sorted([r for r in fl if r['len'] >= 150 and r['H']], key=lambda r: r['pae'])[:12]
for r in vh:
    print(f"  {r['gene']:12s} {r['len']:4d}aa  PAE {r['pae']:.2f}  reprod {r['rep']:.3f}")
