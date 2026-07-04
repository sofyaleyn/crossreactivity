import boltz_api, inspect, pkgutil, os, json

out = []
def p(*a): out.append(" ".join(str(x) for x in a))

p("PKG PATH:", os.path.dirname(boltz_api.__file__))
p("VERSION:", getattr(boltz_api, "__version__", "?"))
p("\n=== top-level exports ===")
p([n for n in dir(boltz_api) if not n.startswith('_')])

from boltz_api import Boltz
p("\n=== Boltz client attributes ===")
p([n for n in dir(Boltz) if not n.startswith('_')])

p("\n=== submodules (predict/structure/binding/param) ===")
mods = []
for m in pkgutil.walk_packages(boltz_api.__path__, 'boltz_api.'):
    if any(k in m.name.lower() for k in ['predict','structure','binding','param','type','resource']):
        mods.append(m.name)
        p(" ", m.name)

# Try to locate the start/create method and its params model
import importlib
def dump_model(cls, depth=0, seen=None):
    seen = seen or set()
    if cls in seen or depth > 4: return
    seen.add(cls)
    fields = getattr(cls, "model_fields", None)
    if not fields: return
    p("  " * depth + f"MODEL {cls.__name__}:")
    for name, f in fields.items():
        ann = getattr(f, "annotation", None)
        req = getattr(f, "is_required", lambda: "?")
        p("  " * depth + f"    {name}: {ann}  required={req() if callable(req) else req}")

for modname in mods:
    try:
        mod = importlib.import_module(modname)
    except Exception as e:
        p("  import fail", modname, e); continue
    for nm, obj in vars(mod).items():
        if inspect.isclass(obj) and hasattr(obj, "model_fields"):
            try: dump_model(obj)
            except Exception as e: p("  dump fail", nm, e)

open("/private/tmp/claude-501/-Users-cheparukhin-crossreactivity/64bf821d-4583-4d7b-92ce-c89685a1d54a/scratchpad/introspect.out","w").write("\n".join(out))
print("wrote introspect.out", len(out), "lines")
