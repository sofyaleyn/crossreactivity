"""Orchestrate the full reference-set build (HANDOFF step 5).

Runs the layer modules in the order from docs/reference-set.md §"Order of
operations", then (once the embedding modules land) builds the ESM-2 index.

    Layer A  ->  anchor injection  ->  background  ->  Layer B  ->  Layer C  ->  merge  ->  embed

Each stage reads/writes the locations declared in `paths.py`. Stages that hit
the network (Layer A UniProt fetch, Layer B AAgAtlas) cache into
`data/reference/raw/` and `data/reference/build/`.

Run:  python -m crossflag.reference.build_set
"""

from . import background, layer_a, layer_b, layer_c, merge, paths


def main(*, include_layer_b: bool = False) -> None:
    paths.ensure_dirs()

    # 1. Layer A base (SURFY + CSPA) with anchor injection; writes build/self_proteins_layer_a.csv
    layer_a.main()

    # 2. Background calibration set sampled from Layer A
    background.main()

    # 3. Layer B autoantigens (AAgAtlas) — not implemented yet; opt-in until it lands
    if include_layer_b:
        layer_b.main()

    # 4. Layer C mimicry seed -> fill sequences -> build/layer_c_mimicry_seed_filled.csv
    layer_c.main()

    # 5. Merge layers into the final data/reference/self_proteins.csv
    merge.main()

    # 6. TODO: build the ESM-2 embedding index over self_proteins.csv into paths.INDEX
    #    (depends on crossflag.embed.antigen, not yet implemented).
    print("Reference set built. Embedding index step pending (crossflag.embed.antigen).")


if __name__ == "__main__":
    main()
