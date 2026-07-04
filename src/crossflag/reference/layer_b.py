"""Build Layer B (autoimmune focus) of the CrossFlag reference set.

TODO: NOT YET IMPLEMENTED. This is the AAgAtlas autoantigen layer described in
docs/reference-set.md §"Layer B — Autoimmune focus (AAgAtlas)". Today the merged
`self_proteins.csv` contains Layer A + the (small) Layer C set only; Layer B is a
known gap.

Planned contract (per reference-set.md):
  1. Download the AAgAtlas gene/disease table (v1: ~1,126 autoantigens).
  2. Resolve each gene symbol -> UniProt ID -> canonical sequence (UniProt REST).
  3. Emit rows with layer_B=True and populated `autoimmune_conditions`.
  4. merge.py folds these into existing Layer-A rows by UniProt ID, or appends
     intracellular autoantigens as new rows with layer_A=False.

Output (when implemented): paths.BUILD / "self_proteins_layer_b.csv".
"""

from . import paths


def main() -> None:
    raise NotImplementedError(
        "Layer B (AAgAtlas autoantigens) is not implemented yet. "
        "See docs/reference-set.md §Layer B for the build contract. "
        f"Output target: {paths.BUILD / 'self_proteins_layer_b.csv'}"
    )


if __name__ == "__main__":
    main()
