# CrossFlag demo — offline, committed-data-only build.
# Requires the `crossflag-spike` conda env (+ matplotlib, pytest); see demo/README.md.

CONDA_RUN = conda run -n crossflag-spike

.PHONY: setup demo test-demo clean-demo

## setup: one-time — install this package (editable) + demo deps into the env.
setup:
	$(CONDA_RUN) pip install -e .
	$(CONDA_RUN) pip install matplotlib pytest

## demo: regenerate demo/dashboard.html + figures + verdict table from committed data.
demo:
	$(CONDA_RUN) python -m crossflag.demo.build

## test-demo: run the demo acceptance tests (Phase 5 gate).
test-demo:
	$(CONDA_RUN) python -m pytest tests/test_demo.py -q

## clean-demo: remove generated demo artifacts (source stays).
clean-demo:
	rm -f demo/dashboard.html demo/verdict_table.json demo/verdict_table.md
	rm -f demo/figures/*.png
