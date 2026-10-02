# P1b supervisor figures

Start with [P1b_supervisor_figures.pdf](P1b_supervisor_figures.pdf): a six-page, landscape figure set. Individual figures are available as 200 dpi PNGs for slides and vector PDFs (field meshes rasterized). [figure_contact_sheet.png](figure_contact_sheet.png) previews all six pages.

| Figure | Content | Evidence type |
|---|---|---|
| 01_domain_and_question | Domain, ice, grain, physical setting and scientific question | To-scale x–y schematic from the submitted case |
| 02_initial_profiles | Temperature, liquid salinity, meltwater fraction and density vs depth | Prescribed analytical profiles; configured nonlinear EOS |
| 03_grain_field_snapshot | Temperature, salinity and tracer around the grain | Actual HDF5 slice, with time and slice position labeled |
| 04_early_grain_motion | Downward displacement and settling speed | Actual particle time series |
| 05_precursor_to_production | P1b, P2 and the conditional larger-domain calculation | Experimental plan; targets are not claims of acceptance |
| 06_matched_control_profiles | Default vs Yang salt-operator control profiles | Actual HDF5 profiles at the same saved time |

The initial capture was made on 2026-09-30 at 16:54 UTC (12:54 EDT). At that capture, the replacement grain run had only its initial full-field output; its particle series reached t = 0.0570475 (1.562 ms). The controls are compared at t = 0.4990127 (13.664 ms). The intended full observation window is t = 30 (821 ms). These are early outputs, not completed precursor results. The authoritative capture times and source paths are in [data/manifest.json](data/manifest.json), which is updated by `--refresh`.

## Interpretation

- Current default-operator cases: grain job 20977443 and matched ice-only job 20977444. Yang ice-only job 20973341 supplies the additional operator reference. The failed Yang grain run is not used in these figures.
- Liquid volume fraction is written as C_L. The grain field plots show pure-liquid cells only (C_L >= 0.999 and C_S <= 0.001). Grey includes the diffuse band and sediment support as well as ice and the grain. The orange contour is C_S = 0.5, not an independently reconstructed geometric particle boundary.
- The initial-profile figure shows the prescribed liquid profiles before diffuse-band weighting. It is not a measured evolved profile. Its density curve is evaluated with the configured nonlinear EOS, relative to the far-field reference state.
- The control salt plot compares default stored s with Yang C_L S: salt amount per cell volume in reference-salinity units. Neither curve is divided by C_L; this plot is not the liquid salinity at the diffuse interface. The tiny Yang regularization capacity is not included in this physical salt-inventory comparison.
- No background current is imposed. The three-equation interface state initializes a freely evolving sublayer; the reference melt rate is not a maintained boundary flux.
- Sc = 70 uses a simulated diffusivity about 26 times the stated physical salt diffusivity. Background-control agreement cannot establish insensitivity of grain-driven transport to the salt operator. Larger-domain dimensions remain to be selected after the precursor checks.

## Reproduce or refresh

Run from this directory:

```bash
python build_figures.py            # redraw from saved small data extracts
python build_figures.py --refresh  # read current outputs, replace extracts and redraw
```

Dependencies: numpy, scipy, matplotlib, h5py and Pillow. The Anvil anaconda 2025.06 Python available in this session provides them.

The script reads simulation files without modifying them. It captures only small slices and profiles, the particle series, the staged case record and a source manifest in `data/`. Replotting the saved capture does not require access to scratch outputs. No full 3-D field is copied or loaded into memory. Refreshing replaces this figure set and its extracts; archive the directory first if you need to preserve a particular capture.
