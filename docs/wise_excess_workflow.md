# WISE W3/W4 Excess Workflow

PHOEBE is useful here as a binary-photosphere model, not as a disk/free-free classifier. It can estimate how much flux the two stars should contribute in a passband, but it does not include a debris disk or optically thin free-free component as a physical emission source.

For the local AllWISE row of `J071951.40-240400.6`, the practical diagnostic is:

1. Fit the stellar mid-infrared baseline with W1 and W2.
2. Compare W3 and W4 to that baseline.
3. If W1/W2 are photospheric but W3/W4 are high, the excess is disk-like.
4. If the excess starts already in K/W1/W2 and follows a smooth power law, free-free emission becomes more plausible.
5. Confirm with quality flags, image inspection, longer-wavelength photometry, and spectra if available.

Run:

```bash
python3 src/wise_excess_diagnostic.py
```

Outputs are written to:

```text
results/test/wise_excess/
```

The first-pass result for `J071951.40-240400.6` is debris-disk-like: W1/W2 sit on the stellar baseline, W3 is a significant excess, and W4 is a very strong excess. The W3/W4 residual ratio is consistent with a cool thermal component around 158 K.

This does not prove a debris disk. W3/W4 alone cannot uniquely separate a cool dust blackbody from every possible free-free model. The important result is narrower: the local WISE photometry does not look like ordinary free-free continuum beginning in the near infrared; it looks like a separate cool infrared component.
