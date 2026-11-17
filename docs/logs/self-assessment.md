# Critical self-assessment

The weakest link is the choice of the smoothing parameter. We report eta0 because the guide
prescribes it as the reference value, and the verdict survives four decades of it, but the
supporting scan is a sensitivity table rather than a statistical selection rule. A cross-validated
or L-curve criterion could move eta by an order of magnitude, which would change the reconstructed
amplitude by about 4 per cent and would not change the pass/fail margin of 0.0024. We did not
implement that criterion, so the map amplitude is reported as conditional on eta0.

The second weakness is the forward operator. It assumes a mobility that does not vary with carrier
density, which is a poor assumption for any disordered film and an untested one for this wafer. If
mobility varies spatially, the map conflates carrier density with transport and the recovered
pattern is only a lower bound on the true non-uniformity. Both points are addressed by adding a
density-dependent mobility term and re-running the same gates; we did not have a parameterisation
that the data could constrain.
