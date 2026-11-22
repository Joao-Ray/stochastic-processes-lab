# Data and code for finite-window Ornstein–Uhlenbeck drift inference

Hangrui Jiao · Ocean University of China

[Download the simulation data and analysis code](OU-simulation-data-and-code.zip)

The archive accompanies *Finite-window uncertainty in Ornstein–Uhlenbeck drift inference: a reproducible benchmark*. It contains synthetic observations where saved, per-attempt inference results including failures, and the code and seeds for reproducing the experiments.

| Study | Records in this version |
| --- | --- |
| Original benchmark and confirmation | 22,000 and 2,000 |
| Three additional benchmark and confirmation batches | 72,000 scenario records from 24,000 distinct latent trajectories |
| Matched stationary-distribution comparison | 18,000 reused design records from 12,000 distinct latent trajectories |
| Original stationary-rank and measurement-error supplements | 16,000 and 9,000 |
| Complete rank-cell study | 500 reused outer paths, 251 cells per path |

Several studies use the same observations. These row counts should not be summed as independent trajectories. The ZIP includes a file guide, reproduction commands, per-file checksums, and the article's 500 complete rank sets. The simulations use no personal or observational data.

The [archive checksum](SHA256.txt) identifies this version. After extraction, run `shasum -a 256 -c SHA256SUMS.txt` to check the individual files.
