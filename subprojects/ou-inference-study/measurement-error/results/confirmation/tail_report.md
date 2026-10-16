# Tail sensitivity of finite drift errors

[Complete report](report.md) · [Every record](replications.csv.gz) · [Full tail audit](tails.json)

No record is removed or winsorized. Each denominator is the finite point-estimate population of that method. Missing boundary/optimizer estimates remain missing and are counted in the main report. This does not estimate unconditional risk of a procedure that always returns a finite drift.

| Design | Method | Available n | Median abs error | 99th percentile | Max abs error | Top-five share of squared error |
| --- | --- | --- | --- | --- | --- | --- |
| T=10,k=1,eta=1,stress=none | conditional | 1000 | 16.9595 | 37.9637 | 94.1424 | 0.0470 |
| T=10,k=1,eta=1,stress=none | stationary_ignorant | 994 | 16.9566 | 38.0345 | 94.2207 | 0.0474 |
| T=10,k=1,eta=1,stress=none | ml | 998 | 0.4937 | 6.2084 | 27.3633 | 0.3579 |
| T=10,k=1,eta=1,stress=none | reml | 973 | 0.4787 | 5.1174 | 26.8458 | 0.4341 |
| T=10,k=2,eta=1,stress=none | conditional | 1000 | 5.6036 | 14.9739 | 24.0000 | 0.0402 |
| T=10,k=2,eta=1,stress=none | stationary_ignorant | 992 | 5.6023 | 14.8848 | 24.0889 | 0.0409 |
| T=10,k=2,eta=1,stress=none | ml | 999 | 0.5675 | 6.0069 | 9.2559 | 0.1459 |
| T=10,k=2,eta=1,stress=none | reml | 983 | 0.5180 | 5.4058 | 8.3214 | 0.1629 |
