# Tail sensitivity of finite drift errors

[Complete report](report.md) · [Every record](replications.csv.gz) · [Full tail audit](tails.json)

No record is removed or winsorized. Each denominator is the finite point-estimate population of that method. Missing boundary/optimizer estimates remain missing and are counted in the main report. This does not estimate unconditional risk of a procedure that always returns a finite drift.

| Design | Method | Available n | Median abs error | 99th percentile | Max abs error | Top-five share of squared error |
| --- | --- | --- | --- | --- | --- | --- |
| T=10,k=1,eta=0.5,stress=none | conditional | 500 | 5.9482 | 15.4617 | 18.7609 | 0.0597 |
| T=10,k=1,eta=0.5,stress=none | stationary_ignorant | 495 | 5.9777 | 15.0221 | 18.8575 | 0.0603 |
| T=10,k=1,eta=0.5,stress=none | ml | 498 | 0.4614 | 3.7902 | 5.5974 | 0.2268 |
| T=10,k=1,eta=0.5,stress=none | reml | 490 | 0.4337 | 3.4560 | 5.1554 | 0.2544 |
| T=10,k=2,eta=0.5,stress=none | conditional | 500 | 1.8908 | 6.2983 | 9.7465 | 0.1062 |
| T=10,k=2,eta=0.5,stress=none | stationary_ignorant | 496 | 1.8529 | 6.3760 | 9.8472 | 0.1048 |
| T=10,k=2,eta=0.5,stress=none | ml | 499 | 0.4708 | 2.7889 | 5.0813 | 0.1941 |
| T=10,k=2,eta=0.5,stress=none | reml | 497 | 0.4279 | 2.4982 | 4.6400 | 0.2203 |
| T=10,k=4,eta=0.5,stress=none | conditional | 500 | 0.7371 | 4.2636 | 5.9837 | 0.1530 |
| T=10,k=4,eta=0.5,stress=none | stationary_ignorant | 498 | 0.7128 | 4.3540 | 5.8645 | 0.1552 |
| T=10,k=4,eta=0.5,stress=none | ml | 499 | 0.4117 | 3.4595 | 4.6954 | 0.1934 |
| T=10,k=4,eta=0.5,stress=none | reml | 496 | 0.4015 | 3.0608 | 4.1674 | 0.2140 |
| T=10,k=1,eta=1,stress=none | conditional | 500 | 17.0347 | 37.8944 | 77.7570 | 0.0737 |
| T=10,k=1,eta=1,stress=none | stationary_ignorant | 499 | 17.0676 | 37.8618 | 77.8503 | 0.0740 |
| T=10,k=1,eta=1,stress=none | ml | 498 | 0.5339 | 5.1154 | 9.7852 | 0.2708 |
| T=10,k=1,eta=1,stress=none | reml | 487 | 0.4709 | 4.5318 | 8.8243 | 0.3079 |
| T=10,k=2,eta=1,stress=none | conditional | 500 | 5.0906 | 15.3558 | 19.0059 | 0.0762 |
| T=10,k=2,eta=1,stress=none | stationary_ignorant | 498 | 5.1323 | 15.4515 | 19.1044 | 0.0762 |
| T=10,k=2,eta=1,stress=none | ml | 500 | 0.4701 | 5.6100 | 11.2828 | 0.3090 |
| T=10,k=2,eta=1,stress=none | reml | 492 | 0.5148 | 5.0095 | 10.3604 | 0.3424 |
| T=10,k=4,eta=1,stress=none | conditional | 499 | 1.6466 | 8.8071 | 10.3829 | 0.1310 |
| T=10,k=4,eta=1,stress=none | stationary_ignorant | 495 | 1.6289 | 8.8254 | 10.4786 | 0.1357 |
| T=10,k=4,eta=1,stress=none | ml | 499 | 0.4601 | 5.7181 | 7.6841 | 0.2312 |
| T=10,k=4,eta=1,stress=none | reml | 487 | 0.4529 | 4.9993 | 18.8050 | 0.4822 |
| T=10,k=1,eta=2,stress=none | conditional | 487 | 34.9709 | 84.4420 | 164.7011 | 0.0749 |
| T=10,k=1,eta=2,stress=none | stationary_ignorant | 484 | 34.8419 | 82.1412 | 162.3706 | 0.0706 |
| T=10,k=1,eta=2,stress=none | ml | 499 | 0.8607 | 45.1982 | 63.5164 | 0.4093 |
| T=10,k=1,eta=2,stress=none | reml | 477 | 0.7010 | 36.0188 | 61.0269 | 0.5504 |
| T=10,k=2,eta=2,stress=none | conditional | 477 | 12.8223 | 40.4514 | 55.5800 | 0.0936 |
| T=10,k=2,eta=2,stress=none | stationary_ignorant | 476 | 12.8427 | 40.5358 | 55.4521 | 0.0937 |
| T=10,k=2,eta=2,stress=none | ml | 488 | 0.6309 | 18.9520 | 36.2762 | 0.4709 |
| T=10,k=2,eta=2,stress=none | reml | 467 | 0.6571 | 17.1051 | 27.3841 | 0.4561 |
| T=10,k=4,eta=2,stress=none | conditional | 480 | 4.5156 | 17.4692 | 23.3356 | 0.0995 |
| T=10,k=4,eta=2,stress=none | stationary_ignorant | 479 | 4.5433 | 17.5747 | 24.0407 | 0.1025 |
| T=10,k=4,eta=2,stress=none | ml | 475 | 0.6425 | 10.0573 | 20.3846 | 0.3369 |
| T=10,k=4,eta=2,stress=none | reml | 451 | 0.6323 | 8.4607 | 12.6782 | 0.3129 |
| T=2,k=1,eta=0.5,stress=none | conditional | 500 | 54.8948 | 157.1958 | 227.0045 | 0.0771 |
| T=2,k=1,eta=0.5,stress=none | stationary_ignorant | 494 | 55.2610 | 157.5610 | 227.4908 | 0.0779 |
| T=2,k=1,eta=0.5,stress=none | ml | 497 | 1.9604 | 20.2335 | 40.5453 | 0.3307 |
| T=2,k=1,eta=0.5,stress=none | reml | 418 | 1.0985 | 18.0337 | 37.3475 | 0.4222 |
| T=30,k=1,eta=0.5,stress=none | conditional | 500 | 1.6929 | 4.2765 | 5.7048 | 0.0573 |
| T=30,k=1,eta=0.5,stress=none | stationary_ignorant | 498 | 1.6883 | 4.2446 | 5.6980 | 0.0558 |
| T=30,k=1,eta=0.5,stress=none | ml | 497 | 0.2791 | 2.1679 | 3.4224 | 0.2078 |
| T=30,k=1,eta=0.5,stress=none | reml | 497 | 0.2677 | 1.9708 | 3.2375 | 0.2168 |
| T=10,k=1,eta=0,stress=none | conditional | 500 | 0.3802 | 2.2414 | 4.3172 | 0.1442 |
| T=10,k=1,eta=0,stress=none | stationary_ignorant | 488 | 0.3606 | 2.1968 | 4.2317 | 0.1536 |
| T=10,k=1,eta=0,stress=none | ml | 497 | 0.3372 | 2.1347 | 4.2317 | 0.1667 |
| T=10,k=1,eta=0,stress=none | reml | 485 | 0.3261 | 1.9206 | 4.0021 | 0.1831 |
| T=10,k=1,eta=0.5,stress=colored | conditional | 500 | 2.0910 | 5.6739 | 6.6803 | 0.0566 |
| T=10,k=1,eta=0.5,stress=colored | stationary_ignorant | 490 | 2.0849 | 5.7150 | 6.7612 | 0.0594 |
| T=10,k=1,eta=0.5,stress=colored | ml | 496 | 1.4655 | 5.4202 | 6.7612 | 0.0910 |
| T=10,k=1,eta=0.5,stress=colored | reml | 499 | 1.2214 | 5.1118 | 6.5179 | 0.1024 |
| T=10,k=1,eta=0.5,stress=fixed5 | conditional | 500 | 2.7875 | 6.3280 | 9.5759 | 0.0623 |
| T=10,k=1,eta=0.5,stress=fixed5 | stationary_ignorant | 492 | 1.8386 | 5.0417 | 8.4620 | 0.0842 |
| T=10,k=1,eta=0.5,stress=fixed5 | ml | 495 | 0.6884 | 0.8864 | 0.9052 | 0.0178 |
| T=10,k=1,eta=0.5,stress=fixed5 | reml | 177 | 0.9438 | 0.9985 | 0.9988 | 0.0350 |
