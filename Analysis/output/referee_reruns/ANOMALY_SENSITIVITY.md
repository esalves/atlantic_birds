# Anomaly-definition sensitivity of the thermal-allometry estimate

Generated: 2026-09-18 16:44

The reported estimate (beta = +0.0212, p = 0.030, n = 7,577) has no generating
script in the repository. Reconstructions from `passer90_climate.rds` give:

```
 anomaly_definition    n    estimate          se         z           p
         A_base_all 7224 0.017455115 0.010040998 1.7383845 0.082143087
      B_base_shared 7224 0.019974029 0.010092030 1.9791884 0.047794798
        C_minus_amt 6945 0.015008813 0.005206365 2.8827814 0.003941809
           D_scaled 7224 0.002039493 0.010974389 0.1858411 0.852569363
 reported_estimate reported_p reported_n
            0.0212       0.03       7577
            0.0212       0.03       7577
            0.0212       0.03       7577
            0.0212       0.03       7577
```

No construction reaches n = 7,577: 353 of the 7,577 shared records have no
`rec_tmean`. The sign is stable across definitions but significance is not.

