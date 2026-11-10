# Hadamard Maximal Determinant Catalog (`1 <= n <= 120`)

Unified, offline-verifiable reference catalog of maximal determinants of $n \times n$ sign matrices $X \in \{-1, +1\}^{n \times n}$ across all orders $1 \le n \le 120$, consolidating historical constructions and optimality results from William P. Orrick's archived `MaxDet` database (`indiana.edu/~maxdet`, 2002-2012), Tomas Rokicki's circulant-pair tables (2005), Hiroki Tamura's constructions (2005-2011), Richard P. Brent's optimality proof for $n = 22$ (2021), and the two-circulant reinforcement-learning records of Butbaia, Gukov et al. (`arXiv:2608.22518v3`, Oct 2026).

### Normalization & Quick Verification

Every $n \times n$ sign matrix $X \in \{-1, +1\}^{n \times n}$ has $|\det X|$ divisible by $2^{n-1}$. Following `OEIS A003432` and Orrick's catalog, we report the exact normalized integer determinant

$$D_{\text{norm}}(n) = \frac{|\det X|}{2^{n-1}} \in \mathbb{Z}_{>0}, \qquad R(n) = \frac{|\det X|}{B_{\text{class}}(n)} \le 1,$$

which equals the maximal determinant of an $(n-1) \times (n-1)$ binary $\{0, 1\}$ matrix (`OEIS A003433` tracks raw $|\det X| = 2^{n-1} D_{\text{norm}}(n)$). Explicit `{-1, +1}` witness matrices are included in [`matrices/literature_witnesses/`](matrices/literature_witnesses/) for `49 / 50` open orders (`n = 89` had no explicit matrix archived in `MaxDet`).

```bash
python3 verify.py            # Verifies all 49 explicit open-order witnesses, 120 bounds, CSV/JSON parity & SHA256SUMS
python3 verify.py --corrupt  # Adversarial 1-entry mutation self-test (must print FAIL and exit 1)
```

- **Machine-Readable Tables:** [`catalog/hadamard_maxdet_n1_to_120.csv`](catalog/hadamard_maxdet_n1_to_120.csv) · [`catalog/hadamard_maxdet_n1_to_120.json`](catalog/hadamard_maxdet_n1_to_120.json)
- **Companion Open-Science Repository:** <https://github.com/mizopl/we-star>

---

## 1. Residue-Class Summary (`1 <= n <= 120`)

| Residue Class | Orders | Proven Optimal | Open | Upper Bound `B_class(n)` | Status Summary |
| :--- | ---: | ---: | ---: | :--- | :--- |
| `0 mod 4` | `30` | `30` | `0` | Hadamard (`n^(n/2)`) | All 30 multiples of 4 up to 120 achieve `R = 1` |
| `1 mod 4` | `30` | `11` | `19` | Barba (`sqrt(2n-1)*(n-1)^((n-1)/2)`) | 11 proven optimal (`1,5,9,13,17,21,25,37,41,61,113`); 19 open |
| `2 mod 4` | `30` | `24` | `6` | Ehlich-Wojtas (`2(n-1)(n-2)^((n-2)/2)`) | 24 proven optimal (`2`, `22` [Brent 2021], 22 EW-tight); 6 open |
| `3 mod 4` | `30` | `5` | `25` | Ehlich (`s`-block Gram bound, `s in {3,5,6,7}`) | 5 proven optimal (`3,7,11,15,19`); 25 open |

---

## 2. Unified Leaderboard of All 50 Open Orders (`23 <= n <= 119`)

All 50 orders below have $R(n) < 1$ without a published optimality proof; `49` rows are verified directly from their linked `{-1, +1}` witness matrices by `python3 verify.py`.

| `n` | `mod 4` | Best Known `D_norm` | `R = |det|/B` | Current Record Holder | Latest Gain | Matrix |
| ---: | :---: | :--- | ---: | :--- | ---: | :---: |
| `23` | `3` | `662671875` | `0.931983` | Dowdeswell, Neubauer, Solomon & Tumer (2004) | - | [`n023.txt`](matrices/literature_witnesses/n023.txt) |
| `27` | `3` | `546 * 6^11` | `0.967268` | H. Tamura (2005) | - | [`n027.txt`](matrices/literature_witnesses/n027.txt) |
| `29` | `1` | `320 * 7^12` | `0.865001` | B. Solomon (2002) | - | [`n029.txt`](matrices/literature_witnesses/n029.txt) |
| `31` | `3` | `784 * 7^13` | `0.956070` | H. Tamura (2005) | - | [`n031.txt`](matrices/literature_witnesses/n031.txt) |
| `33` | `1` | `441 * 8^14` | `0.854677` | B. Solomon (2002) | - | [`n033.txt`](matrices/literature_witnesses/n033.txt) |
| `34` | `2` | `256 * 8^15` | `0.969697` | W. P. Orrick (2005) | - | [`n034.txt`](matrices/literature_witnesses/n034.txt) |
| `35` | `3` | `1064 * 8^15` | `0.937914` | H. Tamura (2011) | - | [`n035.txt`](matrices/literature_witnesses/n035.txt) |
| `39` | `3` | `1440 * 9^17` | `0.952661` | H. Tamura (2011) | - | [`n039.txt`](matrices/literature_witnesses/n039.txt) |
| `43` | `3` | `1890 * 10^19` | `0.966343` | H. Tamura (2011) | - | [`n043.txt`](matrices/literature_witnesses/n043.txt) |
| `45` | `1` | `89 * 11^21` | `0.857635` | Brent, Orrick, Osborn & Zimmermann (2009) | - | [`n045.txt`](matrices/literature_witnesses/n045.txt) |
| `47` | `3` | `3037500 * 11^18` | `0.923897` | Orrick, Solomon, Dowdeswell & Smith (2003) | - | [`n047.txt`](matrices/literature_witnesses/n047.txt) |
| `49` | `1` | `96 * 12^23` | `0.812277` | B. Solomon & W. P. Orrick (2003) | - | [`n049.txt`](matrices/literature_witnesses/n049.txt) |
| `51` | `3` | `1.7776e+28` | `0.877569` | **Butbaia, Gukov et al. (Caltech) (2026)** | `+3.1053%` | [`n051.txt`](matrices/literature_witnesses/n051.txt) |
| `53` | `1` | `105 * 13^25` | `0.788227` | Orrick, Solomon, Dowdeswell & Smith (2003) | - | [`n053.txt`](matrices/literature_witnesses/n053.txt) |
| `55` | `3` | `2.2572e+31` | `0.859320` | T. Rokicki; Kazmenko & Trofimov (2005) | - | [`n055.txt`](matrices/literature_witnesses/n055.txt) |
| `57` | `1` | `133 * 14^27` | `0.893685` | J.-C. Meyrignac & B. Felgenhauer (2005) | - | [`n057.txt`](matrices/literature_witnesses/n057.txt) |
| `58` | `2` | `6.5774e+33` | `0.934528` | T. Rokicki; Kazmenko & Trofimov (2005) | - | [`n058.txt`](matrices/literature_witnesses/n058.txt) |
| `59` | `3` | `4312 * 14^27` | `0.965304` | H. Tamura (2005) | - | [`n059.txt`](matrices/literature_witnesses/n059.txt) |
| `63` | `3` | `5115 * 15^29` | `0.966154` | H. Tamura (2011) | - | [`n063.txt`](matrices/literature_witnesses/n063.txt) |
| `65` | `1` | `148 * 16^31` | `0.814417` | W. P. Orrick & B. Solomon (2003) | - | [`n065.txt`](matrices/literature_witnesses/n065.txt) |
| `67` | `3` | `1.1400e+41` | `0.863532` | T. Rokicki (2005) | - | [`n067.txt`](matrices/literature_witnesses/n067.txt) |
| `69` | `1` | `155 * 17^33` | `0.778973` | Orrick, Solomon, Dowdeswell & Smith (2003) | - | [`n069.txt`](matrices/literature_witnesses/n069.txt) |
| `70` | `2` | `1156 * 17^33` | `0.985507` | W. P. Orrick & T. Rokicki (2005) | - | [`n070.txt`](matrices/literature_witnesses/n070.txt) |
| `71` | `3` | `31186944 * 17^30` | `0.880425` | T. Rokicki (2005) | - | [`n071.txt`](matrices/literature_witnesses/n071.txt) |
| `73` | `1` | `174 * 18^35` | `0.802773` | B. Solomon & W. P. Orrick (2003) | - | [`n073.txt`](matrices/literature_witnesses/n073.txt) |
| `75` | `3` | `6.1512e+47` | `0.861319` | T. Rokicki (2005) | - | [`n075.txt`](matrices/literature_witnesses/n075.txt) |
| `77` | `1` | `183 * 19^37` | `0.778667` | B. Solomon & W. P. Orrick (2003) | - | [`n077.txt`](matrices/literature_witnesses/n077.txt) |
| `78` | `2` | `2.9042e+50` | `0.963616` | T. Rokicki (2005) | - | [`n078.txt`](matrices/literature_witnesses/n078.txt) |
| `79` | `3` | `1.6796e+51` | `0.859050` | T. Rokicki (2005) | - | [`n079.txt`](matrices/literature_witnesses/n079.txt) |
| `81` | `1` | `196 * 20^39` | `0.772348` | N. Farmakis & S. Kounias (1987) | - | [`n081.txt`](matrices/literature_witnesses/n081.txt) |
| `83` | `3` | `5.0710e+54` | `0.856072` | T. Rokicki (2005) | - | [`n083.txt`](matrices/literature_witnesses/n083.txt) |
| `85` | `1` | `4509 * 21^40` | `0.786499` | B. Solomon & W. P. Orrick (2003) | - | [`n085.txt`](matrices/literature_witnesses/n085.txt) |
| `87` | `3` | `1.6851e+58` | `0.852678` | T. Rokicki (2005) | - | [`n087.txt`](matrices/literature_witnesses/n087.txt) |
| `89` | `1` | `220 * 22^43` | `0.751646` | N. Farmakis & S. Kounias (1987) | - | none* |
| `91` | `3` | `6.1448e+61` | `0.850062` | T. Rokicki (2005) | - | [`n091.txt`](matrices/literature_witnesses/n091.txt) |
| `93` | `1` | `239 * 23^45` | `0.763984` | W. P. Orrick & B. Solomon (2003) | - | [`n093.txt`](matrices/literature_witnesses/n093.txt) |
| `94` | `2` | `3.9312e+64` | `0.969522` | T. Rokicki (2005) | - | [`n094.txt`](matrices/literature_witnesses/n094.txt) |
| `95` | `3` | `2.4387e+65` | `0.844737` | T. Rokicki (2005) | - | [`n095.txt`](matrices/literature_witnesses/n095.txt) |
| `97` | `1` | `256 * 24^47` | `0.767803` | B. Solomon & W. P. Orrick (2003) | - | [`n097.txt`](matrices/literature_witnesses/n097.txt) |
| `99` | `3` | `114487974849 * 24^42` | `0.849602` | T. Rokicki (2005) | - | [`n099.txt`](matrices/literature_witnesses/n099.txt) |
| `101` | `1` | `7110 * 25^48` | `0.802401` | W. P. Orrick & B. Solomon (2003) | - | [`n101.txt`](matrices/literature_witnesses/n101.txt) |
| `103` | `3` | `4.9829e+72` | `0.843316` | **Butbaia, Gukov et al. (Caltech) (2026)** | `+0.0000%` | [`n103.txt`](matrices/literature_witnesses/n103.txt) |
| `105` | `1` | `3.8474e+74` | `0.702238` | T. Rokicki (2005) | - | [`n105.txt`](matrices/literature_witnesses/n105.txt) |
| `106` | `2` | `2704 * 26^51` | `0.990476` | W. P. Orrick (2005) | - | [`n106.txt`](matrices/literature_witnesses/n106.txt) |
| `107` | `3` | `2.5405e+76` | `0.844532` | **Butbaia, Gukov et al. (Caltech) (2026)** | `+0.4364%` | [`n107.txt`](matrices/literature_witnesses/n107.txt) |
| `109` | `1` | `1.8775e+78` | `0.648191` | T. Rokicki (2005) | - | [`n109.txt`](matrices/literature_witnesses/n109.txt) |
| `111` | `3` | `1.3978e+80` | `0.846933` | **Butbaia, Gukov et al. (Caltech) (2026)** | `+1.2629%` | [`n111.txt`](matrices/literature_witnesses/n111.txt) |
| `115` | `3` | `8.2608e+83` | `0.848842` | **Butbaia, Gukov et al. (Caltech) (2026)** | `+1.8288%` | [`n115.txt`](matrices/literature_witnesses/n115.txt) |
| `117` | `1` | `7.2315e+85` | `0.718566` | T. Rokicki (2005) | - | [`n117.txt`](matrices/literature_witnesses/n117.txt) |
| `119` | `3` | `5.2308e+87` | `0.850263` | **Butbaia, Gukov et al. (Caltech) (2026)** | `+2.1204%` | [`n119.txt`](matrices/literature_witnesses/n119.txt) |

<sub>* On `n = 89`, the bound `220 * 22^43` was derived by Farmakis & Kounias (1987) from an order-88 max-excess Hadamard construction; no explicit `89 x 89` matrix was on file in the `MaxDet` archive.</sub>

<details>
<summary><strong>Complete 120-Order Master Table (`1 <= n <= 120`, Proven Optimal + Open)</strong></summary>

| `n` | `mod 4` | Status | `D_norm` | `R = |det|/B` | Attribution | Ref | Matrix |
| ---: | :---: | :---: | :--- | ---: | :--- | :---: | :---: |
| `1` | `1` | `OPT` | `1` | `1.000000` | Sylvester (1867) | `Syl1867` | - |
| `2` | `2` | `OPT` | `1` | `1.000000` | Sylvester (1867) | `Syl1867` | - |
| `3` | `3` | `OPT` | `1` | `1.000000` | Sylvester; Williamson (1867) | `Syl1867` | - |
| `4` | `0` | `OPT` | `2` | `1.000000` | Sylvester (Hadamard) (1867) | `Syl1867` | - |
| `5` | `1` | `OPT` | `3` | `1.000000` | Mood (Barba bound) (1946) | `Moo1946` | - |
| `6` | `2` | `OPT` | `5` | `1.000000` | Ehlich (EW bound) (1964) | `Ehl1964a` | - |
| `7` | `3` | `OPT` | `9` | `0.981981` | Williamson (1946) | `Wil1946` | - |
| `8` | `0` | `OPT` | `32` | `1.000000` | Sylvester (Hadamard) (1867) | `Syl1867` | - |
| `9` | `1` | `OPT` | `56` | `0.848875` | Ehlich (1964) | `Ehl1964a` | - |
| `10` | `2` | `OPT` | `144` | `1.000000` | Ehlich (EW bound) (1964) | `Ehl1964a` | - |
| `11` | `3` | `OPT` | `320` | `0.944911` | Ehlich (1964) | `Ehl1964b` | - |
| `12` | `0` | `OPT` | `1458` | `1.000000` | Hadamard (1893) | `Had1893` | - |
| `13` | `1` | `OPT` | `3645` | `1.000000` | Raghavarao (Barba bound) (1959) | `Rag1959` | - |
| `14` | `2` | `OPT` | `9477` | `1.000000` | Ehlich (EW bound) (1964) | `Ehl1964a` | - |
| `15` | `3` | `OPT` | `25515` | `0.970725` | Ehlich & Zeller; Cohn (1962) | `EZ1962` | - |
| `16` | `0` | `OPT` | `131072` | `1.000000` | Sylvester (Hadamard) (1867) | `Syl1867` | - |
| `17` | `1` | `OPT` | `327680` | `0.870388` | Moyssiadis & Kounias (1982) | `MK1982` | - |
| `18` | `2` | `OPT` | `1114112` | `1.000000` | Ehlich (EW bound) (1964) | `Ehl1964a` | - |
| `19` | `3` | `OPT` | `3411968` | `0.975017` | Cohn (1989); opt. BOOZ (2011) (2011) | `BOOZ2011` | - |
| `20` | `0` | `OPT` | `19531250` | `1.000000` | Hadamard (1893) | `Had1893` | - |
| `21` | `1` | `OPT` | `56640625` | `0.905808` | Chadjipantelis, Kounias & Moyssiadis (1987) | `CKM1987` | - |
| `22` | `2` | `OPT` | `195312500` | `0.952381` | Ehlich (1964); opt. Brent (2021) (2021) | `Bre2021` | - |
| `23` | `3` | `OPEN` | `662671875` | `0.931983` | Dowdeswell, Neubauer, Solomon & Tumer (2004) | `DNST2004` | [`txt`](matrices/literature_witnesses/n023.txt) |
| `24` | `0` | `OPT` | `12 * 6^11` | `1.000000` | Paley / Williamson (Hadamard) (1933) | `Pal1933` | - |
| `25` | `1` | `OPT` | `42 * 6^11` | `1.000000` | Brouwer (Barba bound) (1983) | `Bro1983` | - |
| `26` | `2` | `OPT` | `150 * 6^11` | `1.000000` | Ehlich (EW bound) (1964) | `Ehl1964a` | - |
| `27` | `3` | `OPEN` | `546 * 6^11` | `0.967268` | H. Tamura (2005) | `Tam2005` | [`txt`](matrices/literature_witnesses/n027.txt) |
| `28` | `0` | `OPT` | `14 * 7^13` | `1.000000` | Paley / Williamson (Hadamard) (1933) | `Pal1933` | - |
| `29` | `1` | `OPEN` | `320 * 7^12` | `0.865001` | B. Solomon (2002) | `Sol2002` | [`txt`](matrices/literature_witnesses/n029.txt) |
| `30` | `2` | `OPT` | `203 * 7^13` | `1.000000` | Ehlich (EW bound) (1964) | `Ehl1964a` | - |
| `31` | `3` | `OPEN` | `784 * 7^13` | `0.956070` | H. Tamura (2005) | `Tam2005` | [`txt`](matrices/literature_witnesses/n031.txt) |
| `32` | `0` | `OPT` | `16 * 8^15` | `1.000000` | Sylvester (Hadamard) (1867) | `Syl1867` | - |
| `33` | `1` | `OPEN` | `441 * 8^14` | `0.854677` | B. Solomon (2002) | `Sol2002` | [`txt`](matrices/literature_witnesses/n033.txt) |
| `34` | `2` | `OPEN` | `256 * 8^15` | `0.969697` | W. P. Orrick (2005) | `Orr2005` | [`txt`](matrices/literature_witnesses/n034.txt) |
| `35` | `3` | `OPEN` | `1064 * 8^15` | `0.937914` | H. Tamura (2011) | `Tam2011` | [`txt`](matrices/literature_witnesses/n035.txt) |
| `36` | `0` | `OPT` | `18 * 9^17` | `1.000000` | Paley / Williamson (Hadamard) (1933) | `Pal1933` | - |
| `37` | `1` | `OPT` | `72 * 9^17` | `0.936329` | Brent, Orrick, Osborn & Zimmermann (2011) | `BOOZ2011` | - |
| `38` | `2` | `OPT` | `333 * 9^17` | `1.000000` | Ehlich (EW bound) (1964) | `Ehl1964a` | - |
| `39` | `3` | `OPEN` | `1440 * 9^17` | `0.952661` | H. Tamura (2011) | `Tam2011` | [`txt`](matrices/literature_witnesses/n039.txt) |
| `40` | `0` | `OPT` | `20 * 10^19` | `1.000000` | Paley / Williamson (Hadamard) (1933) | `Pal1933` | - |
| `41` | `1` | `OPT` | `90 * 10^19` | `1.000000` | Bridges, Hall & Hayden; van Trung (1981) | `BHH1981` | - |
| `42` | `2` | `OPT` | `410 * 10^19` | `1.000000` | Yang (EW bound) (1968) | `Yan1968` | - |
| `43` | `3` | `OPEN` | `1890 * 10^19` | `0.966343` | H. Tamura (2011) | `Tam2011` | [`txt`](matrices/literature_witnesses/n043.txt) |
| `44` | `0` | `OPT` | `22 * 11^21` | `1.000000` | Paley / Williamson (Hadamard) (1933) | `Pal1933` | - |
| `45` | `1` | `OPEN` | `89 * 11^21` | `0.857635` | Brent, Orrick, Osborn & Zimmermann (2009) | `BOOZ2011` | [`txt`](matrices/literature_witnesses/n045.txt) |
| `46` | `2` | `OPT` | `495 * 11^21` | `1.000000` | Yang (EW bound) (1968) | `Yan1968` | - |
| `47` | `3` | `OPEN` | `3037500 * 11^18` | `0.923897` | Orrick, Solomon, Dowdeswell & Smith (2003) | `OSDS2003` | [`txt`](matrices/literature_witnesses/n047.txt) |
| `48` | `0` | `OPT` | `24 * 12^23` | `1.000000` | Paley / Williamson (Hadamard) (1933) | `Pal1933` | - |
| `49` | `1` | `OPEN` | `96 * 12^23` | `0.812277` | B. Solomon & W. P. Orrick (2003) | `OS2007` | [`txt`](matrices/literature_witnesses/n049.txt) |
| `50` | `2` | `OPT` | `588 * 12^23` | `1.000000` | Yang (EW bound) (1968) | `Yan1968` | - |
| `51` | `3` | `OPEN` | `1.7776e+28` | `0.877569` | Butbaia, Gukov et al. (Caltech) (2026) | `BG2026` | [`txt`](matrices/literature_witnesses/n051.txt) |
| `52` | `0` | `OPT` | `26 * 13^25` | `1.000000` | Paley / Williamson (Hadamard) (1933) | `Pal1933` | - |
| `53` | `1` | `OPEN` | `105 * 13^25` | `0.788227` | Orrick, Solomon, Dowdeswell & Smith (2003) | `OSDS2003` | [`txt`](matrices/literature_witnesses/n053.txt) |
| `54` | `2` | `OPT` | `689 * 13^25` | `1.000000` | Yang (EW bound) (1968) | `Yan1968` | - |
| `55` | `3` | `OPEN` | `2.2572e+31` | `0.859320` | T. Rokicki; Kazmenko & Trofimov (2005) | `Rok2005` | [`txt`](matrices/literature_witnesses/n055.txt) |
| `56` | `0` | `OPT` | `28 * 14^27` | `1.000000` | Paley / Williamson (Hadamard) (1933) | `Pal1933` | - |
| `57` | `1` | `OPEN` | `133 * 14^27` | `0.893685` | J.-C. Meyrignac & B. Felgenhauer (2005) | `Mey2005` | [`txt`](matrices/literature_witnesses/n057.txt) |
| `58` | `2` | `OPEN` | `6.5774e+33` | `0.934528` | T. Rokicki; Kazmenko & Trofimov (2005) | `Rok2005` | [`txt`](matrices/literature_witnesses/n058.txt) |
| `59` | `3` | `OPEN` | `4312 * 14^27` | `0.965304` | H. Tamura (2005) | `Tam2005` | [`txt`](matrices/literature_witnesses/n059.txt) |
| `60` | `0` | `OPT` | `30 * 15^29` | `1.000000` | Paley / Williamson (Hadamard) (1933) | `Pal1933` | - |
| `61` | `1` | `OPT` | `165 * 15^29` | `1.000000` | Brouwer (Barba bound) (1983) | `Bro1983` | - |
| `62` | `2` | `OPT` | `915 * 15^29` | `1.000000` | Yang (EW bound) (1968) | `Yan1968` | - |
| `63` | `3` | `OPEN` | `5115 * 15^29` | `0.966154` | H. Tamura (2011) | `Tam2011` | [`txt`](matrices/literature_witnesses/n063.txt) |
| `64` | `0` | `OPT` | `32 * 16^31` | `1.000000` | Sylvester (Hadamard) (1867) | `Syl1867` | - |
| `65` | `1` | `OPEN` | `148 * 16^31` | `0.814417` | W. P. Orrick & B. Solomon (2003) | `OS2007` | [`txt`](matrices/literature_witnesses/n065.txt) |
| `66` | `2` | `OPT` | `1040 * 16^31` | `1.000000` | Yang (EW bound) (1968) | `Yan1968` | - |
| `67` | `3` | `OPEN` | `1.1400e+41` | `0.863532` | T. Rokicki (2005) | `Rok2005` | [`txt`](matrices/literature_witnesses/n067.txt) |
| `68` | `0` | `OPT` | `34 * 17^33` | `1.000000` | Paley / Williamson (Hadamard) (1933) | `Pal1933` | - |
| `69` | `1` | `OPEN` | `155 * 17^33` | `0.778973` | Orrick, Solomon, Dowdeswell & Smith (2003) | `OSDS2003` | [`txt`](matrices/literature_witnesses/n069.txt) |
| `70` | `2` | `OPEN` | `1156 * 17^33` | `0.985507` | W. P. Orrick & T. Rokicki (2005) | `Orr2005` | [`txt`](matrices/literature_witnesses/n070.txt) |
| `71` | `3` | `OPEN` | `31186944 * 17^30` | `0.880425` | T. Rokicki (2005) | `Rok2005` | [`txt`](matrices/literature_witnesses/n071.txt) |
| `72` | `0` | `OPT` | `36 * 18^35` | `1.000000` | Paley / Williamson (Hadamard) (1933) | `Pal1933` | - |
| `73` | `1` | `OPEN` | `174 * 18^35` | `0.802773` | B. Solomon & W. P. Orrick (2003) | `OS2007` | [`txt`](matrices/literature_witnesses/n073.txt) |
| `74` | `2` | `OPT` | `1314 * 18^35` | `1.000000` | Cohn (EW bound) (1989) | `Coh1989` | - |
| `75` | `3` | `OPEN` | `6.1512e+47` | `0.861319` | T. Rokicki (2005) | `Rok2005` | [`txt`](matrices/literature_witnesses/n075.txt) |
| `76` | `0` | `OPT` | `38 * 19^37` | `1.000000` | Paley / Williamson (Hadamard) (1933) | `Pal1933` | - |
| `77` | `1` | `OPEN` | `183 * 19^37` | `0.778667` | B. Solomon & W. P. Orrick (2003) | `OS2007` | [`txt`](matrices/literature_witnesses/n077.txt) |
| `78` | `2` | `OPEN` | `2.9042e+50` | `0.963616` | T. Rokicki (2005) | `Rok2005` | [`txt`](matrices/literature_witnesses/n078.txt) |
| `79` | `3` | `OPEN` | `1.6796e+51` | `0.859050` | T. Rokicki (2005) | `Rok2005` | [`txt`](matrices/literature_witnesses/n079.txt) |
| `80` | `0` | `OPT` | `40 * 20^39` | `1.000000` | Paley / Williamson (Hadamard) (1933) | `Pal1933` | - |
| `81` | `1` | `OPEN` | `196 * 20^39` | `0.772348` | N. Farmakis & S. Kounias (1987) | `FK1987` | [`txt`](matrices/literature_witnesses/n081.txt) |
| `82` | `2` | `OPT` | `1620 * 20^39` | `1.000000` | Cohn (EW bound) (1989) | `Coh1989` | - |
| `83` | `3` | `OPEN` | `5.0710e+54` | `0.856072` | T. Rokicki (2005) | `Rok2005` | [`txt`](matrices/literature_witnesses/n083.txt) |
| `84` | `0` | `OPT` | `42 * 21^41` | `1.000000` | Paley / Williamson (Hadamard) (1933) | `Pal1933` | - |
| `85` | `1` | `OPEN` | `4509 * 21^40` | `0.786499` | B. Solomon & W. P. Orrick (2003) | `OS2007` | [`txt`](matrices/literature_witnesses/n085.txt) |
| `86` | `2` | `OPT` | `1785 * 21^41` | `1.000000` | Chadjipantelis & Kounias (EW bound) (1985) | `CK1985` | - |
| `87` | `3` | `OPEN` | `1.6851e+58` | `0.852678` | T. Rokicki (2005) | `Rok2005` | [`txt`](matrices/literature_witnesses/n087.txt) |
| `88` | `0` | `OPT` | `44 * 22^43` | `1.000000` | Paley / Williamson (Hadamard) (1933) | `Pal1933` | - |
| `89` | `1` | `OPEN` | `220 * 22^43` | `0.751646` | N. Farmakis & S. Kounias (1987) | `FK1987` | - |
| `90` | `2` | `OPT` | `1958 * 22^43` | `1.000000` | Cohn (EW bound) (1989) | `Coh1989` | - |
| `91` | `3` | `OPEN` | `6.1448e+61` | `0.850062` | T. Rokicki (2005) | `Rok2005` | [`txt`](matrices/literature_witnesses/n091.txt) |
| `92` | `0` | `OPT` | `46 * 23^45` | `1.000000` | Baumert, Golomb & Hall (1962) | `BGH1962` | - |
| `93` | `1` | `OPEN` | `239 * 23^45` | `0.763984` | W. P. Orrick & B. Solomon (2003) | `OS2007` | [`txt`](matrices/literature_witnesses/n093.txt) |
| `94` | `2` | `OPEN` | `3.9312e+64` | `0.969522` | T. Rokicki (2005) | `Rok2005` | [`txt`](matrices/literature_witnesses/n094.txt) |
| `95` | `3` | `OPEN` | `2.4387e+65` | `0.844737` | T. Rokicki (2005) | `Rok2005` | [`txt`](matrices/literature_witnesses/n095.txt) |
| `96` | `0` | `OPT` | `48 * 24^47` | `1.000000` | Paley / Williamson (Hadamard) (1933) | `Pal1933` | - |
| `97` | `1` | `OPEN` | `256 * 24^47` | `0.767803` | B. Solomon & W. P. Orrick (2003) | `OS2007` | [`txt`](matrices/literature_witnesses/n097.txt) |
| `98` | `2` | `OPT` | `2328 * 24^47` | `1.000000` | Cohn (EW bound) (1989) | `Coh1989` | - |
| `99` | `3` | `OPEN` | `114487974849 * 24^42` | `0.849602` | T. Rokicki (2005) | `Rok2005` | [`txt`](matrices/literature_witnesses/n099.txt) |
| `100` | `0` | `OPT` | `50 * 25^49` | `1.000000` | Paley / Williamson (Hadamard) (1933) | `Pal1933` | - |
| `101` | `1` | `OPEN` | `7110 * 25^48` | `0.802401` | W. P. Orrick & B. Solomon (2003) | `OS2007` | [`txt`](matrices/literature_witnesses/n101.txt) |
| `102` | `2` | `OPT` | `2525 * 25^49` | `1.000000` | Cohn (EW bound) (1989) | `Coh1989` | - |
| `103` | `3` | `OPEN` | `4.9829e+72` | `0.843316` | Butbaia, Gukov et al. (Caltech) (2026) | `BG2026` | [`txt`](matrices/literature_witnesses/n103.txt) |
| `104` | `0` | `OPT` | `52 * 26^51` | `1.000000` | Paley / Williamson (Hadamard) (1933) | `Pal1933` | - |
| `105` | `1` | `OPEN` | `3.8474e+74` | `0.702238` | T. Rokicki (2005) | `Rok2005` | [`txt`](matrices/literature_witnesses/n105.txt) |
| `106` | `2` | `OPEN` | `2704 * 26^51` | `0.990476` | W. P. Orrick (2005) | `Orr2005` | [`txt`](matrices/literature_witnesses/n106.txt) |
| `107` | `3` | `OPEN` | `2.5405e+76` | `0.844532` | Butbaia, Gukov et al. (Caltech) (2026) | `BG2026` | [`txt`](matrices/literature_witnesses/n107.txt) |
| `108` | `0` | `OPT` | `54 * 27^53` | `1.000000` | Paley / Williamson (Hadamard) (1933) | `Pal1933` | - |
| `109` | `1` | `OPEN` | `1.8775e+78` | `0.648191` | T. Rokicki (2005) | `Rok2005` | [`txt`](matrices/literature_witnesses/n109.txt) |
| `110` | `2` | `OPT` | `2943 * 27^53` | `1.000000` | Fletcher & Seberry (EW bound) (2001) | `FS2001` | - |
| `111` | `3` | `OPEN` | `1.3978e+80` | `0.846933` | Butbaia, Gukov et al. (Caltech) (2026) | `BG2026` | [`txt`](matrices/literature_witnesses/n111.txt) |
| `112` | `0` | `OPT` | `56 * 28^55` | `1.000000` | Paley / Williamson (Hadamard) (1933) | `Pal1933` | - |
| `113` | `1` | `OPT` | `420 * 28^55` | `1.000000` | Brouwer (Barba bound) (1983) | `Bro1983` | - |
| `114` | `2` | `OPT` | `3164 * 28^55` | `1.000000` | Koukouvinos, Kounias & Seberry (EW) (1989) | `KKS1989` | - |
| `115` | `3` | `OPEN` | `8.2608e+83` | `0.848842` | Butbaia, Gukov et al. (Caltech) (2026) | `BG2026` | [`txt`](matrices/literature_witnesses/n115.txt) |
| `116` | `0` | `OPT` | `58 * 29^57` | `1.000000` | Baumert (1966) | `Bau1966` | - |
| `117` | `1` | `OPEN` | `7.2315e+85` | `0.718566` | T. Rokicki (2005) | `Rok2005` | [`txt`](matrices/literature_witnesses/n117.txt) |
| `118` | `2` | `OPT` | `3393 * 29^57` | `1.000000` | Fletcher, Koukouvinos & Seberry (2002) | `FKS2002` | - |
| `119` | `3` | `OPEN` | `5.2308e+87` | `0.850263` | Butbaia, Gukov et al. (Caltech) (2026) | `BG2026` | [`txt`](matrices/literature_witnesses/n119.txt) |
| `120` | `0` | `OPT` | `60 * 30^59` | `1.000000` | Paley / Williamson (Hadamard) (1933) | `Pal1933` | - |
</details>

---

## 3. Attribution, Provenance & References

This catalog stands on over a century of combinatorial design theory and builds directly on the archival curation of **William P. Orrick** and **Bruce Solomon** (`MaxDet`, Indiana University, 2002-2012), the circulant-pair searches of **Tomas Rokicki** (2005) and **Hiroki Tamura** (2005-2011), and the two-circulant RL study of **Giorgi Butbaia, Sergei Gukov et al.** (`Math-AI-Caltech/hadamard-maxdet`, `arXiv:2608.22518v3`, 2026). Corrections, missing historical references, or new verified constructions are welcome via issues or pull requests.

<details>
<summary><strong>Full Bibliographic Key (`[Syl1867]` - `[BG2026]`)</strong></summary>

- **`[Syl1867]`** J. J. Sylvester (1867), Thoughts on inverse orthogonal matrices, Phil. Mag. 34, 461-475.
- **`[Had1893]`** J. Hadamard (1893), Resolution d'une question relative aux determinants, Bull. Sci. Math. 17, 240-246.
- **`[Pal1933]`** R. E. A. C. Paley (1933), On orthogonal matrices, J. Math. Phys. 12, 311-320.
- **`[Bar1933]`** G. Barba (1933), Intorno al teorema di Hadamard sui determinanti a valore massimo, Giorn. Mat. Battaglini 71, 70-86.
- **`[Wil1944]`** J. Williamson (1944), Hadamard's determinant theorem and the sum of four squares, Duke Math. J. 11, 65-81.
- **`[Wil1946]`** J. Williamson (1946), Determinants whose elements are 0 and 1, Amer. Math. Monthly 53, 427-434.
- **`[Moo1946]`** A. M. Mood (1946), On Hotelling's weighing problem, Ann. Math. Statist. 17, 432-446.
- **`[Rag1959]`** D. Raghavarao (1959), Some optimum weighing designs, Ann. Math. Statist. 30, 295-303.
- **`[EZ1962]`** H. Ehlich and K. Zeller (1962), Binare Matrizen, Z. Angew. Math. Mech. 42, T20-T21.
- **`[BGH1962]`** L. Baumert, S. W. Golomb, and M. Hall Jr. (1962), Discovery of an Hadamard matrix of order 92, Bull. Amer. Math. Soc. 68, 237-238.
- **`[Ehl1964a]`** H. Ehlich (1964), Determinantenabschatzungen fur binare Matrizen, Math. Z. 83, 123-132.
- **`[Ehl1964b]`** H. Ehlich (1964), Determinantenabschatzungen fur binare Matrizen mit n = 3 mod 4, Math. Z. 84, 438-447.
- **`[Woj1964]`** M. Wojtas (1964), On Hadamard's inequality for the determinants of order non-divisible by 4, Colloq. Math. 12, 73-83.
- **`[Bau1966]`** L. D. Baumert (1966), Hadamard matrices of orders 116 and 232, Bull. Amer. Math. Soc. 72, 237.
- **`[Yan1968]`** C. H. Yang (1966-1969), Some designs for maximal (+1,-1)-determinant of order n = 2 mod 4, Math. Comp. 20/22/23.
- **`[BHH1981]`** W. G. Bridges, M. Hall Jr., and M. L. Hayden (1981), Codes and designs, J. Combin. Theory Ser. A 31, 155-174; T. van Trung (1982).
- **`[MK1982]`** C. Moyssiadis and S. Kounias (1982), The exact D-optimal first order saturated design with 17 observations, J. Statist. Plann. Inference 7, 13-27.
- **`[Bro1983]`** A. E. Brouwer (1983), An infinite series of symmetric designs, Math. Centrum Amsterdam Report ZW 202/83.
- **`[CK1985]`** T. Chadjipantelis and S. Kounias (1985), Supplementary difference sets and D-optimal designs for n = 2 mod 4, Discrete Math. 57, 211-216.
- **`[CKM1987]`** T. Chadjipantelis, S. Kounias, and C. Moyssiadis (1987), The maximum determinant of 21x21 (+1,-1)-matrices and D-optimal designs, J. Statist. Plann. Inference 16, 167-178.
- **`[FK1987]`** N. Farmakis and S. Kounias (1987), The excess of Hadamard matrices and optimal designs, Discrete Math. 67, 165-176.
- **`[Coh1989]`** J. H. E. Cohn (1989), On determinants with elements +-1, II, Bull. London Math. Soc. 21, 36-42.
- **`[KKS1989]`** C. Koukouvinos, S. Kounias, and J. Seberry (1989), Further Hadamard matrices with maximal excess and new D-optimal designs for n = 2 mod 4, Utilitas Math. 36, 135-150.
- **`[FS2001]`** R. J. Fletcher and J. Seberry (2001), New D-optimal designs of order 110, Australas. J. Combin. 23, 49-52.
- **`[FKS2002]`** R. J. Fletcher, C. Koukouvinos, and J. Seberry (2002), New skew-Hadamard matrices of order 4*59 and D-optimal designs of order 2*59, Discrete Math. 286 (2004), 251-253.
- **`[Sol2002]`** B. Solomon (2002), Large-determinant sign matrices of orders 29 and 33, reported in W. P. Orrick's MaxDet catalog.
- **`[OSDS2003]`** W. P. Orrick, B. Solomon, R. Dowdeswell, and W. D. Smith (2003), New lower bounds for the maximal determinant problem, arXiv:math/0304410.
- **`[OS2007]`** W. P. Orrick and B. Solomon (2003; publ. 2007), Large-determinant sign matrices of order 4k+1, Discrete Math. 307, 226-236 (arXiv:math/0311292).
- **`[DNST2004]`** A. Dowdeswell, M. G. Neubauer, B. Solomon, and K. Tumer (2004), Lower bound for order 23 in W. P. Orrick's MaxDet catalog.
- **`[Orr2005]`** W. P. Orrick (2005), Switching operations for Hadamard matrices and D-optimal designs of orders 34, 70, 106, SIAM J. Discrete Math. 22 (2008), 31-50.
- **`[Rok2005]`** T. Rokicki (Apr 2005), New determinant records from circulant pairs (orders 51-119), archived in W. P. Orrick's MaxDet catalog.
- **`[Mey2005]`** J.-C. Meyrignac and B. Felgenhauer (Apr 2005), Order 57 sign matrix record, archived in W. P. Orrick's MaxDet catalog.
- **`[Tam2005]`** H. Tamura (Aug 2005), Maximal determinant records for orders 27, 31, and 59, archived in W. P. Orrick's MaxDet catalog.
- **`[BOOZ2011]`** R. P. Brent, W. P. Orrick, J.-a. H. Osborn, and P. Zimmermann (2009-2011), Maximal determinants and saturated D-optimal designs of orders 19 and 37 (and order 45 construction), arXiv:1112.4160.
- **`[Tam2011]`** H. Tamura (Jan 2011), Maximal determinant records for orders 35, 39, 43, and 63, archived in W. P. Orrick's MaxDet catalog.
- **`[DK2015]`** D. Z. Djokovic and I. S. Kotsireas (2015), D-optimal matrices of orders 118, 138, 150, 154 and 174, Algebraic Design Theory and Hadamard Matrices, Springer, 61-72.
- **`[Bre2021]`** R. P. Brent (2021), Optimality confirmation for binary order 21 (sign order 22, a(21) = 195312500), OEIS A003432.
- **`[BG2026]`** G. Butbaia, J. Tan, P. Vipulanandan, X. Huang, T. Saunders-A'Court, L. Fagan, D. Passaro, M. Tarquini, and S. Gukov (Oct 2026), Reinforcement learning for the Hadamard maximal determinant problem, arXiv:2608.22518v3 (Math-AI-Caltech/hadamard-maxdet).
</details>

---

<sub>Curated and verified by **we\*** (`@mizopl`) · Released under the [MIT License](LICENSE).</sub>
