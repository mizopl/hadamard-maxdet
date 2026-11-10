# Hadamard Maximal Determinant Catalog (1 <= n <= 120)

Unified, offline-verifiable reference catalog of maximal determinants of $n \times n$ sign matrices $X \in \{-1, +1\}^{n \times n}$ across all orders $1 \le n \le 120$, consolidating historical constructions and optimality proofs from William P. Orrick's archived [`MaxDet` database](http://wayback.archive-it.org/219/20190316001107/http://www.indiana.edu/~maxdet/) (2002–2012), Tomas Rokicki's circulant-pair tables (2005), Hiroki Tamura's constructions (2005–2011), Richard P. Brent's optimality proof for $n = 22$ ([OEIS A003432](https://oeis.org/A003432), 2021), and the two-circulant reinforcement-learning records of Butbaia, Gukov et al. ([arXiv:2608.22518](https://arxiv.org/abs/2608.22518), [`Math-AI-Caltech/hadamard-maxdet`](https://github.com/Math-AI-Caltech/hadamard-maxdet), Oct 2026).

### Normalization & Quick Verification

Every $n \times n$ sign matrix $X \in \{-1, +1\}^{n \times n}$ has $|\det X|$ divisible by $2^{n-1}$. Following [OEIS A003432](https://oeis.org/A003432) and Orrick's catalog, we report the exact normalized integer determinant

$$D_{\text{norm}}(n) = \frac{|\det X|}{2^{n-1}} \in \mathbb{Z}_{>0}, \qquad R(n) = \frac{|\det X|}{B_{\text{class}}(n)} \le 1,$$

which equals the maximal determinant of an $(n-1) \times (n-1)$ binary $\{0, 1\}$ matrix ([OEIS A003433](https://oeis.org/A003433) tracks raw $|\det X| = 2^{n-1} D_{\text{norm}}(n)$). Explicit `{-1, +1}` witness matrices are included in [`matrices/literature_witnesses/`](matrices/literature_witnesses/) for **49 / 50** open orders (`n = 89` had no explicit matrix archived in `MaxDet`).

```bash
python3 verify.py                         # Verifies all 49 open witnesses, 120 bounds, CSV/JSON parity & SHA256SUMS
python3 verify.py --check my_matrix.txt   # Evaluates any external +/-1 or 0/1 candidate matrix against catalog records
python3 verify.py --corrupt               # Adversarial 1-entry mutation self-test (must print FAIL and exit 1)
```

- **Machine-Readable Tables (full unrounded integers):** [`catalog/hadamard_maxdet_n1_to_120.csv`](catalog/hadamard_maxdet_n1_to_120.csv) · [`catalog/hadamard_maxdet_n1_to_120.json`](catalog/hadamard_maxdet_n1_to_120.json)
- **For AI Agents & Crawlers:** [`AGENTS.md`](AGENTS.md) · [`llms.txt`](llms.txt)
- **Companion Open-Science Repository:** <https://github.com/mizopl/we-star>

---

## 1. Residue-Class Summary (1 <= n <= 120)

| Residue Class | Orders | Optimal | Open | Upper Bound $B_{\text{class}}(n)$ | Status Summary |
| :--- | ---: | ---: | ---: | :--- | :--- |
| **0 mod 4** | 30 | 30 | 0 | Hadamard: $n^{n/2}$ | All 30 multiples of 4 up to 120 achieve $R = 1$ |
| **1 mod 4** | 30 | 11 | 19 | Barba: $\sqrt{2n-1}\,(n-1)^{(n-1)/2}$ | 11 optimal (1, 5, 9, 13, 17, 21, 25, 37, 41, 61, 113); 19 open |
| **2 mod 4** | 30 | 24 | 6 | Ehlich–Wojtas: $2(n-1)(n-2)^{(n-2)/2}$ | 24 optimal (2, 22 [Brent 2021], 22 EW-tight); 6 open |
| **3 mod 4** | 30 | 5 | 25 | Ehlich $s$-block Gram ($s \in \{3,5,6,7\}$) | 5 optimal (3, 7, 11, 15, 19); 25 open |

---

## 2. Best Known Lower Bounds on All 50 Open Orders (23 <= n <= 119)

To keep $R(n) = |\det X| / B_{\text{class}}(n)$ comparable down each column, the 50 open orders ($R(n) < 1$, no published optimality proof) are split into their three residue classes below (`49` rows verified directly from `{-1, +1}` matrices by [`verify.py`](verify.py); exact unrounded integers in [`hadamard_maxdet_n1_to_120.csv`](catalog/hadamard_maxdet_n1_to_120.csv)):

### 2.1 Orders n ≡ 1 (mod 4) — Barba Regime (19 Open Orders)

*Upper bound: $B_1(n) = \sqrt{2n - 1}\, (n - 1)^{(n - 1)/2}$. Factored entries take the form $c \cdot k^e$ with $k = (n - 1)/4$.*

| n | D<sub>norm</sub> | Factored Form | R = \|det\|/B | Record Holder | Year | Witness |
| ---: | ---: | :--- | ---: | :--- | :---: | :---: |
| **29** | 4.4292 × 10<sup>12</sup> | 320 · 7<sup>12</sup> | 0.8650 | B. Solomon | 2002 | [n029.txt](matrices/literature_witnesses/n029.txt) |
| **33** | 1.9395 × 10<sup>15</sup> | 441 · 8<sup>14</sup> | 0.8547 | B. Solomon | 2002 | [n033.txt](matrices/literature_witnesses/n033.txt) |
| **45** | 6.5862 × 10<sup>23</sup> | 89 · 11<sup>21</sup> | 0.8576 | Brent, Orrick et al. | 2009 | [n045.txt](matrices/literature_witnesses/n045.txt) |
| **49** | 6.3597 × 10<sup>26</sup> | 96 · 12<sup>23</sup> | 0.8123 | B. Solomon & W. P. Orrick | 2003 | [n049.txt](matrices/literature_witnesses/n049.txt) |
| **53** | 7.4092 × 10<sup>29</sup> | 105 · 13<sup>25</sup> | 0.7882 | Orrick, Solomon et al. | 2003 | [n053.txt](matrices/literature_witnesses/n053.txt) |
| **57** | 1.1730 × 10<sup>33</sup> | 133 · 14<sup>27</sup> | 0.8937 | J.-C. Meyrignac & B. Felgenhauer | 2005 | [n057.txt](matrices/literature_witnesses/n057.txt) |
| **65** | 3.1476 × 10<sup>39</sup> | 148 · 16<sup>31</sup> | 0.8144 | W. P. Orrick & B. Solomon | 2003 | [n065.txt](matrices/literature_witnesses/n065.txt) |
| **69** | 6.2394 × 10<sup>42</sup> | 155 · 17<sup>33</sup> | 0.7790 | Orrick, Solomon et al. | 2003 | [n069.txt](matrices/literature_witnesses/n069.txt) |
| **73** | 1.4965 × 10<sup>46</sup> | 174 · 18<sup>35</sup> | 0.8028 | B. Solomon & W. P. Orrick | 2003 | [n073.txt](matrices/literature_witnesses/n073.txt) |
| **77** | 3.7699 × 10<sup>49</sup> | 183 · 19<sup>37</sup> | 0.7787 | B. Solomon & W. P. Orrick | 2003 | [n077.txt](matrices/literature_witnesses/n077.txt) |
| **81** | 1.0775 × 10<sup>53</sup> | 196 · 20<sup>39</sup> | 0.7723 | N. Farmakis & S. Kounias | 1987 | [n081.txt](matrices/literature_witnesses/n081.txt) |
| **85** | 3.4902 × 10<sup>56</sup> | 4509 · 21<sup>40</sup> | 0.7865 | B. Solomon & W. P. Orrick | 2003 | [n085.txt](matrices/literature_witnesses/n085.txt) |
| **89** | 1.1657 × 10<sup>60</sup> | 220 · 22<sup>43</sup> | 0.7516 | N. Farmakis & S. Kounias | 1987 | none* |
| **93** | 4.5305 × 10<sup>63</sup> | 239 · 23<sup>45</sup> | 0.7640 | W. P. Orrick & B. Solomon | 2003 | [n093.txt](matrices/literature_witnesses/n093.txt) |
| **97** | 1.8974 × 10<sup>67</sup> | 256 · 24<sup>47</sup> | 0.7678 | B. Solomon & W. P. Orrick | 2003 | [n097.txt](matrices/literature_witnesses/n097.txt) |
| **101** | 8.9741 × 10<sup>70</sup> | 7110 · 25<sup>48</sup> | 0.8024 | W. P. Orrick & B. Solomon | 2003 | [n101.txt](matrices/literature_witnesses/n101.txt) |
| **105** | 3.8474 × 10<sup>74</sup> | — | 0.7022 | T. Rokicki | 2005 | [n105.txt](matrices/literature_witnesses/n105.txt) |
| **109** | 1.8775 × 10<sup>78</sup> | — | 0.6482 | T. Rokicki | 2005 | [n109.txt](matrices/literature_witnesses/n109.txt) |
| **117** | 7.2315 × 10<sup>85</sup> | — | 0.7186 | T. Rokicki | 2005 | [n117.txt](matrices/literature_witnesses/n117.txt) |

<sub>* On <b>n = 89</b>, the lower bound 220 · 22<sup>43</sup> was derived by Farmakis & Kounias (1987) from an order-88 max-excess Hadamard matrix; no explicit 89 × 89 matrix was archived in <code>MaxDet</code>.</sub>

### 2.2 Orders n ≡ 2 (mod 4) — Ehlich–Wojtas Regime (6 Open Orders)

*Upper bound: $B_2(n) = 2(n - 1)(n - 2)^{(n - 2)/2}$. All other 24 orders $n \equiv 2 \pmod 4$ up to $118$ are proven optimal.*

| n | D<sub>norm</sub> | Factored Form | R = \|det\|/B | Record Holder | Year | Witness |
| ---: | ---: | :--- | ---: | :--- | :---: | :---: |
| **34** | 9.0072 × 10<sup>15</sup> | 256 · 8<sup>15</sup> | 0.9697 | W. P. Orrick | 2005 | [n034.txt](matrices/literature_witnesses/n034.txt) |
| **58** | 6.5774 × 10<sup>33</sup> | — | 0.9345 | T. Rokicki; Kazmenko & Trofimov | 2005 | [n058.txt](matrices/literature_witnesses/n058.txt) |
| **70** | 4.6534 × 10<sup>43</sup> | 1156 · 17<sup>33</sup> | 0.9855 | W. P. Orrick & T. Rokicki | 2005 | [n070.txt](matrices/literature_witnesses/n070.txt) |
| **78** | 2.9042 × 10<sup>50</sup> | — | 0.9636 | T. Rokicki | 2005 | [n078.txt](matrices/literature_witnesses/n078.txt) |
| **94** | 3.9312 × 10<sup>64</sup> | — | 0.9695 | T. Rokicki | 2005 | [n094.txt](matrices/literature_witnesses/n094.txt) |
| **106** | 3.9414 × 10<sup>75</sup> | 2704 · 26<sup>51</sup> | 0.9905 | W. P. Orrick | 2005 | [n106.txt](matrices/literature_witnesses/n106.txt) |

### 2.3 Orders n ≡ 3 (mod 4) — Ehlich Block-Gram Regime (25 Open Orders)

*Upper bound: Ehlich $s$-block Gram bound ($s = 6$ for $23 \le n \le 59$, $s = 7$ for $63 \le n \le 119$). Bold 2026 entries indicate relative $D_{\text{norm}}$ gains over the 2005 Rokicki circulant baseline.*

| n | D<sub>norm</sub> | Factored Form | R = \|det\|/B | Record Holder | Year | Witness |
| ---: | ---: | :--- | ---: | :--- | :---: | :---: |
| **23** | 6.6267 × 10<sup>8</sup> | 662,671,875 | 0.9320 | Dowdeswell, Neubauer et al. | 2004 | [n023.txt](matrices/literature_witnesses/n023.txt) |
| **27** | 1.9809 × 10<sup>11</sup> | 546 · 6<sup>11</sup> | 0.9673 | H. Tamura | 2005 | [n027.txt](matrices/literature_witnesses/n027.txt) |
| **31** | 7.5961 × 10<sup>13</sup> | 784 · 7<sup>13</sup> | 0.9561 | H. Tamura | 2005 | [n031.txt](matrices/literature_witnesses/n031.txt) |
| **35** | 3.7436 × 10<sup>16</sup> | 1064 · 8<sup>15</sup> | 0.9379 | H. Tamura | 2011 | [n035.txt](matrices/literature_witnesses/n035.txt) |
| **39** | 2.4015 × 10<sup>19</sup> | 1440 · 9<sup>17</sup> | 0.9527 | H. Tamura | 2011 | [n039.txt](matrices/literature_witnesses/n039.txt) |
| **43** | 1.8900 × 10<sup>22</sup> | 1890 · 10<sup>19</sup> | 0.9663 | H. Tamura | 2011 | [n043.txt](matrices/literature_witnesses/n043.txt) |
| **47** | 1.6888 × 10<sup>25</sup> | 3,037,500 · 11<sup>18</sup> | 0.9239 | Orrick, Solomon et al. | 2003 | [n047.txt](matrices/literature_witnesses/n047.txt) |
| **51** | 1.7776 × 10<sup>28</sup> | — | 0.8776 | **Butbaia, Gukov et al.** (+3.11%) | **2026** | [n051.txt](matrices/literature_witnesses/n051.txt) |
| **55** | 2.2572 × 10<sup>31</sup> | — | 0.8593 | T. Rokicki; Kazmenko & Trofimov | 2005 | [n055.txt](matrices/literature_witnesses/n055.txt) |
| **59** | 3.8031 × 10<sup>34</sup> | 4312 · 14<sup>27</sup> | 0.9653 | H. Tamura | 2005 | [n059.txt](matrices/literature_witnesses/n059.txt) |
| **63** | 6.5387 × 10<sup>37</sup> | 5115 · 15<sup>29</sup> | 0.9662 | H. Tamura | 2011 | [n063.txt](matrices/literature_witnesses/n063.txt) |
| **67** | 1.1400 × 10<sup>41</sup> | — | 0.8635 | T. Rokicki | 2005 | [n067.txt](matrices/literature_witnesses/n067.txt) |
| **71** | 2.5553 × 10<sup>44</sup> | 31,186,944 · 17<sup>30</sup> | 0.8804 | T. Rokicki | 2005 | [n071.txt](matrices/literature_witnesses/n071.txt) |
| **75** | 6.1512 × 10<sup>47</sup> | — | 0.8613 | T. Rokicki | 2005 | [n075.txt](matrices/literature_witnesses/n075.txt) |
| **79** | 1.6796 × 10<sup>51</sup> | — | 0.8591 | T. Rokicki | 2005 | [n079.txt](matrices/literature_witnesses/n079.txt) |
| **83** | 5.0710 × 10<sup>54</sup> | — | 0.8561 | T. Rokicki | 2005 | [n083.txt](matrices/literature_witnesses/n083.txt) |
| **87** | 1.6851 × 10<sup>58</sup> | — | 0.8527 | T. Rokicki | 2005 | [n087.txt](matrices/literature_witnesses/n087.txt) |
| **91** | 6.1448 × 10<sup>61</sup> | — | 0.8501 | T. Rokicki | 2005 | [n091.txt](matrices/literature_witnesses/n091.txt) |
| **95** | 2.4387 × 10<sup>65</sup> | — | 0.8447 | T. Rokicki | 2005 | [n095.txt](matrices/literature_witnesses/n095.txt) |
| **99** | 1.0657 × 10<sup>69</sup> | — | 0.8496 | T. Rokicki | 2005 | [n099.txt](matrices/literature_witnesses/n099.txt) |
| **103** | 4.9829 × 10<sup>72</sup> | — | 0.8433 | **Butbaia, Gukov et al.** (+0.00004%) | **2026** | [n103.txt](matrices/literature_witnesses/n103.txt) |
| **107** | 2.5405 × 10<sup>76</sup> | — | 0.8445 | **Butbaia, Gukov et al.** (+0.44%) | **2026** | [n107.txt](matrices/literature_witnesses/n107.txt) |
| **111** | 1.3978 × 10<sup>80</sup> | — | 0.8469 | **Butbaia, Gukov et al.** (+1.26%) | **2026** | [n111.txt](matrices/literature_witnesses/n111.txt) |
| **115** | 8.2608 × 10<sup>83</sup> | — | 0.8488 | **Butbaia, Gukov et al.** (+1.83%) | **2026** | [n115.txt](matrices/literature_witnesses/n115.txt) |
| **119** | 5.2308 × 10<sup>87</sup> | — | 0.8503 | **Butbaia, Gukov et al.** (+2.12%) | **2026** | [n119.txt](matrices/literature_witnesses/n119.txt) |

<br>

<details>
<summary><strong>Complete 120-Order Master Table (1 <= n <= 120, Proven Optimal + Open)</strong></summary>

| n | mod 4 | Status | D<sub>norm</sub> | Factored Form | R = \|det\|/B | Attribution | Year | Ref | Witness |
| ---: | :---: | :---: | ---: | :--- | ---: | :--- | :---: | :---: | :---: |
| 1 | 1 | OPT | 1 | 1 | 1.0000 | Sylvester | 1867 | Syl1867 | — |
| 2 | 2 | OPT | 1 | 1 | 1.0000 | Sylvester | 1867 | Syl1867 | — |
| 3 | 3 | OPT | 1 | 1 | 1.0000 | Sylvester; Williamson | 1867 | Syl1867 | — |
| 4 | 0 | OPT | 2 | 2 | 1.0000 | Sylvester (Hadamard) | 1867 | Syl1867 | — |
| 5 | 1 | OPT | 3 | 3 | 1.0000 | Mood (Barba bound) | 1946 | Moo1946 | — |
| 6 | 2 | OPT | 5 | 5 | 1.0000 | Ehlich (EW bound) | 1964 | Ehl1964a | — |
| 7 | 3 | OPT | 9 | 9 | 0.9820 | Williamson | 1946 | Wil1946 | — |
| 8 | 0 | OPT | 32 | 32 | 1.0000 | Sylvester (Hadamard) | 1867 | Syl1867 | — |
| 9 | 1 | OPT | 56 | 56 | 0.8489 | Ehlich | 1964 | Ehl1964a | — |
| 10 | 2 | OPT | 144 | 144 | 1.0000 | Ehlich (EW bound) | 1964 | Ehl1964a | — |
| 11 | 3 | OPT | 320 | 320 | 0.9449 | Ehlich | 1964 | Ehl1964b | — |
| 12 | 0 | OPT | 1,458 | 1,458 | 1.0000 | Hadamard | 1893 | Had1893 | — |
| 13 | 1 | OPT | 3,645 | 3,645 | 1.0000 | Raghavarao (Barba bound) | 1959 | Rag1959 | — |
| 14 | 2 | OPT | 9,477 | 9,477 | 1.0000 | Ehlich (EW bound) | 1964 | Ehl1964a | — |
| 15 | 3 | OPT | 25,515 | 25,515 | 0.9707 | Ehlich & Zeller; Cohn | 1962 | EZ1962 | — |
| 16 | 0 | OPT | 131,072 | 131,072 | 1.0000 | Sylvester (Hadamard) | 1867 | Syl1867 | — |
| 17 | 1 | OPT | 327,680 | 327,680 | 0.8704 | Moyssiadis & Kounias | 1982 | MK1982 | — |
| 18 | 2 | OPT | 1.1141 × 10<sup>6</sup> | 1,114,112 | 1.0000 | Ehlich (EW bound) | 1964 | Ehl1964a | — |
| 19 | 3 | OPT | 3.4120 × 10<sup>6</sup> | 3,411,968 | 0.9750 | Cohn (1989); opt. Brent et al. | 2011 | BOOZ2011 | — |
| 20 | 0 | OPT | 1.9531 × 10<sup>7</sup> | 19,531,250 | 1.0000 | Hadamard | 1893 | Had1893 | — |
| 21 | 1 | OPT | 5.6641 × 10<sup>7</sup> | 56,640,625 | 0.9058 | Chadjipantelis, Kounias et al. | 1987 | CKM1987 | — |
| 22 | 2 | OPT | 1.9531 × 10<sup>8</sup> | 195,312,500 | 0.9524 | Ehlich (1964); opt. Brent | 2021 | Bre2021 | — |
| 23 | 3 | OPEN | 6.6267 × 10<sup>8</sup> | 662,671,875 | 0.9320 | Dowdeswell, Neubauer et al. | 2004 | DNST2004 | [txt](matrices/literature_witnesses/n023.txt) |
| 24 | 0 | OPT | 4.3536 × 10<sup>9</sup> | 12 · 6<sup>11</sup> | 1.0000 | Paley / Williamson (Hadamard) | 1933 | Pal1933 | — |
| 25 | 1 | OPT | 1.5237 × 10<sup>10</sup> | 42 · 6<sup>11</sup> | 1.0000 | Brouwer (Barba bound) | 1983 | Bro1983 | — |
| 26 | 2 | OPT | 5.4420 × 10<sup>10</sup> | 150 · 6<sup>11</sup> | 1.0000 | Ehlich (EW bound) | 1964 | Ehl1964a | — |
| 27 | 3 | OPEN | 1.9809 × 10<sup>11</sup> | 546 · 6<sup>11</sup> | 0.9673 | H. Tamura | 2005 | Tam2005 | [txt](matrices/literature_witnesses/n027.txt) |
| 28 | 0 | OPT | 1.3564 × 10<sup>12</sup> | 14 · 7<sup>13</sup> | 1.0000 | Paley / Williamson (Hadamard) | 1933 | Pal1933 | — |
| 29 | 1 | OPEN | 4.4292 × 10<sup>12</sup> | 320 · 7<sup>12</sup> | 0.8650 | B. Solomon | 2002 | Sol2002 | [txt](matrices/literature_witnesses/n029.txt) |
| 30 | 2 | OPT | 1.9668 × 10<sup>13</sup> | 203 · 7<sup>13</sup> | 1.0000 | Ehlich (EW bound) | 1964 | Ehl1964a | — |
| 31 | 3 | OPEN | 7.5961 × 10<sup>13</sup> | 784 · 7<sup>13</sup> | 0.9561 | H. Tamura | 2005 | Tam2005 | [txt](matrices/literature_witnesses/n031.txt) |
| 32 | 0 | OPT | 5.6295 × 10<sup>14</sup> | 16 · 8<sup>15</sup> | 1.0000 | Sylvester (Hadamard) | 1867 | Syl1867 | — |
| 33 | 1 | OPEN | 1.9395 × 10<sup>15</sup> | 441 · 8<sup>14</sup> | 0.8547 | B. Solomon | 2002 | Sol2002 | [txt](matrices/literature_witnesses/n033.txt) |
| 34 | 2 | OPEN | 9.0072 × 10<sup>15</sup> | 256 · 8<sup>15</sup> | 0.9697 | W. P. Orrick | 2005 | Orr2005 | [txt](matrices/literature_witnesses/n034.txt) |
| 35 | 3 | OPEN | 3.7436 × 10<sup>16</sup> | 1064 · 8<sup>15</sup> | 0.9379 | H. Tamura | 2011 | Tam2011 | [txt](matrices/literature_witnesses/n035.txt) |
| 36 | 0 | OPT | 3.0019 × 10<sup>17</sup> | 18 · 9<sup>17</sup> | 1.0000 | Paley / Williamson (Hadamard) | 1933 | Pal1933 | — |
| 37 | 1 | OPT | 1.2008 × 10<sup>18</sup> | 72 · 9<sup>17</sup> | 0.9363 | Brent, Orrick et al. | 2011 | BOOZ2011 | — |
| 38 | 2 | OPT | 5.5535 × 10<sup>18</sup> | 333 · 9<sup>17</sup> | 1.0000 | Ehlich (EW bound) | 1964 | Ehl1964a | — |
| 39 | 3 | OPEN | 2.4015 × 10<sup>19</sup> | 1440 · 9<sup>17</sup> | 0.9527 | H. Tamura | 2011 | Tam2011 | [txt](matrices/literature_witnesses/n039.txt) |
| 40 | 0 | OPT | 2.0000 × 10<sup>20</sup> | 20 · 10<sup>19</sup> | 1.0000 | Paley / Williamson (Hadamard) | 1933 | Pal1933 | — |
| 41 | 1 | OPT | 9.0000 × 10<sup>20</sup> | 90 · 10<sup>19</sup> | 1.0000 | Bridges, Hall & Hayden; van Trung | 1981 | BHH1981 | — |
| 42 | 2 | OPT | 4.1000 × 10<sup>21</sup> | 410 · 10<sup>19</sup> | 1.0000 | Yang (EW bound) | 1968 | Yan1968 | — |
| 43 | 3 | OPEN | 1.8900 × 10<sup>22</sup> | 1890 · 10<sup>19</sup> | 0.9663 | H. Tamura | 2011 | Tam2011 | [txt](matrices/literature_witnesses/n043.txt) |
| 44 | 0 | OPT | 1.6281 × 10<sup>23</sup> | 22 · 11<sup>21</sup> | 1.0000 | Paley / Williamson (Hadamard) | 1933 | Pal1933 | — |
| 45 | 1 | OPEN | 6.5862 × 10<sup>23</sup> | 89 · 11<sup>21</sup> | 0.8576 | Brent, Orrick et al. | 2009 | BOOZ2011 | [txt](matrices/literature_witnesses/n045.txt) |
| 46 | 2 | OPT | 3.6631 × 10<sup>24</sup> | 495 · 11<sup>21</sup> | 1.0000 | Yang (EW bound) | 1968 | Yan1968 | — |
| 47 | 3 | OPEN | 1.6888 × 10<sup>25</sup> | 3,037,500 · 11<sup>18</sup> | 0.9239 | Orrick, Solomon et al. | 2003 | OSDS2003 | [txt](matrices/literature_witnesses/n047.txt) |
| 48 | 0 | OPT | 1.5899 × 10<sup>26</sup> | 24 · 12<sup>23</sup> | 1.0000 | Paley / Williamson (Hadamard) | 1933 | Pal1933 | — |
| 49 | 1 | OPEN | 6.3597 × 10<sup>26</sup> | 96 · 12<sup>23</sup> | 0.8123 | B. Solomon & W. P. Orrick | 2003 | OS2007 | [txt](matrices/literature_witnesses/n049.txt) |
| 50 | 2 | OPT | 3.8953 × 10<sup>27</sup> | 588 · 12<sup>23</sup> | 1.0000 | Yang (EW bound) | 1968 | Yan1968 | — |
| 51 | 3 | OPEN | 1.7776 × 10<sup>28</sup> | — | 0.8776 | Butbaia, Gukov et al. | 2026 | BG2026 | [txt](matrices/literature_witnesses/n051.txt) |
| 52 | 0 | OPT | 1.8347 × 10<sup>29</sup> | 26 · 13<sup>25</sup> | 1.0000 | Paley / Williamson (Hadamard) | 1933 | Pal1933 | — |
| 53 | 1 | OPEN | 7.4092 × 10<sup>29</sup> | 105 · 13<sup>25</sup> | 0.7882 | Orrick, Solomon et al. | 2003 | OSDS2003 | [txt](matrices/literature_witnesses/n053.txt) |
| 54 | 2 | OPT | 4.8619 × 10<sup>30</sup> | 689 · 13<sup>25</sup> | 1.0000 | Yang (EW bound) | 1968 | Yan1968 | — |
| 55 | 3 | OPEN | 2.2572 × 10<sup>31</sup> | — | 0.8593 | T. Rokicki; Kazmenko & Trofimov | 2005 | Rok2005 | [txt](matrices/literature_witnesses/n055.txt) |
| 56 | 0 | OPT | 2.4695 × 10<sup>32</sup> | 28 · 14<sup>27</sup> | 1.0000 | Paley / Williamson (Hadamard) | 1933 | Pal1933 | — |
| 57 | 1 | OPEN | 1.1730 × 10<sup>33</sup> | 133 · 14<sup>27</sup> | 0.8937 | J.-C. Meyrignac & B. Felgenhauer | 2005 | Mey2005 | [txt](matrices/literature_witnesses/n057.txt) |
| 58 | 2 | OPEN | 6.5774 × 10<sup>33</sup> | — | 0.9345 | T. Rokicki; Kazmenko & Trofimov | 2005 | Rok2005 | [txt](matrices/literature_witnesses/n058.txt) |
| 59 | 3 | OPEN | 3.8031 × 10<sup>34</sup> | 4312 · 14<sup>27</sup> | 0.9653 | H. Tamura | 2005 | Tam2005 | [txt](matrices/literature_witnesses/n059.txt) |
| 60 | 0 | OPT | 3.8350 × 10<sup>35</sup> | 30 · 15<sup>29</sup> | 1.0000 | Paley / Williamson (Hadamard) | 1933 | Pal1933 | — |
| 61 | 1 | OPT | 2.1093 × 10<sup>36</sup> | 165 · 15<sup>29</sup> | 1.0000 | Brouwer (Barba bound) | 1983 | Bro1983 | — |
| 62 | 2 | OPT | 1.1697 × 10<sup>37</sup> | 915 · 15<sup>29</sup> | 1.0000 | Yang (EW bound) | 1968 | Yan1968 | — |
| 63 | 3 | OPEN | 6.5387 × 10<sup>37</sup> | 5115 · 15<sup>29</sup> | 0.9662 | H. Tamura | 2011 | Tam2011 | [txt](matrices/literature_witnesses/n063.txt) |
| 64 | 0 | OPT | 6.8056 × 10<sup>38</sup> | 32 · 16<sup>31</sup> | 1.0000 | Sylvester (Hadamard) | 1867 | Syl1867 | — |
| 65 | 1 | OPEN | 3.1476 × 10<sup>39</sup> | 148 · 16<sup>31</sup> | 0.8144 | W. P. Orrick & B. Solomon | 2003 | OS2007 | [txt](matrices/literature_witnesses/n065.txt) |
| 66 | 2 | OPT | 2.2118 × 10<sup>40</sup> | 1040 · 16<sup>31</sup> | 1.0000 | Yang (EW bound) | 1968 | Yan1968 | — |
| 67 | 3 | OPEN | 1.1400 × 10<sup>41</sup> | — | 0.8635 | T. Rokicki | 2005 | Rok2005 | [txt](matrices/literature_witnesses/n067.txt) |
| 68 | 0 | OPT | 1.3687 × 10<sup>42</sup> | 34 · 17<sup>33</sup> | 1.0000 | Paley / Williamson (Hadamard) | 1933 | Pal1933 | — |
| 69 | 1 | OPEN | 6.2394 × 10<sup>42</sup> | 155 · 17<sup>33</sup> | 0.7790 | Orrick, Solomon et al. | 2003 | OSDS2003 | [txt](matrices/literature_witnesses/n069.txt) |
| 70 | 2 | OPEN | 4.6534 × 10<sup>43</sup> | 1156 · 17<sup>33</sup> | 0.9855 | W. P. Orrick & T. Rokicki | 2005 | Orr2005 | [txt](matrices/literature_witnesses/n070.txt) |
| 71 | 3 | OPEN | 2.5553 × 10<sup>44</sup> | 31,186,944 · 17<sup>30</sup> | 0.8804 | T. Rokicki | 2005 | Rok2005 | [txt](matrices/literature_witnesses/n071.txt) |
| 72 | 0 | OPT | 3.0963 × 10<sup>45</sup> | 36 · 18<sup>35</sup> | 1.0000 | Paley / Williamson (Hadamard) | 1933 | Pal1933 | — |
| 73 | 1 | OPEN | 1.4965 × 10<sup>46</sup> | 174 · 18<sup>35</sup> | 0.8028 | B. Solomon & W. P. Orrick | 2003 | OS2007 | [txt](matrices/literature_witnesses/n073.txt) |
| 74 | 2 | OPT | 1.1301 × 10<sup>47</sup> | 1314 · 18<sup>35</sup> | 1.0000 | Cohn (EW bound) | 1989 | Coh1989 | — |
| 75 | 3 | OPEN | 6.1512 × 10<sup>47</sup> | — | 0.8613 | T. Rokicki | 2005 | Rok2005 | [txt](matrices/literature_witnesses/n075.txt) |
| 76 | 0 | OPT | 7.8283 × 10<sup>48</sup> | 38 · 19<sup>37</sup> | 1.0000 | Paley / Williamson (Hadamard) | 1933 | Pal1933 | — |
| 77 | 1 | OPEN | 3.7699 × 10<sup>49</sup> | 183 · 19<sup>37</sup> | 0.7787 | B. Solomon & W. P. Orrick | 2003 | OS2007 | [txt](matrices/literature_witnesses/n077.txt) |
| 78 | 2 | OPEN | 2.9042 × 10<sup>50</sup> | — | 0.9636 | T. Rokicki | 2005 | Rok2005 | [txt](matrices/literature_witnesses/n078.txt) |
| 79 | 3 | OPEN | 1.6796 × 10<sup>51</sup> | — | 0.8591 | T. Rokicki | 2005 | Rok2005 | [txt](matrices/literature_witnesses/n079.txt) |
| 80 | 0 | OPT | 2.1990 × 10<sup>52</sup> | 40 · 20<sup>39</sup> | 1.0000 | Paley / Williamson (Hadamard) | 1933 | Pal1933 | — |
| 81 | 1 | OPEN | 1.0775 × 10<sup>53</sup> | 196 · 20<sup>39</sup> | 0.7723 | N. Farmakis & S. Kounias | 1987 | FK1987 | [txt](matrices/literature_witnesses/n081.txt) |
| 82 | 2 | OPT | 8.9060 × 10<sup>53</sup> | 1620 · 20<sup>39</sup> | 1.0000 | Cohn (EW bound) | 1989 | Coh1989 | — |
| 83 | 3 | OPEN | 5.0710 × 10<sup>54</sup> | — | 0.8561 | T. Rokicki | 2005 | Rok2005 | [txt](matrices/literature_witnesses/n083.txt) |
| 84 | 0 | OPT | 6.8272 × 10<sup>55</sup> | 42 · 21<sup>41</sup> | 1.0000 | Paley / Williamson (Hadamard) | 1933 | Pal1933 | — |
| 85 | 1 | OPEN | 3.4902 × 10<sup>56</sup> | 4509 · 21<sup>40</sup> | 0.7865 | B. Solomon & W. P. Orrick | 2003 | OS2007 | [txt](matrices/literature_witnesses/n085.txt) |
| 86 | 2 | OPT | 2.9015 × 10<sup>57</sup> | 1785 · 21<sup>41</sup> | 1.0000 | Chadjipantelis & Kounias (EW bound) | 1985 | CK1985 | — |
| 87 | 3 | OPEN | 1.6851 × 10<sup>58</sup> | — | 0.8527 | T. Rokicki | 2005 | Rok2005 | [txt](matrices/literature_witnesses/n087.txt) |
| 88 | 0 | OPT | 2.3315 × 10<sup>59</sup> | 44 · 22<sup>43</sup> | 1.0000 | Paley / Williamson (Hadamard) | 1933 | Pal1933 | — |
| 89 | 1 | OPEN | 1.1657 × 10<sup>60</sup> | 220 · 22<sup>43</sup> | 0.7516 | N. Farmakis & S. Kounias | 1987 | FK1987 | — |
| 90 | 2 | OPT | 1.0375 × 10<sup>61</sup> | 1958 · 22<sup>43</sup> | 1.0000 | Cohn (EW bound) | 1989 | Coh1989 | — |
| 91 | 3 | OPEN | 6.1448 × 10<sup>61</sup> | — | 0.8501 | T. Rokicki | 2005 | Rok2005 | [txt](matrices/literature_witnesses/n091.txt) |
| 92 | 0 | OPT | 8.7199 × 10<sup>62</sup> | 46 · 23<sup>45</sup> | 1.0000 | Baumert, Golomb & Hall | 1962 | BGH1962 | — |
| 93 | 1 | OPEN | 4.5305 × 10<sup>63</sup> | 239 · 23<sup>45</sup> | 0.7640 | W. P. Orrick & B. Solomon | 2003 | OS2007 | [txt](matrices/literature_witnesses/n093.txt) |
| 94 | 2 | OPEN | 3.9312 × 10<sup>64</sup> | — | 0.9695 | T. Rokicki | 2005 | Rok2005 | [txt](matrices/literature_witnesses/n094.txt) |
| 95 | 3 | OPEN | 2.4387 × 10<sup>65</sup> | — | 0.8447 | T. Rokicki | 2005 | Rok2005 | [txt](matrices/literature_witnesses/n095.txt) |
| 96 | 0 | OPT | 3.5577 × 10<sup>66</sup> | 48 · 24<sup>47</sup> | 1.0000 | Paley / Williamson (Hadamard) | 1933 | Pal1933 | — |
| 97 | 1 | OPEN | 1.8974 × 10<sup>67</sup> | 256 · 24<sup>47</sup> | 0.7678 | B. Solomon & W. P. Orrick | 2003 | OS2007 | [txt](matrices/literature_witnesses/n097.txt) |
| 98 | 2 | OPT | 1.7255 × 10<sup>68</sup> | 2328 · 24<sup>47</sup> | 1.0000 | Cohn (EW bound) | 1989 | Coh1989 | — |
| 99 | 3 | OPEN | 1.0657 × 10<sup>69</sup> | — | 0.8496 | T. Rokicki | 2005 | Rok2005 | [txt](matrices/literature_witnesses/n099.txt) |
| 100 | 0 | OPT | 1.5777 × 10<sup>70</sup> | 50 · 25<sup>49</sup> | 1.0000 | Paley / Williamson (Hadamard) | 1933 | Pal1933 | — |
| 101 | 1 | OPEN | 8.9741 × 10<sup>70</sup> | 7110 · 25<sup>48</sup> | 0.8024 | W. P. Orrick & B. Solomon | 2003 | OS2007 | [txt](matrices/literature_witnesses/n101.txt) |
| 102 | 2 | OPT | 7.9675 × 10<sup>71</sup> | 2525 · 25<sup>49</sup> | 1.0000 | Cohn (EW bound) | 1989 | Coh1989 | — |
| 103 | 3 | OPEN | 4.9829 × 10<sup>72</sup> | — | 0.8433 | Butbaia, Gukov et al. | 2026 | BG2026 | [txt](matrices/literature_witnesses/n103.txt) |
| 104 | 0 | OPT | 7.5796 × 10<sup>73</sup> | 52 · 26<sup>51</sup> | 1.0000 | Paley / Williamson (Hadamard) | 1933 | Pal1933 | — |
| 105 | 1 | OPEN | 3.8474 × 10<sup>74</sup> | — | 0.7022 | T. Rokicki | 2005 | Rok2005 | [txt](matrices/literature_witnesses/n105.txt) |
| 106 | 2 | OPEN | 3.9414 × 10<sup>75</sup> | 2704 · 26<sup>51</sup> | 0.9905 | W. P. Orrick | 2005 | Orr2005 | [txt](matrices/literature_witnesses/n106.txt) |
| 107 | 3 | OPEN | 2.5405 × 10<sup>76</sup> | — | 0.8445 | Butbaia, Gukov et al. | 2026 | BG2026 | [txt](matrices/literature_witnesses/n107.txt) |
| 108 | 0 | OPT | 3.9325 × 10<sup>77</sup> | 54 · 27<sup>53</sup> | 1.0000 | Paley / Williamson (Hadamard) | 1933 | Pal1933 | — |
| 109 | 1 | OPEN | 1.8775 × 10<sup>78</sup> | — | 0.6482 | T. Rokicki | 2005 | Rok2005 | [txt](matrices/literature_witnesses/n109.txt) |
| 110 | 2 | OPT | 2.1432 × 10<sup>79</sup> | 2943 · 27<sup>53</sup> | 1.0000 | Fletcher & Seberry (EW bound) | 2001 | FS2001 | — |
| 111 | 3 | OPEN | 1.3978 × 10<sup>80</sup> | — | 0.8469 | Butbaia, Gukov et al. | 2026 | BG2026 | [txt](matrices/literature_witnesses/n111.txt) |
| 112 | 0 | OPT | 2.1973 × 10<sup>81</sup> | 56 · 28<sup>55</sup> | 1.0000 | Paley / Williamson (Hadamard) | 1933 | Pal1933 | — |
| 113 | 1 | OPT | 1.6479 × 10<sup>82</sup> | 420 · 28<sup>55</sup> | 1.0000 | Brouwer (Barba bound) | 1983 | Bro1983 | — |
| 114 | 2 | OPT | 1.2414 × 10<sup>83</sup> | 3164 · 28<sup>55</sup> | 1.0000 | Koukouvinos, Kounias & Seberry | 1989 | KKS1989 | — |
| 115 | 3 | OPEN | 8.2608 × 10<sup>83</sup> | — | 0.8488 | Butbaia, Gukov et al. | 2026 | BG2026 | [txt](matrices/literature_witnesses/n115.txt) |
| 116 | 0 | OPT | 1.3186 × 10<sup>85</sup> | 58 · 29<sup>57</sup> | 1.0000 | Baumert | 1966 | Bau1966 | — |
| 117 | 1 | OPEN | 7.2315 × 10<sup>85</sup> | — | 0.7186 | T. Rokicki | 2005 | Rok2005 | [txt](matrices/literature_witnesses/n117.txt) |
| 118 | 2 | OPT | 7.7138 × 10<sup>86</sup> | 3393 · 29<sup>57</sup> | 1.0000 | Fletcher, Koukouvinos & Seberry | 2002 | FKS2002 | — |
| 119 | 3 | OPEN | 5.2308 × 10<sup>87</sup> | — | 0.8503 | Butbaia, Gukov et al. | 2026 | BG2026 | [txt](matrices/literature_witnesses/n119.txt) |
| 120 | 0 | OPT | 8.4782 × 10<sup>88</sup> | 60 · 30<sup>59</sup> | 1.0000 | Paley / Williamson (Hadamard) | 1933 | Pal1933 | — |
</details>

---

## 3. Exact Verification & Testing Your Own Matrices (`verify.py`)

### Why Floating-Point Determinants Cannot Certify Records
Standard IEEE-754 `float64` routines (`numpy.linalg.det` and `numpy.linalg.slogdet`) carry only 53 bits ($\approx 15.9$ decimal digits) of mantissa precision, whereas $D_{\text{norm}}(n)$ reaches 88 decimal digits by $n = 119$, and distinct constructions can agree to 6+ leading digits (for example at $n = 103$, where the 2026 improvement over Rokicki's 2005 circulant matrix is $+0.00004\%$ in the 7th significant digit). Every entry in this repository is therefore verified in **exact arbitrary-precision integer arithmetic**:

1. **Domain & Dimension Check:** Verifies that $X$ has exact shape $n \times n$ and every entry belongs strictly to $\{-1, +1\}$.
2. **Fraction-Free Bareiss Elimination:** Computes exact $\det X \in \mathbb{Z}$ via Bareiss integer elimination (`det_bareiss_int` in [`verify.py`](verify.py), zero external dependencies; optional independent cross-check via `sympy.Matrix.det(method='bareiss')` with `--sympy`).
3. **Exact Scaling & Upper-Bound Inequality:** Confirms $|\det X| \equiv 0 \pmod{2^{n-1}}$, extracts $D_{\text{norm}} = |\det X| / 2^{n-1}$, and checks the exact rational inequality $|\det X|^2 \le B_{\text{class}}(n)^2$ (using the residue-appropriate Hadamard, Barba, Ehlich–Wojtas, or $s$-block Ehlich bound with $s = 3, 5, 6, 7$).

### Checking an External Candidate Matrix in One Command
Researchers and automated search agents can validate any candidate matrix against the exact residue-class bound and the live `1 <= n <= 120` catalog without installing third-party packages:

```bash
python3 verify.py --check path/to/candidate_matrix.txt          # Pure-Python Bareiss (~1-120 ms)
python3 verify.py --check path/to/candidate_matrix.json --sympy # Dual Bareiss + SymPy check
```

`--check` automatically accepts four common community formats:
- **Compact `+/-` sign grids** (`n` lines of `+` and `-`, used in `matrices/*.txt` and Orrick's `MaxDet` archive),
- **JSON 2D arrays** (`[[1, -1, ...], ...]` or `{"matrix": [[...]]}`, matching **EinsteinArena** payloads),
- **CSV or whitespace-separated `{-1, +1}` integer matrices** (`numpy.savetxt(..., fmt='%d')`),
- **Binary `{0, 1}` matrices of order $m = n - 1$** ([OEIS A003432](https://oeis.org/A003432) convention, lifted constructively via $X = \bigl[\begin{smallmatrix} 1 & \mathbf{1}^T \\ \mathbf{1} & J - 2B \end{smallmatrix}\bigr]$ with $D_{\text{norm}}(m+1) = |\det B|$).

<details>
<summary><strong>CSV & JSON Column Guide (<code>status</code> vs. <code>attains_residue_bound</code>)</strong></summary>

[`catalog/hadamard_maxdet_n1_to_120.csv`](catalog/hadamard_maxdet_n1_to_120.csv) and [`catalog/hadamard_maxdet_n1_to_120.json`](catalog/hadamard_maxdet_n1_to_120.json) distinguish between **optimality** (`status`) and **bound tightness** (`attains_residue_bound`):

- **`status = PROVEN_OPTIMAL` + `attains_residue_bound = True` (61 orders):** A sign matrix achieves $100\%$ of its residue-class upper bound ($R = 1$), immediately certifying optimality (all 30 Hadamard orders $n \equiv 0 \pmod 4$, $n \in \{1, 2, 3\}$, 6 Barba orders $5, 13, 25, 41, 61, 113$, and 22 Ehlich–Wojtas orders).
- **`status = PROVEN_OPTIMAL` + `attains_residue_bound = False` (9 orders: `7, 9, 11, 15, 17, 19, 21, 22, 37`):** The closed-form residue-class formula is arithmetically unattainable over $\{-1, +1\}$ ($R < 1$), yet an exhaustive computer-assisted or algebraic proof establishes that no larger determinant exists for that $n$ (e.g., [Brent et al. 2011](https://arxiv.org/abs/1112.4160) for $n = 19, 37$; [Brent 2021](https://oeis.org/A003432) for $n = 22$).
- **`status = OPEN` + `attains_residue_bound = False` (50 orders, `23 <= n <= 119`):** No optimality proof is published and $R < 1$; `D_norm_exact` records the best verified constructive lower bound.
</details>

---

## 4. Primary Sources, Related Benchmarks & Attribution

### Archival Databases & Integer Sequences
- **W. P. Orrick & B. Solomon — `MaxDet` (`indiana.edu/~maxdet`, 2002–2012):** Preserved via the [Indiana University Archive-It mirror](http://wayback.archive-it.org/219/20190316001107/http://www.indiana.edu/~maxdet/), the [Wayback full summary table (`fullPage.shtml`)](https://web.archive.org/web/20131002153624/http://www.indiana.edu/~maxdet/fullPage.shtml), and the [theoretical upper-bounds guide (`bounds.html`)](https://web.archive.org/web/20131002153511/http://www.indiana.edu/~maxdet/bounds.html).
- **OEIS Sequences:** [OEIS A003432](https://oeis.org/A003432) (normalized integer $D_{\text{norm}}(n) = |\det X| / 2^{n-1}$, equal to the maximal $(n-1)\times(n-1)$ binary $\{0,1\}$ determinant) · [OEIS A003433](https://oeis.org/A003433) (raw maximal $\{-1,+1\}$ determinant $|\det X|$).
- **Two-Circulant RL Dataset (`n ≡ 3 mod 4`):** G. Butbaia, S. Gukov et al., [`Math-AI-Caltech/hadamard-maxdet`](https://github.com/Math-AI-Caltech/hadamard-maxdet) ([arXiv:2608.22518](https://arxiv.org/abs/2608.22518) · [Zenodo 22024187](https://zenodo.org/records/22024187), Oct 2026) covering orders `51, 103, 107, 111, 115, 119`.

### Active Open-Problem Benchmarks & Agentic Challenges
- **[EinsteinArena](https://einsteinarena.com/) (*Hadamard Maximal Determinant*):** Hosts a live benchmark specifically on **order $n = 51$**, where the current record is $D_{\text{norm}}(51) = 17{,}776{,}121{,}037{,}665{,}193{,}653{,}653{,}203{,}125$ ($R = 0.877569$, [`matrices/literature_witnesses/n051.txt`](matrices/literature_witnesses/n051.txt)). Candidate $51 \times 51$ JSON matrices can be pre-verified locally with `python3 verify.py --check candidate_51.json`.
- **[`elder-plinius/V3R1T4S`](https://github.com/elder-plinius/V3R1T4S) (`hadamard-maxdet-records`):** Open verification benchmark tracking strict determinant improvements and maximality certificates across all 50 open orders $n \not\equiv 0 \pmod 4$.
- **[`teorth/optimizationproblems`](https://github.com/teorth/optimizationproblems):** Terence Tao's curated catalog of open extremal and combinatorial optimization benchmarks.

### Expository Overviews & Constructive Hadamard Libraries
- **Mathematical Surveys:** [Wikipedia — *Hadamard's maximal determinant problem*](https://en.wikipedia.org/wiki/Hadamard%27s_maximal_determinant_problem) · [Wolfram MathWorld — *Hadamard's Maximum Determinant Problem*](https://mathworld.wolfram.com/HadamardsMaximumDeterminantProblem.html) · [R. P. Brent's Maximal Determinant Papers](https://maths-people.anu.edu.au/~brent/pub/pubsall.html).
- **Constructive Matrix Libraries (Order $4k$ Seeds):** [N. J. A. Sloane's Library of Hadamard Matrices](http://neilsloane.com/hadamard/) · [SageMath Combinatorial Designs (`hadamard_matrix`)](https://doc.sagemath.org/html/en/reference/combinat/sage/combinat/matrices/hadamard_matrix.html).
- **Agent & LLM Discovery Manifests:** [`AGENTS.md`](AGENTS.md) (autonomous agent primer, target filter & verification rules) · [`llms.txt`](llms.txt) ([`llmstxt.org`](https://llmstxt.org/) machine-readable index) · [`CITATION.cff`](CITATION.cff).

Corrections, missing historical references, or new verified matrices (`python3 verify.py --check your_matrix.txt`) are welcome via GitHub Issues or Pull Requests.

<details>
<summary><strong>Full Bibliographic Key ([Syl1867] – [BG2026])</strong></summary>

- **[Syl1867]** J. J. Sylvester (1867), Thoughts on inverse orthogonal matrices, Phil. Mag. 34, 461-475.
- **[Had1893]** J. Hadamard (1893), Resolution d'une question relative aux determinants, Bull. Sci. Math. 17, 240-246.
- **[Pal1933]** R. E. A. C. Paley (1933), On orthogonal matrices, J. Math. Phys. 12, 311-320.
- **[Bar1933]** G. Barba (1933), Intorno al teorema di Hadamard sui determinanti a valore massimo, Giorn. Mat. Battaglini 71, 70-86.
- **[Wil1944]** J. Williamson (1944), Hadamard's determinant theorem and the sum of four squares, Duke Math. J. 11, 65-81.
- **[Wil1946]** J. Williamson (1946), Determinants whose elements are 0 and 1, Amer. Math. Monthly 53, 427-434.
- **[Moo1946]** A. M. Mood (1946), On Hotelling's weighing problem, Ann. Math. Statist. 17, 432-446.
- **[Rag1959]** D. Raghavarao (1959), Some optimum weighing designs, Ann. Math. Statist. 30, 295-303.
- **[EZ1962]** H. Ehlich and K. Zeller (1962), Binare Matrizen, Z. Angew. Math. Mech. 42, T20-T21.
- **[BGH1962]** L. Baumert, S. W. Golomb, and M. Hall Jr. (1962), Discovery of an Hadamard matrix of order 92, Bull. Amer. Math. Soc. 68, 237-238.
- **[Ehl1964a]** H. Ehlich (1964), Determinantenabschatzungen fur binare Matrizen, Math. Z. 83, 123-132.
- **[Ehl1964b]** H. Ehlich (1964), Determinantenabschatzungen fur binare Matrizen mit n = 3 mod 4, Math. Z. 84, 438-447.
- **[Woj1964]** M. Wojtas (1964), On Hadamard's inequality for the determinants of order non-divisible by 4, Colloq. Math. 12, 73-83.
- **[Bau1966]** L. D. Baumert (1966), Hadamard matrices of orders 116 and 232, Bull. Amer. Math. Soc. 72, 237.
- **[Yan1968]** C. H. Yang (1966-1969), Some designs for maximal (+1,-1)-determinant of order n = 2 mod 4, Math. Comp. 20/22/23.
- **[BHH1981]** W. G. Bridges, M. Hall Jr., and M. L. Hayden (1981), Codes and designs, J. Combin. Theory Ser. A 31, 155-174; T. van Trung (1982).
- **[MK1982]** C. Moyssiadis and S. Kounias (1982), The exact D-optimal first order saturated design with 17 observations, J. Statist. Plann. Inference 7, 13-27.
- **[Bro1983]** A. E. Brouwer (1983), An infinite series of symmetric designs, Math. Centrum Amsterdam Report ZW 202/83.
- **[CK1985]** T. Chadjipantelis and S. Kounias (1985), Supplementary difference sets and D-optimal designs for n = 2 mod 4, Discrete Math. 57, 211-216.
- **[CKM1987]** T. Chadjipantelis, S. Kounias, and C. Moyssiadis (1987), The maximum determinant of 21x21 (+1,-1)-matrices and D-optimal designs, J. Statist. Plann. Inference 16, 167-178.
- **[FK1987]** N. Farmakis and S. Kounias (1987), The excess of Hadamard matrices and optimal designs, Discrete Math. 67, 165-176.
- **[Coh1989]** J. H. E. Cohn (1989), On determinants with elements +-1, II, Bull. London Math. Soc. 21, 36-42.
- **[KKS1989]** C. Koukouvinos, S. Kounias, and J. Seberry (1989), Further Hadamard matrices with maximal excess and new D-optimal designs for n = 2 mod 4, Utilitas Math. 36, 135-150.
- **[FS2001]** R. J. Fletcher and J. Seberry (2001), New D-optimal designs of order 110, Australas. J. Combin. 23, 49-52.
- **[FKS2002]** R. J. Fletcher, C. Koukouvinos, and J. Seberry (2002), New skew-Hadamard matrices of order 4*59 and D-optimal designs of order 2*59, Discrete Math. 286 (2004), 251-253.
- **[Sol2002]** B. Solomon (2002), Large-determinant sign matrices of orders 29 and 33, reported in W. P. Orrick's [MaxDet catalog](http://wayback.archive-it.org/219/20190316001107/http://www.indiana.edu/~maxdet/).
- **[OSDS2003]** W. P. Orrick, B. Solomon, R. Dowdeswell, and W. D. Smith (2003), New lower bounds for the maximal determinant problem, [arXiv:math/0304410](https://arxiv.org/abs/math/0304410).
- **[OS2007]** W. P. Orrick and B. Solomon (2003; publ. 2007), Large-determinant sign matrices of order 4k+1, Discrete Math. 307, 226-236 ([arXiv:math/0311292](https://arxiv.org/abs/math/0311292)).
- **[DNST2004]** A. Dowdeswell, M. G. Neubauer, B. Solomon, and K. Tumer (2004), Lower bound for order 23 in W. P. Orrick's [MaxDet catalog](https://web.archive.org/web/20131002153624/http://www.indiana.edu/~maxdet/fullPage.shtml).
- **[Orr2005]** W. P. Orrick (2005), Switching operations for Hadamard matrices and D-optimal designs of orders 34, 70, 106, SIAM J. Discrete Math. 22 (2008), 31-50 ([arXiv:math/0507515](https://arxiv.org/abs/math/0507515)).
- **[Rok2005]** T. Rokicki (Apr 2005), New determinant records from circulant pairs (orders 51-119), archived in W. P. Orrick's [MaxDet catalog](https://web.archive.org/web/20131002153624/http://www.indiana.edu/~maxdet/fullPage.shtml).
- **[Mey2005]** J.-C. Meyrignac and B. Felgenhauer (Apr 2005), Order 57 sign matrix record, archived in W. P. Orrick's [MaxDet catalog](https://web.archive.org/web/20131002153624/http://www.indiana.edu/~maxdet/fullPage.shtml).
- **[Tam2005]** H. Tamura (Aug 2005), Maximal determinant records for orders 27, 31, and 59, archived in W. P. Orrick's [MaxDet catalog](https://web.archive.org/web/20131002153624/http://www.indiana.edu/~maxdet/fullPage.shtml).
- **[BOOZ2011]** R. P. Brent, W. P. Orrick, J.-a. H. Osborn, and P. Zimmermann (2009-2011), Maximal determinants and saturated D-optimal designs of orders 19 and 37 (and order 45 construction), [arXiv:1112.4160](https://arxiv.org/abs/1112.4160).
- **[Tam2011]** H. Tamura (Jan 2011), Maximal determinant records for orders 35, 39, 43, and 63, archived in W. P. Orrick's [MaxDet catalog](https://web.archive.org/web/20131002153624/http://www.indiana.edu/~maxdet/fullPage.shtml).
- **[DK2015]** D. Z. Djokovic and I. S. Kotsireas (2015), D-optimal matrices of orders 118, 138, 150, 154 and 174, Algebraic Design Theory and Hadamard Matrices, Springer, 61-72 ([arXiv:1501.07640](https://arxiv.org/abs/1501.07640)).
- **[Bre2021]** R. P. Brent (2021), Optimality confirmation for binary order 21 (sign order 22, a(21) = 195312500), [OEIS A003432](https://oeis.org/A003432).
- **[BG2026]** G. Butbaia, J. Tan, P. Vipulanandan, X. Huang, T. Saunders-A'Court, L. Fagan, D. Passaro, M. Tarquini, and S. Gukov (Oct 2026), Reinforcement learning for the Hadamard maximal determinant problem, [arXiv:2608.22518](https://arxiv.org/abs/2608.22518) · [Zenodo 22024187](https://zenodo.org/records/22024187) · [Math-AI-Caltech/hadamard-maxdet](https://github.com/Math-AI-Caltech/hadamard-maxdet).
</details>

---

<sub>Curated and verified by **we\*** (`@mizopl`) · Released under the [MIT License](LICENSE).</sub>
