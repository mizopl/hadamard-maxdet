# Hadamard Maximal Determinant Catalog & New Constructive Records (1 <= n <= 120)

Unified, offline-verifiable catalog of maximal determinants of $n \times n$ sign matrices $X \in \{-1, +1\}^{n \times n}$ across all orders $1 \le n \le 120$, combining historical records from William P. Orrick's archived [`MaxDet` database](http://wayback.archive-it.org/219/20190316001107/http://www.indiana.edu/~maxdet/) (2002–2012), Tomas Rokicki's circulant constructions (2005), Hiroki Tamura's constructions (2005–2011), Chasiotis, Kounias, and Farmakis's optimality proof for $n = 22$ (*Discrete Math.*, 2018; [OEIS A003432](https://oeis.org/A003432)), Browne et al.'s survey ([arXiv:2104.06756](https://arxiv.org/abs/2104.06756), 2021), the two-circulant reinforcement-learning records of Butbaia, Gukov et al. ([arXiv:2608.22518](https://arxiv.org/abs/2608.22518), [`Math-AI-Caltech/hadamard-maxdet`](https://github.com/Math-AI-Caltech/hadamard-maxdet), Oct 2026), and **8 new constructive records** ($n = 81, 89, 93, 97, 101, 105, 109, 117$) established by **we\*** (`@mizopl`, Nov 2026).

<details>
<summary><strong>Historical Background, Practical Significance &amp; Plain-Language Engineering Guide (<code>*</code> Legend)</strong></summary>
<br>
<table>
<tr>
<td>
<img src="assets/jacques_hadamard.jpg" alt="Jacques Hadamard (1865–1963)" title="Jacques Hadamard (1865–1963), Archives École Polytechnique" width="105" align="right">

<b>Historical Origin &amp; Practical Significance.</b> Posed by <a href="https://en.wikipedia.org/wiki/Hadamard%27s_maximal_determinant_problem">Jacques Hadamard (1893)</a> following <a href="https://en.wikipedia.org/wiki/Sylvester%27s_construction">J. J. Sylvester's (1867)</a> dyadic constructions, the maximal determinant problem asks for the largest possible |det <i>X</i>| among <i>n</i> &times; <i>n</i> sign matrices <i>X</i> &isin; {&minus;1, +1}<sup><i>n</i>&times;<i>n</i></sup>. Hadamard's orthogonal ceiling |det <i>X</i>| &le; <i>n</i><sup><i>n</i>/2</sup> is attainable only for <i>n</i> = 1, 2, and multiples of 4 (where existence for every <i>n</i> &equiv; 0 mod 4 remains the open Hadamard conjecture, currently settled up to <i>n</i> &lt; 668); on all remaining orders <i>n</i> &equiv; 1, 2, 3 (mod 4), exact column orthogonality is arithmetically impossible and sharper bounds proved by Barba (1933), Ehlich (1964), and Wojtas (1964) govern the extremal geometry. In applied statistics and metrology, maximal-determinant &plusmn;1 matrices are <b><i>D</i>-optimal weighing and two-level factorial designs</b> (Hotelling 1944; Kiefer): maximizing |det <i>X</i>| minimizes the generalized covariance det(<i>X</i><sup>T</sup><i>X</i>)<sup>&minus;1</sup> = |det <i>X</i>|<sup>&minus;2</sup> and the joint confidence-ellipsoid volume when estimating <i>n</i> unknown parameters under bounded &plusmn;1 probes. Through the Williamson {0, 1} bijection, the same matrices maximize multiplex signal-to-noise gain (Fellgett's advantage) in <b>Hadamard-transform spectroscopy and coded-aperture imaging</b> when detector geometry fixes <i>n</i> &not;&equiv; 0 (mod 4), supply low-sidelobe sequences for <b>radar and spread-spectrum communications</b>, and calibrate <b>complete-pivoting growth bounds</b> in numerical linear algebra.

</td>
</tr>
<tr>
<td>

<b>Plain-Language Intuition &amp; Why <code>*</code> Marks Structured Constructions <i>(a layperson's note)</i>.</b> <i>Written from a non-specialist perspective:</i> imagine weighing <i>n</i> objects simultaneously on a two-pan balance (<code>+1</code> = left pan, <code>&minus;1</code> = right pan) or shining light through an <i>n</i>-slit mask across <i>n</i> steps instead of measuring one slot at a time. A larger determinant |det <i>X</i>| keeps those <i>n</i> combined measurements as mutually independent as discrete &plusmn;1 switches allow, shrinking reconstruction noise without paying for extra runs. In real hardware, however, <b>not all orders or constructions are equally practical</b>: whenever engineers can freely choose <i>n</i>, they pick <b>powers of two</b> (<i>n</i> = 2<sup><i>m</i></sup>, marked <code>*</code> on <i>n</i> in the 120-order table — fast Walsh&ndash;Hadamard butterflies using only additions and subtractions) or <b>multiples of 4</b> (<code>0 mod 4*</code>, exact orthogonality <i>R</i> = 1). The 50 open orders (<i>n</i> &equiv; 1, 2, 3 mod 4) become valuable when physical detector channels, clinical cohorts, or industrial batches lock <i>n</i> so padding dummy runs up to 4<i>k</i> is impossible. Even at a fixed open order, <b>block-circulant witnesses (marked with <code>*</code> in the tables below)</b> carry a major engineering advantage over unstructured search grids: because their core consists of cyclic shifts of short seed vectors, an optical spectrometer can step a <b>single etched slit-strip</b> instead of swapping <i>n</i> separate masks, while an FPGA or DSP stores <i>O</i>(<i>n</i>) bits instead of <i>n</i><sup>2</sup> and decodes signals in <i>O</i>(<i>n</i> log <i>n</i>) time via FFT.

</td>
</tr>
</table>
</details>

### Normalization & Quick Verification

Every $n \times n$ sign matrix $X \in \{-1, +1\}^{n \times n}$ has $|\det X|$ divisible by $2^{n-1}$. Following [OEIS A003432](https://oeis.org/A003432) and Orrick's catalog, we report the exact normalized integer determinant

$$D_{\text{norm}}(n) = \frac{|\det X|}{2^{n-1}} \in \mathbb{Z}_{>0}, \qquad R(n) = \frac{|\det X|}{B_{\text{class}}(n)} \le 1,$$

which equals the maximal determinant of an $(n-1) \times (n-1)$ binary $\{0, 1\}$ matrix ([OEIS A003433](https://oeis.org/A003433) tracks raw $|\det X| = 2^{n-1} D_{\text{norm}}(n)$). With our explicit $89 \times 89$ construction (`n = 89` previously lacked an archived matrix), **all 50 open orders (23 <= n <= 119) now have explicit, machine-verified `{-1, +1}` witness matrices** in [`matrices/`](matrices/).

```bash
python3 verify.py --sympy                 # Verifies all 8 new we* records (Bareiss + SymPy), all 50 open witnesses & SHA256SUMS
python3 verify.py --check my_matrix.txt   # Evaluates any external +/-1 or 0/1 candidate matrix against catalog records
python3 verify.py --corrupt               # Adversarial 1-entry mutation self-test (must print FAIL and exit 1)
```

- **Upstream V3R1T4S Verification Package:** [`elder-plinius/V3R1T4S` PR #1053](https://github.com/elder-plinius/V3R1T4S/pull/1053)
- **Companion Open-Science Repository:** <https://github.com/mizopl/we-star>
- **Machine-Readable Tables (full unrounded integers & BibTeX):** [`catalog/hadamard_maxdet_n1_to_120.csv`](catalog/hadamard_maxdet_n1_to_120.csv) · [`catalog/hadamard_maxdet_n1_to_120.json`](catalog/hadamard_maxdet_n1_to_120.json) · [`catalog/new_records.json`](catalog/new_records.json) · [`catalog/references.bib`](catalog/references.bib)
- **For AI Agents & Crawlers:** [`AGENTS.md`](AGENTS.md) · [`llms.txt`](llms.txt)

---

## 1. New Constructive Records (`we*`, Nov 2026 — 8 Orders Improved)

| n | mod 4 | Previous D<sub>norm</sub> | New D<sub>norm</sub> (`we*` 2026) | Gain | R<sub>prev</sub> → R<sub>new</sub> | Previous Holder | Witness |
| ---: | :---: | :--- | :--- | ---: | :---: | :--- | :---: |
| **81** | 1 | 196 · 20<sup>39</sup> | **200 · 20<sup>39</sup>** | **+2.04%** | 0.7723 → **0.7881** | N. Farmakis & S. Kounias (1987) | [n081.txt](matrices/we_star_2026/n081.txt) |
| **89** | 1 | 220 · 22<sup>43</sup> | **226 · 22<sup>43</sup>** | **+2.73%** | 0.7516 → **0.7721** | N. Farmakis & S. Kounias (1987)† | [n089.txt](matrices/we_star_2026/n089.txt) |
| **93** | 1 | 239 · 23<sup>45</sup> | **243 · 23<sup>45</sup>** | **+1.67%** | 0.7640 → **0.7768** | W. P. Orrick & B. Solomon (2003) | [n093.txt](matrices/we_star_2026/n093.txt) |
| **97** | 1 | 256 · 24<sup>47</sup> | **264 · 24<sup>47</sup>** | **+3.12%** | 0.7678 → **0.7918** | B. Solomon & W. P. Orrick (2003) | [n097.txt](matrices/we_star_2026/n097.txt) |
| **101** | 1 | 7110 · 25<sup>48</sup> | **7626 · 25<sup>48</sup>** | **+7.26%** | 0.8024 → **0.8606** | W. P. Orrick & B. Solomon (2003) | [n101.txt](matrices/we_star_2026/n101.txt) |
| **105** | 1 | 3.8474 × 10<sup>74</sup> | **7413 · 26<sup>50</sup>** | **+8.02%** | 0.7022 → **0.7585** | T. Rokicki (2005) | [n105.txt](matrices/we_star_2026/n105.txt) |
| **109** | 1 | 1.8775 × 10<sup>78</sup> | **309 · 27<sup>53</sup>** | **+19.86%** | 0.6482 → **0.7769** | T. Rokicki (2005) | [n109.txt](matrices/we_star_2026/n109.txt) |
| **117** | 1 | 7.2315 × 10<sup>85</sup> | **335 · 29<sup>57</sup>** | **+5.32%** | 0.7186 → **0.7568** | T. Rokicki (2005) | [n117.txt](matrices/we_star_2026/n117.txt) |

<sub>† On <b>n = 89</b>, the value 220 · 22<sup>43</sup> was derived by Farmakis & Kounias (1987) from an order-88 max-excess Hadamard construction without an archived explicit 89 × 89 matrix in <code>MaxDet</code>; <code>matrices/we_star_2026/n089.txt</code> provides an explicit {-1, +1} witness achieving 226 · 22<sup>43</sup>.</sub>

<details>
<summary><strong>Construction methods & exact unrounded integers for the 8 new records</strong></summary>

#### Algebraic construction summary
- **n = 81** (+2.0408%, R: 0.772348 → 0.788110): Closed-Four Hadamard switching on H_80 + [OS-1] triplet border (`Switch_Hcat_s18_L1=640_trip=(40, 56, 59)`). Formula: `Closed-Four Hadamard switching walk on catalog Hadamard core H_80 -> inequivalent H_80 + [OS-1] high-excess 3-row border + 2-column Schur-4 zero-excess alignment`.
- **n = 89** (+2.7273%, R: 0.751646 → 0.772145): Closed-Four switching on Kronecker core H_44 x H_2 + [OS-1] border (`Switch_Kron_44x2_s13_L1=728_n0=12_trip=(42, 53, 58)`). Formula: `Closed-Four Hadamard switching walk on Kronecker Paley-43 x H_2 -> inequivalent H_88 + [OS-1] high-excess 3-row border + 2-column Schur-4 zero-excess alignment`.
- **n = 93** (+1.6736%, R: 0.763984 → 0.776770): Asymmetric Goethals-Seidel core (A^T != A) + [OS-1] R-chirality border (`OS1_M11_AsymGS_23_9_3_1_1_q5_d-1_L1=788_n0=3`). Formula: `Asymmetric Goethals-Seidel Hadamard core H_92 (4x23 circulant blocks, row sums [9, 3, 1, -1]) + [OS-1] 3-row excess border on triplet (65, 85, 91)`.
- **n = 97** (+3.1250%, R: 0.767803 → 0.791797): Asymmetric Goethals-Seidel core (A^T != A) + [OS-1] R-chirality border (`OS1_M11_AsymGS_24_8_4_4_0_q261_d+1_L1=864_n0=6`). Formula: `Asymmetric Goethals-Seidel Hadamard core H_96 (4x24 circulant blocks, row sums [8, 4, 4, 0]) + [OS-1] 3-row excess border on triplet (3, 11, 19) + Schur-4 border polish (3 border flips)`.
- **n = 101** (+7.2574%, R: 0.802401 → 0.860635): Coupled two-sided border flip on doubly 3-normalized [OS-2] Goethals-Seidel core (`SMv2_Orrick_n101_s0`). Formula: `Doubly 3-normalized [OS-2] Hadamard core H_100 + coupled rank-2 border/core unfreezing (6 sign flips from [OS2007] breaking the frozen core orthogonality barrier)`.
- **n = 105** (+8.0164%, R: 0.702238 → 0.758532): Asymmetric Goethals-Seidel core (A^T != A) + [OS-1] R-chirality border (`WarmSM_OS1_M14_AsymGS_26_8_6_2_0_q0_sl0_t100_d+1_L1=920_n0=16_r1`). Formula: `Asymmetric Goethals-Seidel Hadamard core H_104 (4x26 circulant blocks, row sums [-6, -4, 4, 6]) + [OS-1] 3-row excess border on triplet (12, 57, 98) + rank-1 Sherman-Morrison unfreezing (11 cell flips)`.
- **n = 109** (+19.8565%, R: 0.648191 → 0.776900): Asymmetric Goethals-Seidel core (A^T != A) + [OS-1] R-chirality border (`OS1_M11_AsymGS_27_9_3_3_3_q49_d+1_L1=1020_n0=3`). Formula: `Asymmetric Goethals-Seidel Hadamard core H_108 (4x27 circulant blocks, row sums [9, 3, 3, 3]) + [OS-1] 3-row excess border on triplet (83, 92, 101)`.
- **n = 117** (+5.3179%, R: 0.718566 → 0.756779): Asymmetric Goethals-Seidel core (A^T != A) + [OS-1] R-chirality border (`OS1_M14_AsymGS_29_7_7_3_3_q8_sl0_t100_d-1_L1=1108_n0=3`). Formula: `Asymmetric Goethals-Seidel Hadamard core H_116 (4x29 circulant blocks, row sums [-7, -7, -3, 3]) + [OS-1] 3-row excess border on triplet (6, 66, 98)`.

#### Exact integer `D_norm = |det X| / 2^(n-1)` values
```text
n =  81 | D_norm_prev = 107752139522048000000000000000000000000000000000000000
        | D_norm_new  = 109951162777600000000000000000000000000000000000000000
n =  89 | D_norm_prev = 1165729954414365496202893614947911393918604879059221018050560
        | D_norm_new  = 1197522589534757282462972531719218068298203193942654318542848
n =  93 | D_norm_prev = 4530545764797772467125412055473257402652070190439124366323586777
        | D_norm_new  = 4606370798518237278290690918326366313156707348438105527266240949
n =  97 | D_norm_prev = 18974411972804596157953446332065856959415742925610486623450687340544
        | D_norm_new  = 19567362346954739787889491529942914989397484892035814330433521319936
n = 101 | D_norm_prev = 89740816577942302983638241586569761487623964058002457022666931152343750
        | D_norm_new  = 96253652211446976449117472621544444599805956386262550950050354003906250
n = 105 | D_norm_prev = 384743747284587846942188483608583600735563257274079045682562818122226665625
        | D_norm_new  = 415586468672646649836878426545094933992421758510370836431909484174705164288
n = 109 | D_norm_prev = 1877483841834222713683125485573890072927341503771574105042244268774869908324352
        | D_norm_new  = 2250287355442438900295757619807994219684406160806013443723720866712282032585703
n = 117 | D_norm_prev = 72314994254613394470813334860641212637154205430658679584420251275641935854599391452615
        | D_norm_new  = 76160657954499973539057578092053129947262150826797050046488183822315692905877125578515
```
</details>

<details>
<summary><strong>Subsection 1B — Hadamard-Inequivalent Matrices Attaining Existing Open Records (19 non-isomorphic classes, incl. n = 51)</strong></summary>

<sub>These constructions attain the exact best-known determinant $D_{\text{norm}} = D_{\text{rec}}$ (and are therefore <b>not</b> counted among the 8 strict determinant improvements above), while establishing provably distinct Hadamard equivalence classes under signed row/column permutations and matrix transposition ($P_1 Y P_2$ and $P_1 Y^T P_2$). Full certificates in <a href="catalog/tied_constructions.json"><code>catalog/tied_constructions.json</code></a>.</sub>

| n | Tied Record Holder | D<sub>norm</sub> | Witness | Invariant Proof (vs Y &amp; Y<sup>T</sup>) | Constructive Formula |
| ---: | :--- | ---: | :---: | :--- | :--- |
| **29** | B. Solomon (2002) | 320 · 7<sup>12</sup> | [tied_n029.txt](matrices/we_star_2026/tied_n029.txt) | `4-profile mod T (5efebc67 != d96603a3)` | Algebraic PaleyII_13 Hadamard core H_28 + doubly 3-normalized [OS-2] border on R=(0, 2, 4... |
| **34** | W. P. Orrick (2005) | 256 · 8<sup>15</sup> | [tied_n034.txt](matrices/we_star_2026/tied_n034.txt) | `4-profile mod T (e1454ea9 != 84bae9ac)` | Constructed 2-bordered 2-circulant ([[H_2, J, K], [J^T, A, B], [K^T, B^T, -A^T]]) (block ... |
| **49** | B. Solomon & W. P. Orrick (2003) | 96 · 12<sup>23</sup> | [tied_n049.txt](matrices/we_star_2026/tied_n049.txt) | `2-point Gram (cbdf9764 != b2180b74)` | Switched/algebraic DisjointClosedFour_Hcat_s0 Hadamard core H_48 + [OS-1] border on tripl... |
| **51** | Butbaia, Gukov et al. (2026) | 1.7776 × 10<sup>28</sup> | [tied_n051.txt](matrices/we_star_2026/tied_n051.txt) | `4-profile mod T (7fa6a09f != f63cd8da)` | Cold-constructed 1-bordered symmetric-skew 2-circulant Type 2 ([[A, B, -j], [B, A^T, j], ... |
| **53** | Orrick, Solomon et al. (2003) | 105 · 13<sup>25</sup> | [tied_n053.txt](matrices/we_star_2026/tied_n053.txt) | `2-point Gram (39dc1c4e != 9fa9751b)` | Cold-constructed Asymmetric Goethals-Seidel Hadamard core H_52 (4x13 blocks) + [OS-1] bor... |
| **55** | T. Rokicki; Kazmenko & Trofimov (2005) | 2.2572 × 10<sup>31</sup> | [tied_n055.txt](matrices/we_star_2026/tied_n055.txt) | `4-profile mod T (498f9db7 != 897c2c09)` | Cold-constructed 1-bordered 2-circulant Type 1 ([[+1, j^T, -j^T], [j, A, B], [-j, B^T, -A... |
| **58** | T. Rokicki; Kazmenko & Trofimov (2005) | 6.5774 × 10<sup>33</sup> | [tied_n058.txt](matrices/we_star_2026/tied_n058.txt) | `4-profile mod T (286ef9da != bc7ddaea)` | Cold-constructed unbordered 2-circulant [[A, B], [-B^T, A^T]] of block size m=29 via scal... |
| **65** | W. P. Orrick & B. Solomon (2003) | 148 · 16<sup>31</sup> | [tied_n065.txt](matrices/we_star_2026/tied_n065.txt) | `2-point Gram (a86dd032 != 66ed6048)` | Switched/algebraic DisjointClosedFour_Hcat_s0 Hadamard core H_64 + [OS-1] border on tripl... |
| **67** | T. Rokicki (2005) | 1.1400 × 10<sup>41</sup> | [tied_n067.txt](matrices/we_star_2026/tied_n067.txt) | `4-profile mod T (8c0b4faf != ad98b8b8)` | Cold-constructed 1-bordered 2-circulant Type 1 ([[+1, j^T, -j^T], [j, A, B], [-j, B^T, -A... |
| **69** | Orrick, Solomon et al. (2003) | 155 · 17<sup>33</sup> | [tied_n069.txt](matrices/we_star_2026/tied_n069.txt) | `2-point Gram (96a76a74 != 73d9af27)` | Cold-constructed Asymmetric Goethals-Seidel Hadamard core H_68 (4x17 blocks) + [OS-1] bor... |
| **70** | W. P. Orrick & T. Rokicki (2005) | 1156 · 17<sup>33</sup> | [tied_n070.txt](matrices/we_star_2026/tied_n070.txt) | `4-profile mod T (f4c9643a != 81fe1930)` | Constructed 2-bordered 2-circulant ([[H_2, J, K], [J^T, A, B], [K^T, B^T, -A^T]]) (block ... |
| **71** | T. Rokicki (2005) | 31,186,944 · 17<sup>30</sup> | [tied_n071.txt](matrices/we_star_2026/tied_n071.txt) | `4-profile mod T (7ba4d457 != 88f85663)` | Constructed 1-bordered 2-circulant Type 3 ([[A, B, -j], [B^T, -A^T, j], [-j^T, -j^T, -1]]... |
| **73** | B. Solomon & W. P. Orrick (2003) | 174 · 18<sup>35</sup> | [tied_n073.txt](matrices/we_star_2026/tied_n073.txt) | `2-point Gram (7fdaf7b1 != 1786d082)` | Switched/algebraic DisjointClosedFour_Hcat_s0 Hadamard core H_72 + [OS-1] border on tripl... |
| **75** | T. Rokicki (2005) | 6.1512 × 10<sup>47</sup> | [tied_n075.txt](matrices/we_star_2026/tied_n075.txt) | `4-profile mod T (f0f13f50 != aa8b36e8)` | Cold-constructed 1-bordered 2-circulant Type 1 ([[+1, j^T, -j^T], [j, A, B], [-j, B^T, -A... |
| **77** | B. Solomon & W. P. Orrick (2003) | 183 · 19<sup>37</sup> | [tied_n077.txt](matrices/we_star_2026/tied_n077.txt) | `2-point Gram (9de1db65 != 51bd9df8)` | Cold-constructed Asymmetric Goethals-Seidel Hadamard core H_76 (4x19 blocks) + [OS-1] bor... |
| **79** | T. Rokicki (2005) | 1.6796 × 10<sup>51</sup> | [tied_n079.txt](matrices/we_star_2026/tied_n079.txt) | `4-profile mod T (56be81f7 != 163d0aa0)` | Cold-constructed 1-bordered 2-circulant Type 1 ([[+1, j^T, -j^T], [j, A, B], [-j, B^T, -A... |
| **83** | T. Rokicki (2005) | 5.0710 × 10<sup>54</sup> | [tied_n083.txt](matrices/we_star_2026/tied_n083.txt) | `4-profile mod T (b2f1ed80 != 8ffe7090)` | Constructed 1-bordered 2-circulant Type 3 ([[A, B, -j], [B^T, -A^T, j], [-j^T, -j^T, -1]]... |
| **87** | T. Rokicki (2005) | 1.6851 × 10<sup>58</sup> | [tied_n087.txt](matrices/we_star_2026/tied_n087.txt) | `4-profile mod T (0dcf6624 != 0c271884)` | Cold-constructed 1-bordered 2-circulant Type 0 ([[-1, j^T, -j^T], [j, A, B], [-j, B^T, -A... |
| **106** | W. P. Orrick (2005) | 2704 · 26<sup>51</sup> | [tied_n106.txt](matrices/we_star_2026/tied_n106.txt) | `2-point Gram (8fc5de8f != a9a15cca)` | Exact 1-cell neutral plateau unfreezing breaking 2-circulant symmetry of 2-bordered 2-cir... |
</details>

<details>
<summary><strong>Independent Algebraic Cold-Reproductions of Existing Open Records (8 witnesses)</strong></summary>

- **n = 78** (confirms T. Rokicki (2005)): [tied_n078.txt](matrices/we_star_2026/tied_n078.txt) — Transposed-chirality 2-circulant (A, B^T) realization of 2-bordered 2-circulant ([[H_2, J, K], [J^T, A...
- **n = 94** (confirms T. Rokicki (2005)): [tied_n094.txt](matrices/we_star_2026/tied_n094.txt) — Transposed-chirality 2-circulant (A, B^T) realization of 2-bordered 2-circulant ([[H_2, J, K], [J^T, A...
- **n = 95** (confirms T. Rokicki (2005)): [tied_n095.txt](matrices/we_star_2026/tied_n095.txt) — Transposed-chirality 2-circulant (A, B^T) realization of 1-bordered 2-circulant Type 0 ([[-1, j^T, -j^...
- **n = 103** (confirms Butbaia, Gukov et al. (Caltech) (2026)): [tied_n103.txt](matrices/we_star_2026/tied_n103.txt) — Transposed-chirality 2-circulant (A, B^T) realization of 1-bordered 2-circulant Type 1 ([[+1, j^T, -j^...
- **n = 107** (confirms Butbaia, Gukov et al. (Caltech) (2026)): [tied_n107.txt](matrices/we_star_2026/tied_n107.txt) — Transposed-chirality 2-circulant (A, B^T) realization of 1-bordered 2-circulant Type 1 ([[+1, j^T, -j^...
- **n = 111** (confirms Butbaia, Gukov et al. (Caltech) (2026)): [tied_n111.txt](matrices/we_star_2026/tied_n111.txt) — Transposed-chirality 2-circulant (A, B^T) realization of 1-bordered 2-circulant Type 3 ([[A, B, -j], [...
- **n = 115** (confirms Butbaia, Gukov et al. (Caltech) (2026)): [tied_n115.txt](matrices/we_star_2026/tied_n115.txt) — Transposed-chirality 2-circulant (A, B^T) realization of 1-bordered 2-circulant Type 1 ([[+1, j^T, -j^...
- **n = 119** (confirms Butbaia, Gukov et al. (Caltech) (2026)): [tied_n119.txt](matrices/we_star_2026/tied_n119.txt) — Transposed-chirality 2-circulant (A, B^T) realization of 1-bordered 2-circulant Type 3 ([[A, B, -j], [...
</details>

---

## 2. Residue-Class Summary (1 <= n <= 120)

| Residue Class | Orders | Optimal | Open | Improved (`we*` 2026) | Upper Bound $B_{\text{class}}(n)$ | Status Summary |
| :--- | ---: | ---: | ---: | ---: | :--- | :--- |
| **0 mod 4**\* | 30 | 30 | 0 | **0** | Hadamard: $n^{n/2}$ | All 30 multiples of 4 achieve $R = 1$ (`*` dyadic Walsh $2^m$: 4, 8, 16, 32, 64) |
| **1 mod 4** | 30 | 11 | 19 | **8** | Barba: $\sqrt{2n-1}\,(n-1)^{(n-1)/2}$ | 11 optimal (1, 5, 9, 13, 17, 21, 25, 37, 41, 61, 113); 19 open |
| **2 mod 4** | 30 | 24 | 6 | **0** | Ehlich–Wojtas: $2(n-1)(n-2)^{(n-2)/2}$ | 24 optimal (2*, 22 [Chasiotis et al. 2018], 22 EW-tight); 6 open |
| **3 mod 4** | 30 | 5 | 25 | **0** | Ehlich $s$-block Gram ($s \in \{3,5,6,7\}$) | 5 optimal (3, 7, 11, 15, 19); 25 open |

---

## 3. Best Known Lower Bounds on All 50 Open Orders (23 <= n <= 119)

Grouped by residue class $n \bmod 4$ so that $R(n) = |\det X| / B_{\text{class}}(n)$ is compared against a single upper-bound formula within each table (all 50 rows verified from linked `{-1, +1}` matrices by [`verify.py`](verify.py); witnesses marked **`*`** have a verified block-circulant core):

### 3.1 Orders n ≡ 1 (mod 4) — Barba Regime (19 Open Orders)

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
| **81** | **1.0995 × 10<sup>53</sup>** | **200 · 20<sup>39</sup>** | **0.7881** | **we\* (@mizopl)** (+2.04%) | **2026** | [n081.txt](matrices/we_star_2026/n081.txt) |
| **85** | 3.4902 × 10<sup>56</sup> | 4509 · 21<sup>40</sup> | 0.7865 | B. Solomon & W. P. Orrick | 2003 | [n085.txt](matrices/literature_witnesses/n085.txt) |
| **89** | **1.1975 × 10<sup>60</sup>** | **226 · 22<sup>43</sup>** | **0.7721** | **we\* (@mizopl)** (+2.73%) | **2026** | [n089.txt](matrices/we_star_2026/n089.txt) |
| **93** | **4.6064 × 10<sup>63</sup>** | **243 · 23<sup>45</sup>** | **0.7768** | **we\* (@mizopl)** (+1.67%) | **2026** | [n093.txt](matrices/we_star_2026/n093.txt) |
| **97** | **1.9567 × 10<sup>67</sup>** | **264 · 24<sup>47</sup>** | **0.7918** | **we\* (@mizopl)** (+3.12%) | **2026** | [n097.txt](matrices/we_star_2026/n097.txt) |
| **101** | **9.6254 × 10<sup>70</sup>** | **7626 · 25<sup>48</sup>** | **0.8606** | **we\* (@mizopl)** (+7.26%) | **2026** | [n101.txt](matrices/we_star_2026/n101.txt) |
| **105** | **4.1559 × 10<sup>74</sup>** | **7413 · 26<sup>50</sup>** | **0.7585** | **we\* (@mizopl)** (+8.02%) | **2026** | [n105.txt](matrices/we_star_2026/n105.txt) |
| **109** | **2.2503 × 10<sup>78</sup>** | **309 · 27<sup>53</sup>** | **0.7769** | **we\* (@mizopl)** (+19.86%) | **2026** | [n109.txt](matrices/we_star_2026/n109.txt) |
| **117** | **7.6161 × 10<sup>85</sup>** | **335 · 29<sup>57</sup>** | **0.7568** | **we\* (@mizopl)** (+5.32%) | **2026** | [n117.txt](matrices/we_star_2026/n117.txt) |

### 3.2 Orders n ≡ 2 (mod 4) — Ehlich–Wojtas Regime (6 Open Orders)

*Upper bound: $B_2(n) = 2(n - 1)(n - 2)^{(n - 2)/2}$. All other 24 orders $n \equiv 2 \pmod 4$ up to $118$ are proven optimal.*

| n | D<sub>norm</sub> | Factored Form | R = \|det\|/B | Record Holder | Year | Witness |
| ---: | ---: | :--- | ---: | :--- | :---: | :---: |
| **34** | 9.0072 × 10<sup>15</sup> | 256 · 8<sup>15</sup> | 0.9697 | W. P. Orrick | 2005 | [n034.txt](matrices/literature_witnesses/n034.txt)* |
| **58** | 6.5774 × 10<sup>33</sup> | — | 0.9345 | T. Rokicki; Kazmenko & Trofimov | 2005 | [n058.txt](matrices/literature_witnesses/n058.txt)* |
| **70** | 4.6534 × 10<sup>43</sup> | 1156 · 17<sup>33</sup> | 0.9855 | W. P. Orrick & T. Rokicki | 2005 | [n070.txt](matrices/literature_witnesses/n070.txt)* |
| **78** | 2.9042 × 10<sup>50</sup> | — | 0.9636 | T. Rokicki | 2005 | [n078.txt](matrices/literature_witnesses/n078.txt)* |
| **94** | 3.9312 × 10<sup>64</sup> | — | 0.9695 | T. Rokicki | 2005 | [n094.txt](matrices/literature_witnesses/n094.txt)* |
| **106** | 3.9414 × 10<sup>75</sup> | 2704 · 26<sup>51</sup> | 0.9905 | W. P. Orrick | 2005 | [n106.txt](matrices/literature_witnesses/n106.txt)* |

### 3.3 Orders n ≡ 3 (mod 4) — Ehlich Block-Gram Regime (25 Open Orders)

*Upper bound: Ehlich $s$-block Gram bound ($s = 6$ for $23 \le n \le 59$, $s = 7$ for $63 \le n \le 119$).*

| n | D<sub>norm</sub> | Factored Form | R = \|det\|/B | Record Holder | Year | Witness |
| ---: | ---: | :--- | ---: | :--- | :---: | :---: |
| **23** | 6.6267 × 10<sup>8</sup> | 662,671,875 | 0.9320 | Dowdeswell, Neubauer et al. | 2004 | [n023.txt](matrices/literature_witnesses/n023.txt) |
| **27** | 1.9809 × 10<sup>11</sup> | 546 · 6<sup>11</sup> | 0.9673 | H. Tamura | 2005 | [n027.txt](matrices/literature_witnesses/n027.txt) |
| **31** | 7.5961 × 10<sup>13</sup> | 784 · 7<sup>13</sup> | 0.9561 | H. Tamura | 2005 | [n031.txt](matrices/literature_witnesses/n031.txt) |
| **35** | 3.7436 × 10<sup>16</sup> | 1064 · 8<sup>15</sup> | 0.9379 | H. Tamura | 2011 | [n035.txt](matrices/literature_witnesses/n035.txt) |
| **39** | 2.4015 × 10<sup>19</sup> | 1440 · 9<sup>17</sup> | 0.9527 | H. Tamura | 2011 | [n039.txt](matrices/literature_witnesses/n039.txt) |
| **43** | 1.8900 × 10<sup>22</sup> | 1890 · 10<sup>19</sup> | 0.9663 | H. Tamura | 2011 | [n043.txt](matrices/literature_witnesses/n043.txt) |
| **47** | 1.6888 × 10<sup>25</sup> | 3,037,500 · 11<sup>18</sup> | 0.9239 | Orrick, Solomon et al. | 2003 | [n047.txt](matrices/literature_witnesses/n047.txt) |
| **51** | 1.7776 × 10<sup>28</sup> | — | 0.8776 | **Butbaia, Gukov et al.** (+3.11%) | **2026** | [n051.txt](matrices/literature_witnesses/n051.txt)* |
| **55** | 2.2572 × 10<sup>31</sup> | — | 0.8593 | T. Rokicki; Kazmenko & Trofimov | 2005 | [n055.txt](matrices/literature_witnesses/n055.txt) |
| **59** | 3.8031 × 10<sup>34</sup> | 4312 · 14<sup>27</sup> | 0.9653 | H. Tamura | 2005 | [n059.txt](matrices/literature_witnesses/n059.txt) |
| **63** | 6.5387 × 10<sup>37</sup> | 5115 · 15<sup>29</sup> | 0.9662 | H. Tamura | 2011 | [n063.txt](matrices/literature_witnesses/n063.txt) |
| **67** | 1.1400 × 10<sup>41</sup> | — | 0.8635 | T. Rokicki | 2005 | [n067.txt](matrices/literature_witnesses/n067.txt)* |
| **71** | 2.5553 × 10<sup>44</sup> | 31,186,944 · 17<sup>30</sup> | 0.8804 | T. Rokicki | 2005 | [n071.txt](matrices/literature_witnesses/n071.txt) |
| **75** | 6.1512 × 10<sup>47</sup> | — | 0.8613 | T. Rokicki | 2005 | [n075.txt](matrices/literature_witnesses/n075.txt)* |
| **79** | 1.6796 × 10<sup>51</sup> | — | 0.8591 | T. Rokicki | 2005 | [n079.txt](matrices/literature_witnesses/n079.txt) |
| **83** | 5.0710 × 10<sup>54</sup> | — | 0.8561 | T. Rokicki | 2005 | [n083.txt](matrices/literature_witnesses/n083.txt)* |
| **87** | 1.6851 × 10<sup>58</sup> | — | 0.8527 | T. Rokicki | 2005 | [n087.txt](matrices/literature_witnesses/n087.txt)* |
| **91** | 6.1448 × 10<sup>61</sup> | — | 0.8501 | T. Rokicki | 2005 | [n091.txt](matrices/literature_witnesses/n091.txt) |
| **95** | 2.4387 × 10<sup>65</sup> | — | 0.8447 | T. Rokicki | 2005 | [n095.txt](matrices/literature_witnesses/n095.txt)* |
| **99** | 1.0657 × 10<sup>69</sup> | — | 0.8496 | T. Rokicki | 2005 | [n099.txt](matrices/literature_witnesses/n099.txt) |
| **103** | 4.9829 × 10<sup>72</sup> | — | 0.8433 | **Butbaia, Gukov et al.** (+0.00004%) | **2026** | [n103.txt](matrices/literature_witnesses/n103.txt)* |
| **107** | 2.5405 × 10<sup>76</sup> | — | 0.8445 | **Butbaia, Gukov et al.** (+0.44%) | **2026** | [n107.txt](matrices/literature_witnesses/n107.txt)* |
| **111** | 1.3978 × 10<sup>80</sup> | — | 0.8469 | **Butbaia, Gukov et al.** (+1.26%) | **2026** | [n111.txt](matrices/literature_witnesses/n111.txt)* |
| **115** | 8.2608 × 10<sup>83</sup> | — | 0.8488 | **Butbaia, Gukov et al.** (+1.83%) | **2026** | [n115.txt](matrices/literature_witnesses/n115.txt)* |
| **119** | 5.2308 × 10<sup>87</sup> | — | 0.8503 | **Butbaia, Gukov et al.** (+2.12%) | **2026** | [n119.txt](matrices/literature_witnesses/n119.txt)* |

<br>

<details>
<summary><strong>Complete 120-Order Master Table (1 <= n <= 120, Proven Optimal + Open)</strong></summary>

<sub><code>*</code> on <code>n</code> marks dyadic Sylvester / Walsh–Hadamard orders (<i>n</i> = 2<sup><i>m</i></sup>, fast <i>O</i>(<i>n</i> log <i>n</i>) butterfly); <code>*</code> on <code>Witness</code> marks open-order matrices with a verified block-circulant core.</sub>

| n | mod 4 | Status | D<sub>norm</sub> | Factored Form | R = \|det\|/B | Attribution | Year | Ref | Witness |
| ---: | :---: | :---: | ---: | :--- | ---: | :--- | :---: | :---: | :---: |
| 1* | 1 | OPT | 1 | 1 | 1.0000 | Sylvester | 1867 | Syl1867 | — |
| 2* | 2 | OPT | 1 | 1 | 1.0000 | Sylvester | 1867 | Syl1867 | — |
| 3 | 3 | OPT | 1 | 1 | 1.0000 | Sylvester; Williamson | 1867 | Syl1867 | — |
| 4* | 0 | OPT | 2 | 2 | 1.0000 | Sylvester (Hadamard) | 1867 | Syl1867 | — |
| 5 | 1 | OPT | 3 | 3 | 1.0000 | Mood (Barba bound) | 1946 | Moo1946 | — |
| 6 | 2 | OPT | 5 | 5 | 1.0000 | Ehlich (EW bound) | 1964 | Ehl1964a | — |
| 7 | 3 | OPT | 9 | 9 | 0.9820 | Williamson | 1946 | Wil1946 | — |
| 8* | 0 | OPT | 32 | 32 | 1.0000 | Sylvester (Hadamard) | 1867 | Syl1867 | — |
| 9 | 1 | OPT | 56 | 56 | 0.8489 | Ehlich | 1964 | Ehl1964a | — |
| 10 | 2 | OPT | 144 | 144 | 1.0000 | Ehlich (EW bound) | 1964 | Ehl1964a | — |
| 11 | 3 | OPT | 320 | 320 | 0.9449 | Ehlich | 1964 | Ehl1964b | — |
| 12 | 0 | OPT | 1,458 | 1,458 | 1.0000 | Hadamard | 1893 | Had1893 | — |
| 13 | 1 | OPT | 3,645 | 3,645 | 1.0000 | Raghavarao (Barba bound) | 1959 | Rag1959 | — |
| 14 | 2 | OPT | 9,477 | 9,477 | 1.0000 | Ehlich (EW bound) | 1964 | Ehl1964a | — |
| 15 | 3 | OPT | 25,515 | 25,515 | 0.9707 | Ehlich & Zeller; opt. Orrick | 2005 | Orr2005b | — |
| 16* | 0 | OPT | 131,072 | 131,072 | 1.0000 | Sylvester (Hadamard) | 1867 | Syl1867 | — |
| 17 | 1 | OPT | 327,680 | 327,680 | 0.8704 | Moyssiadis & Kounias | 1982 | MK1982 | — |
| 18 | 2 | OPT | 1.1141 × 10<sup>6</sup> | 1,114,112 | 1.0000 | Ehlich (EW bound) | 1964 | Ehl1964a | — |
| 19 | 3 | OPT | 3.4120 × 10<sup>6</sup> | 3,411,968 | 0.9750 | Cohn (1989); opt. Brent et al. | 2011 | BOOZ2011 | — |
| 20 | 0 | OPT | 1.9531 × 10<sup>7</sup> | 19,531,250 | 1.0000 | Hadamard | 1893 | Had1893 | — |
| 21 | 1 | OPT | 5.6641 × 10<sup>7</sup> | 56,640,625 | 0.9058 | Chadjipantelis, Kounias et al. | 1987 | CKM1987 | — |
| 22 | 2 | OPT | 1.9531 × 10<sup>8</sup> | 195,312,500 | 0.9524 | Chasiotis, Kounias & Farmakis | 2018 | CKF2018 | — |
| 23 | 3 | OPEN | 6.6267 × 10<sup>8</sup> | 662,671,875 | 0.9320 | Dowdeswell, Neubauer et al. | 2004 | DNST2004 | [txt](matrices/literature_witnesses/n023.txt) |
| 24 | 0 | OPT | 4.3536 × 10<sup>9</sup> | 12 · 6<sup>11</sup> | 1.0000 | Paley / Williamson (Hadamard) | 1933 | Pal1933 | — |
| 25 | 1 | OPT | 1.5237 × 10<sup>10</sup> | 42 · 6<sup>11</sup> | 1.0000 | Brouwer (Barba bound) | 1983 | Bro1983 | — |
| 26 | 2 | OPT | 5.4420 × 10<sup>10</sup> | 150 · 6<sup>11</sup> | 1.0000 | Ehlich (EW bound) | 1964 | Ehl1964a | — |
| 27 | 3 | OPEN | 1.9809 × 10<sup>11</sup> | 546 · 6<sup>11</sup> | 0.9673 | H. Tamura | 2005 | Tam2005 | [txt](matrices/literature_witnesses/n027.txt) |
| 28 | 0 | OPT | 1.3564 × 10<sup>12</sup> | 14 · 7<sup>13</sup> | 1.0000 | Paley / Williamson (Hadamard) | 1933 | Pal1933 | — |
| 29 | 1 | OPEN | 4.4292 × 10<sup>12</sup> | 320 · 7<sup>12</sup> | 0.8650 | B. Solomon | 2002 | Sol2002 | [txt](matrices/literature_witnesses/n029.txt) |
| 30 | 2 | OPT | 1.9668 × 10<sup>13</sup> | 203 · 7<sup>13</sup> | 1.0000 | Ehlich (EW bound) | 1964 | Ehl1964a | — |
| 31 | 3 | OPEN | 7.5961 × 10<sup>13</sup> | 784 · 7<sup>13</sup> | 0.9561 | H. Tamura | 2005 | Tam2005 | [txt](matrices/literature_witnesses/n031.txt) |
| 32* | 0 | OPT | 5.6295 × 10<sup>14</sup> | 16 · 8<sup>15</sup> | 1.0000 | Sylvester (Hadamard) | 1867 | Syl1867 | — |
| 33 | 1 | OPEN | 1.9395 × 10<sup>15</sup> | 441 · 8<sup>14</sup> | 0.8547 | B. Solomon | 2002 | Sol2002 | [txt](matrices/literature_witnesses/n033.txt) |
| 34 | 2 | OPEN | 9.0072 × 10<sup>15</sup> | 256 · 8<sup>15</sup> | 0.9697 | W. P. Orrick | 2005 | Orr2005 | [txt](matrices/literature_witnesses/n034.txt)* |
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
| 51 | 3 | OPEN | 1.7776 × 10<sup>28</sup> | — | 0.8776 | Butbaia, Gukov et al. | 2026 | BG2026 | [txt](matrices/literature_witnesses/n051.txt)* |
| 52 | 0 | OPT | 1.8347 × 10<sup>29</sup> | 26 · 13<sup>25</sup> | 1.0000 | Paley / Williamson (Hadamard) | 1933 | Pal1933 | — |
| 53 | 1 | OPEN | 7.4092 × 10<sup>29</sup> | 105 · 13<sup>25</sup> | 0.7882 | Orrick, Solomon et al. | 2003 | OSDS2003 | [txt](matrices/literature_witnesses/n053.txt) |
| 54 | 2 | OPT | 4.8619 × 10<sup>30</sup> | 689 · 13<sup>25</sup> | 1.0000 | Yang (EW bound) | 1968 | Yan1968 | — |
| 55 | 3 | OPEN | 2.2572 × 10<sup>31</sup> | — | 0.8593 | T. Rokicki; Kazmenko & Trofimov | 2005 | Rok2005 | [txt](matrices/literature_witnesses/n055.txt) |
| 56 | 0 | OPT | 2.4695 × 10<sup>32</sup> | 28 · 14<sup>27</sup> | 1.0000 | Paley / Williamson (Hadamard) | 1933 | Pal1933 | — |
| 57 | 1 | OPEN | 1.1730 × 10<sup>33</sup> | 133 · 14<sup>27</sup> | 0.8937 | J.-C. Meyrignac & B. Felgenhauer | 2005 | Mey2005 | [txt](matrices/literature_witnesses/n057.txt) |
| 58 | 2 | OPEN | 6.5774 × 10<sup>33</sup> | — | 0.9345 | T. Rokicki; Kazmenko & Trofimov | 2005 | Rok2005 | [txt](matrices/literature_witnesses/n058.txt)* |
| 59 | 3 | OPEN | 3.8031 × 10<sup>34</sup> | 4312 · 14<sup>27</sup> | 0.9653 | H. Tamura | 2005 | Tam2005 | [txt](matrices/literature_witnesses/n059.txt) |
| 60 | 0 | OPT | 3.8350 × 10<sup>35</sup> | 30 · 15<sup>29</sup> | 1.0000 | Paley / Williamson (Hadamard) | 1933 | Pal1933 | — |
| 61 | 1 | OPT | 2.1093 × 10<sup>36</sup> | 165 · 15<sup>29</sup> | 1.0000 | Brouwer (Barba bound) | 1983 | Bro1983 | — |
| 62 | 2 | OPT | 1.1697 × 10<sup>37</sup> | 915 · 15<sup>29</sup> | 1.0000 | Yang (EW bound) | 1968 | Yan1968 | — |
| 63 | 3 | OPEN | 6.5387 × 10<sup>37</sup> | 5115 · 15<sup>29</sup> | 0.9662 | H. Tamura | 2011 | Tam2011 | [txt](matrices/literature_witnesses/n063.txt) |
| 64* | 0 | OPT | 6.8056 × 10<sup>38</sup> | 32 · 16<sup>31</sup> | 1.0000 | Sylvester (Hadamard) | 1867 | Syl1867 | — |
| 65 | 1 | OPEN | 3.1476 × 10<sup>39</sup> | 148 · 16<sup>31</sup> | 0.8144 | W. P. Orrick & B. Solomon | 2003 | OS2007 | [txt](matrices/literature_witnesses/n065.txt) |
| 66 | 2 | OPT | 2.2118 × 10<sup>40</sup> | 1040 · 16<sup>31</sup> | 1.0000 | Yang (EW bound) | 1968 | Yan1968 | — |
| 67 | 3 | OPEN | 1.1400 × 10<sup>41</sup> | — | 0.8635 | T. Rokicki | 2005 | Rok2005 | [txt](matrices/literature_witnesses/n067.txt)* |
| 68 | 0 | OPT | 1.3687 × 10<sup>42</sup> | 34 · 17<sup>33</sup> | 1.0000 | Paley / Williamson (Hadamard) | 1933 | Pal1933 | — |
| 69 | 1 | OPEN | 6.2394 × 10<sup>42</sup> | 155 · 17<sup>33</sup> | 0.7790 | Orrick, Solomon et al. | 2003 | OSDS2003 | [txt](matrices/literature_witnesses/n069.txt) |
| 70 | 2 | OPEN | 4.6534 × 10<sup>43</sup> | 1156 · 17<sup>33</sup> | 0.9855 | W. P. Orrick & T. Rokicki | 2005 | Orr2005 | [txt](matrices/literature_witnesses/n070.txt)* |
| 71 | 3 | OPEN | 2.5553 × 10<sup>44</sup> | 31,186,944 · 17<sup>30</sup> | 0.8804 | T. Rokicki | 2005 | Rok2005 | [txt](matrices/literature_witnesses/n071.txt) |
| 72 | 0 | OPT | 3.0963 × 10<sup>45</sup> | 36 · 18<sup>35</sup> | 1.0000 | Paley / Williamson (Hadamard) | 1933 | Pal1933 | — |
| 73 | 1 | OPEN | 1.4965 × 10<sup>46</sup> | 174 · 18<sup>35</sup> | 0.8028 | B. Solomon & W. P. Orrick | 2003 | OS2007 | [txt](matrices/literature_witnesses/n073.txt) |
| 74 | 2 | OPT | 1.1301 × 10<sup>47</sup> | 1314 · 18<sup>35</sup> | 1.0000 | Cohn; Koukouvinos, Kounias & Seberry | 1989 | Coh1989 | — |
| 75 | 3 | OPEN | 6.1512 × 10<sup>47</sup> | — | 0.8613 | T. Rokicki | 2005 | Rok2005 | [txt](matrices/literature_witnesses/n075.txt)* |
| 76 | 0 | OPT | 7.8283 × 10<sup>48</sup> | 38 · 19<sup>37</sup> | 1.0000 | Paley / Williamson (Hadamard) | 1933 | Pal1933 | — |
| 77 | 1 | OPEN | 3.7699 × 10<sup>49</sup> | 183 · 19<sup>37</sup> | 0.7787 | B. Solomon & W. P. Orrick | 2003 | OS2007 | [txt](matrices/literature_witnesses/n077.txt) |
| 78 | 2 | OPEN | 2.9042 × 10<sup>50</sup> | — | 0.9636 | T. Rokicki | 2005 | Rok2005 | [txt](matrices/literature_witnesses/n078.txt)* |
| 79 | 3 | OPEN | 1.6796 × 10<sup>51</sup> | — | 0.8591 | T. Rokicki | 2005 | Rok2005 | [txt](matrices/literature_witnesses/n079.txt) |
| 80 | 0 | OPT | 2.1990 × 10<sup>52</sup> | 40 · 20<sup>39</sup> | 1.0000 | Paley / Williamson (Hadamard) | 1933 | Pal1933 | — |
| 81 | 1 | OPEN | 1.0995 × 10<sup>53</sup> | 200 · 20<sup>39</sup> | 0.7881 | **we\* (@mizopl)** | 2026 | WE2026 | [txt](matrices/we_star_2026/n081.txt) |
| 82 | 2 | OPT | 8.9060 × 10<sup>53</sup> | 1620 · 20<sup>39</sup> | 1.0000 | Cohn; Koukouvinos, Kounias & Seberry | 1989 | Coh1989 | — |
| 83 | 3 | OPEN | 5.0710 × 10<sup>54</sup> | — | 0.8561 | T. Rokicki | 2005 | Rok2005 | [txt](matrices/literature_witnesses/n083.txt)* |
| 84 | 0 | OPT | 6.8272 × 10<sup>55</sup> | 42 · 21<sup>41</sup> | 1.0000 | Paley / Williamson (Hadamard) | 1933 | Pal1933 | — |
| 85 | 1 | OPEN | 3.4902 × 10<sup>56</sup> | 4509 · 21<sup>40</sup> | 0.7865 | B. Solomon & W. P. Orrick | 2003 | OS2007 | [txt](matrices/literature_witnesses/n085.txt) |
| 86 | 2 | OPT | 2.9015 × 10<sup>57</sup> | 1785 · 21<sup>41</sup> | 1.0000 | Chadjipantelis & Kounias (EW bound) | 1985 | CK1985 | — |
| 87 | 3 | OPEN | 1.6851 × 10<sup>58</sup> | — | 0.8527 | T. Rokicki | 2005 | Rok2005 | [txt](matrices/literature_witnesses/n087.txt)* |
| 88 | 0 | OPT | 2.3315 × 10<sup>59</sup> | 44 · 22<sup>43</sup> | 1.0000 | Paley / Williamson (Hadamard) | 1933 | Pal1933 | — |
| 89 | 1 | OPEN | 1.1975 × 10<sup>60</sup> | 226 · 22<sup>43</sup> | 0.7721 | **we\* (@mizopl)** | 2026 | WE2026 | [txt](matrices/we_star_2026/n089.txt) |
| 90 | 2 | OPT | 1.0375 × 10<sup>61</sup> | 1958 · 22<sup>43</sup> | 1.0000 | Cohn; Koukouvinos, Kounias & Seberry | 1989 | Coh1989 | — |
| 91 | 3 | OPEN | 6.1448 × 10<sup>61</sup> | — | 0.8501 | T. Rokicki | 2005 | Rok2005 | [txt](matrices/literature_witnesses/n091.txt) |
| 92 | 0 | OPT | 8.7199 × 10<sup>62</sup> | 46 · 23<sup>45</sup> | 1.0000 | Baumert, Golomb & Hall | 1962 | BGH1962 | — |
| 93 | 1 | OPEN | 4.6064 × 10<sup>63</sup> | 243 · 23<sup>45</sup> | 0.7768 | **we\* (@mizopl)** | 2026 | WE2026 | [txt](matrices/we_star_2026/n093.txt) |
| 94 | 2 | OPEN | 3.9312 × 10<sup>64</sup> | — | 0.9695 | T. Rokicki | 2005 | Rok2005 | [txt](matrices/literature_witnesses/n094.txt)* |
| 95 | 3 | OPEN | 2.4387 × 10<sup>65</sup> | — | 0.8447 | T. Rokicki | 2005 | Rok2005 | [txt](matrices/literature_witnesses/n095.txt)* |
| 96 | 0 | OPT | 3.5577 × 10<sup>66</sup> | 48 · 24<sup>47</sup> | 1.0000 | Paley / Williamson (Hadamard) | 1933 | Pal1933 | — |
| 97 | 1 | OPEN | 1.9567 × 10<sup>67</sup> | 264 · 24<sup>47</sup> | 0.7918 | **we\* (@mizopl)** | 2026 | WE2026 | [txt](matrices/we_star_2026/n097.txt) |
| 98 | 2 | OPT | 1.7255 × 10<sup>68</sup> | 2328 · 24<sup>47</sup> | 1.0000 | Djokovic (EW bound) | 1997 | Djo1997 | — |
| 99 | 3 | OPEN | 1.0657 × 10<sup>69</sup> | — | 0.8496 | T. Rokicki | 2005 | Rok2005 | [txt](matrices/literature_witnesses/n099.txt) |
| 100 | 0 | OPT | 1.5777 × 10<sup>70</sup> | 50 · 25<sup>49</sup> | 1.0000 | Paley / Williamson (Hadamard) | 1933 | Pal1933 | — |
| 101 | 1 | OPEN | 9.6254 × 10<sup>70</sup> | 7626 · 25<sup>48</sup> | 0.8606 | **we\* (@mizopl)** | 2026 | WE2026 | [txt](matrices/we_star_2026/n101.txt) |
| 102 | 2 | OPT | 7.9675 × 10<sup>71</sup> | 2525 · 25<sup>49</sup> | 1.0000 | Cohn (EW bound) | 1992 | Coh1992 | — |
| 103 | 3 | OPEN | 4.9829 × 10<sup>72</sup> | — | 0.8433 | Butbaia, Gukov et al. | 2026 | BG2026 | [txt](matrices/literature_witnesses/n103.txt)* |
| 104 | 0 | OPT | 7.5796 × 10<sup>73</sup> | 52 · 26<sup>51</sup> | 1.0000 | Paley / Williamson (Hadamard) | 1933 | Pal1933 | — |
| 105 | 1 | OPEN | 4.1559 × 10<sup>74</sup> | 7413 · 26<sup>50</sup> | 0.7585 | **we\* (@mizopl)** | 2026 | WE2026 | [txt](matrices/we_star_2026/n105.txt) |
| 106 | 2 | OPEN | 3.9414 × 10<sup>75</sup> | 2704 · 26<sup>51</sup> | 0.9905 | W. P. Orrick | 2005 | Orr2005 | [txt](matrices/literature_witnesses/n106.txt)* |
| 107 | 3 | OPEN | 2.5405 × 10<sup>76</sup> | — | 0.8445 | Butbaia, Gukov et al. | 2026 | BG2026 | [txt](matrices/literature_witnesses/n107.txt)* |
| 108 | 0 | OPT | 3.9325 × 10<sup>77</sup> | 54 · 27<sup>53</sup> | 1.0000 | Paley / Williamson (Hadamard) | 1933 | Pal1933 | — |
| 109 | 1 | OPEN | 2.2503 × 10<sup>78</sup> | 309 · 27<sup>53</sup> | 0.7769 | **we\* (@mizopl)** | 2026 | WE2026 | [txt](matrices/we_star_2026/n109.txt) |
| 110 | 2 | OPT | 2.1432 × 10<sup>79</sup> | 2943 · 27<sup>53</sup> | 1.0000 | Fletcher & Seberry (EW bound) | 2001 | FS2001 | — |
| 111 | 3 | OPEN | 1.3978 × 10<sup>80</sup> | — | 0.8469 | Butbaia, Gukov et al. | 2026 | BG2026 | [txt](matrices/literature_witnesses/n111.txt)* |
| 112 | 0 | OPT | 2.1973 × 10<sup>81</sup> | 56 · 28<sup>55</sup> | 1.0000 | Paley / Williamson (Hadamard) | 1933 | Pal1933 | — |
| 113 | 1 | OPT | 1.6479 × 10<sup>82</sup> | 420 · 28<sup>55</sup> | 1.0000 | Brouwer (Barba bound) | 1983 | Bro1983 | — |
| 114 | 2 | OPT | 1.2414 × 10<sup>83</sup> | 3164 · 28<sup>55</sup> | 1.0000 | Koukouvinos, Kounias & Seberry | 1991 | KKS1991 | — |
| 115 | 3 | OPEN | 8.2608 × 10<sup>83</sup> | — | 0.8488 | Butbaia, Gukov et al. | 2026 | BG2026 | [txt](matrices/literature_witnesses/n115.txt)* |
| 116 | 0 | OPT | 1.3186 × 10<sup>85</sup> | 58 · 29<sup>57</sup> | 1.0000 | Baumert | 1966 | Bau1966 | — |
| 117 | 1 | OPEN | 7.6161 × 10<sup>85</sup> | 335 · 29<sup>57</sup> | 0.7568 | **we\* (@mizopl)** | 2026 | WE2026 | [txt](matrices/we_star_2026/n117.txt) |
| 118 | 2 | OPT | 7.7138 × 10<sup>86</sup> | 3393 · 29<sup>57</sup> | 1.0000 | Fletcher, Koukouvinos & Seberry | 2004 | FKS2004 | — |
| 119 | 3 | OPEN | 5.2308 × 10<sup>87</sup> | — | 0.8503 | Butbaia, Gukov et al. | 2026 | BG2026 | [txt](matrices/literature_witnesses/n119.txt)* |
| 120 | 0 | OPT | 8.4782 × 10<sup>88</sup> | 60 · 30<sup>59</sup> | 1.0000 | Paley / Williamson (Hadamard) | 1933 | Pal1933 | — |
</details>

---

## 4. Exact Verification & Testing Your Own Matrices (`verify.py`)

### Why Floating-Point Determinants Cannot Certify Records
Standard IEEE-754 `float64` routines (`numpy.linalg.det` and `numpy.linalg.slogdet`) carry only 53 bits ($\approx 15.9$ decimal digits) of mantissa precision, whereas $D_{\text{norm}}(n)$ reaches 88 decimal digits by $n = 119$, and distinct constructions can agree to 6+ leading digits (for example at $n = 103$, where the 2026 improvement over Rokicki's 2005 circulant matrix is $+0.00004\%$ in the 7th significant digit). Every entry in this repository is therefore verified in **exact arbitrary-precision integer arithmetic**:

1. **Domain & Dimension Check:** Verifies that $X$ has exact shape $n \times n$ and every entry belongs strictly to $\{-1, +1\}$.
2. **Fraction-Free Bareiss Elimination:** Computes exact $\det X \in \mathbb{Z}$ via Bareiss integer elimination ([`verifiers/determinant_and_bounds.py`](verifiers/determinant_and_bounds.py), zero external dependencies; optional independent cross-check via `sympy.Matrix.det(method='bareiss')` with `--sympy`).
3. **Exact Scaling & Upper-Bound Inequality:** Confirms $|\det X| \equiv 0 \pmod{2^{n-1}}$, extracts $D_{\text{norm}} = |\det X| / 2^{n-1}$, and checks the exact rational inequality $|\det X|^2 \le B_{\text{class}}(n)^2$ (using the residue-appropriate Hadamard, Barba, Ehlich–Wojtas, or $s$-block Ehlich bound with $s = 3, 5, 6, 7$).
4. **Exact $\mathbb{F}_2$ Sign-Gauge Block-Circulant Audit (`*`):** Solves the exact binary linear system $s_r[i] \oplus s_c[j] \oplus c[k] = \mathbf{1}_{B_{i,j}=-1}$ over $\mathbb{F}_2$ ([`verifiers/circulant_structure.py`](verifiers/circulant_structure.py)) to certify 1-, 2-, and 4-block circulant/back-circulant cores even when archived matrices are normalized by row/column sign flips.

### Modular `verifiers/` Directory (Template Reuse & Co-Maintainer Toolkit)
All verification and table-maintenance code lives openly in [`verifiers/`](verifiers/) (full workflow in [`CONTRIBUTING.md`](CONTRIBUTING.md)):
- [`verifiers/determinant_and_bounds.py`](verifiers/determinant_and_bounds.py) — standalone Bareiss determinant & exact rational residue-bound engine.
- [`verifiers/circulant_structure.py`](verifiers/circulant_structure.py) — standalone $\mathbb{F}_2$ block-circulant classifier (`python3 verifiers/circulant_structure.py --inspect matrices/literature_witnesses/n051.txt`).
- [`verifiers/candidate_audit.py`](verifiers/candidate_audit.py) — multi-format external candidate parser & Hadamard-equivalence Gram fingerprinter.
- [`verifiers/admin_update_catalog.py`](verifiers/admin_update_catalog.py) — co-maintainer CLI for ingesting verified new records into JSON/CSV/`matrices/` and refreshing `SHA256SUMS`.

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
- **`status = PROVEN_OPTIMAL` + `attains_residue_bound = False` (9 orders: `7, 9, 11, 15, 17, 19, 21, 22, 37`):** The closed-form residue-class formula is arithmetically unattainable over $\{-1, +1\}$ ($R < 1$), yet an exhaustive computer-assisted or algebraic proof establishes that no larger determinant exists for that $n$ (e.g., [Orrick 2005](https://arxiv.org/abs/math/0401179) for $n = 15$; [Brent et al. 2011](https://arxiv.org/abs/1112.4160) for $n = 19, 37$; Chasiotis, Kounias & Farmakis 2018 for $n = 22$).
- **`status = OPEN` + `attains_residue_bound = False` (50 orders, `23 <= n <= 119`):** No optimality proof is published and $R < 1$; `D_norm_exact` records the best verified constructive lower bound.
</details>

---

## 5. Primary Sources, Related Benchmarks & Attribution

### Archival Databases & Integer Sequences
- **W. P. Orrick & B. Solomon — `MaxDet` (`indiana.edu/~maxdet`, 2002–2012):** Preserved via the [Indiana University Archive-It mirror](http://wayback.archive-it.org/219/20190316001107/http://www.indiana.edu/~maxdet/), the [Wayback full summary table (`fullPage.shtml`)](https://web.archive.org/web/20131002153624/http://www.indiana.edu/~maxdet/fullPage.shtml), and the [theoretical upper-bounds guide (`bounds.html`)](https://web.archive.org/web/20131002153511/http://www.indiana.edu/~maxdet/bounds.html).
- **OEIS Sequences:** [OEIS A003432](https://oeis.org/A003432) (normalized integer $D_{\text{norm}}(n) = |\det X| / 2^{n-1}$, equal to the maximal $(n-1)\times(n-1)$ binary $\{0,1\}$ determinant) · [OEIS A003433](https://oeis.org/A003433) (raw maximal $\{-1,+1\}$ determinant $|\det X|$).
- **Two-Circulant RL Dataset (`n ≡ 3 mod 4`):** G. Butbaia, S. Gukov et al., [`Math-AI-Caltech/hadamard-maxdet`](https://github.com/Math-AI-Caltech/hadamard-maxdet) ([arXiv:2608.22518](https://arxiv.org/abs/2608.22518) · [Zenodo 22024187](https://zenodo.org/records/22024187), Oct 2026) covering orders `51, 103, 107, 111, 115, 119`.

### Active Open-Problem Benchmarks & Agentic Challenges
- **[EinsteinArena](https://einsteinarena.com/) (*Hadamard Maximal Determinant*):** Hosts a live benchmark specifically on **order $n = 51$**, where the current record is $D_{\text{norm}}(51) = 17{,}776{,}121{,}037{,}665{,}193{,}653{,}653{,}203{,}125$ ($R = 0.877569$, [`matrices/literature_witnesses/n051.txt`](matrices/literature_witnesses/n051.txt)). Candidate $51 \times 51$ JSON matrices can be pre-verified locally with `python3 verify.py --check candidate_51.json`.
- **[`elder-plinius/V3R1T4S`](https://github.com/elder-plinius/V3R1T4S) (`hadamard-maxdet-records`):** Open verification benchmark tracking strict determinant improvements and maximality certificates across all 50 open orders $n \not\equiv 0 \pmod 4$.
- **[`teorth/optimizationproblems`](https://github.com/teorth/optimizationproblems):** Terence Tao's curated catalog of open extremal and combinatorial optimization benchmarks.

### Expository Overviews & Constructive Hadamard Libraries
- **Mathematical Surveys:** P. Browne, R. Egan, F. Hegarty, and P. Ó Catháin, [*A survey of the Hadamard maximal determinant problem*](https://arxiv.org/abs/2104.06756) (*Electron. J. Combin.* 28(4), P4.41, 2021) · [Wikipedia — *Hadamard's maximal determinant problem*](https://en.wikipedia.org/wiki/Hadamard%27s_maximal_determinant_problem) · [Wolfram MathWorld — *Hadamard's Maximum Determinant Problem*](https://mathworld.wolfram.com/HadamardsMaximumDeterminantProblem.html) · [R. P. Brent's Maximal Determinant Papers](https://maths-people.anu.edu.au/~brent/pub/pubsall.html).
- **Constructive Matrix Libraries (Order $4k$ Seeds):** [N. J. A. Sloane's Library of Hadamard Matrices](http://neilsloane.com/hadamard/) · [SageMath Combinatorial Designs (`hadamard_matrix`)](https://doc.sagemath.org/html/en/reference/combinat/sage/combinat/matrices/hadamard_matrix.html).
- **Agent, Maintainer & LLM Guides:** [`AGENTS.md`](AGENTS.md) (autonomous agent primer & verification rules) · [`CONTRIBUTING.md`](CONTRIBUTING.md) (PR, template-reuse & co-maintainer guide) · [`llms.txt`](llms.txt) ([`llmstxt.org`](https://llmstxt.org/) index) · [`CITATION.cff`](CITATION.cff).

Corrections, missing historical references, verified new matrices (`python3 verify.py --check your_matrix.txt`), or inquiries regarding editorial co-stewardship (see [`CONTRIBUTING.md`](CONTRIBUTING.md)) are welcome via GitHub Issues or Pull Requests.

<details>
<summary><strong>Full Bibliographic Key ([Syl1867] – [WE2026], BibTeX: <code>catalog/references.bib</code>)</strong></summary>

Machine-readable BibTeX entries for all references below are provided in [`catalog/references.bib`](catalog/references.bib):

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
- **[Bro1983]** A. E. Brouwer (1983), An infinite series of symmetric designs, Math. Centrum Amsterdam Report ZW 202/83; independently A. L. Whiteman (1990), Ars Combin. 30, 23-26.
- **[CK1985]** T. Chadjipantelis and S. Kounias (1985), Supplementary difference sets and D-optimal designs for n = 2 mod 4, Discrete Math. 57, 211-216.
- **[CKM1987]** T. Chadjipantelis, S. Kounias, and C. Moyssiadis (1987), The maximum determinant of 21x21 (+1,-1)-matrices and D-optimal designs, J. Statist. Plann. Inference 16, 167-178.
- **[FK1987]** N. Farmakis and S. Kounias (1987), The excess of Hadamard matrices and optimal designs, Discrete Math. 67, 165-176.
- **[Coh1989]** J. H. E. Cohn (1989), On determinants with elements +-1, II, Bull. London Math. Soc. 21, 36-42; C. Koukouvinos, S. Kounias, and J. Seberry (1989), Utilitas Math. 36, 135-150.
- **[KKS1991]** C. Koukouvinos, S. Kounias, and J. Seberry (1991), Supplementary difference sets and optimal designs, Discrete Math. 88, 49-58.
- **[Coh1992]** J. H. E. Cohn (1992), A D-optimal design of order 102, Discrete Math. 102, 61-65.
- **[Djo1997]** D. Z. Djokovic (1997), Some new D-optimal designs, Australas. J. Combin. 15, 221-231.
- **[FS2001]** R. J. Fletcher and J. Seberry (2001), New D-optimal designs of order 110, Australas. J. Combin. 23, 49-52.
- **[FKS2004]** R. J. Fletcher, C. Koukouvinos, and J. Seberry (2004), New skew-Hadamard matrices of order 4*59 and new D-optimal designs of order 2*59, Discrete Math. 286, 251-253.
- **[Sol2002]** B. Solomon (2002), Large-determinant sign matrices of orders 29 and 33, reported in W. P. Orrick's [MaxDet catalog](http://wayback.archive-it.org/219/20190316001107/http://www.indiana.edu/~maxdet/).
- **[OSDS2003]** W. P. Orrick, B. Solomon, R. Dowdeswell, and W. D. Smith (2003), New lower bounds for the maximal determinant problem, [arXiv:math/0304410](https://arxiv.org/abs/math/0304410).
- **[OS2007]** W. P. Orrick and B. Solomon (2003; publ. 2007), Large-determinant sign matrices of order 4k+1, Discrete Math. 307, 226-236 ([arXiv:math/0311292](https://arxiv.org/abs/math/0311292)).
- **[DNST2004]** A. Dowdeswell, M. G. Neubauer, B. Solomon, and K. Tumer (2004), Lower bound for order 23 in W. P. Orrick's [MaxDet catalog](https://web.archive.org/web/20131002153624/http://www.indiana.edu/~maxdet/fullPage.shtml).
- **[Orr2005]** W. P. Orrick (2005), Switching operations for Hadamard matrices and D-optimal designs of orders 34, 70, 106, SIAM J. Discrete Math. 22 (2008), 31-50 ([arXiv:math/0507515](https://arxiv.org/abs/math/0507515)).
- **[Orr2005b]** W. P. Orrick (2005), The maximal {-1, 1}-determinant of order 15, Metrika 62, 195-219 ([arXiv:math/0401179](https://arxiv.org/abs/math/0401179)).
- **[Rok2005]** T. Rokicki (Apr 2005), New determinant records from circulant pairs (orders 51-119), archived in W. P. Orrick's [MaxDet catalog](https://web.archive.org/web/20131002153624/http://www.indiana.edu/~maxdet/fullPage.shtml).
- **[Mey2005]** J.-C. Meyrignac and B. Felgenhauer (Apr 2005), Order 57 sign matrix record, archived in W. P. Orrick's [MaxDet catalog](https://web.archive.org/web/20131002153624/http://www.indiana.edu/~maxdet/fullPage.shtml).
- **[Tam2005]** H. Tamura (Aug 2005), Maximal determinant records for orders 27, 31, and 59, archived in W. P. Orrick's [MaxDet catalog](https://web.archive.org/web/20131002153624/http://www.indiana.edu/~maxdet/fullPage.shtml).
- **[BOOZ2011]** R. P. Brent, W. P. Orrick, J.-a. H. Osborn, and P. Zimmermann (2009-2011), Maximal determinants and saturated D-optimal designs of orders 19 and 37 (and order 45 construction), [arXiv:1112.4160](https://arxiv.org/abs/1112.4160).
- **[Tam2011]** H. Tamura (Jan 2011), Maximal determinant records for orders 35, 39, 43, and 63, archived in W. P. Orrick's [MaxDet catalog](https://web.archive.org/web/20131002153624/http://www.indiana.edu/~maxdet/fullPage.shtml).
- **[DK2015]** D. Z. Djokovic and I. S. Kotsireas (2012-2015), New results on D-optimal matrices, J. Combin. Des. 20, 278-289 ([arXiv:1103.3626](https://arxiv.org/abs/1103.3626)); Algebraic Design Theory and Hadamard Matrices, Springer, 61-72 ([arXiv:1501.07640](https://arxiv.org/abs/1501.07640)).
- **[CKF2018]** V. Chasiotis, S. Kounias, and N. Farmakis (2018), The D-optimal saturated designs of order 22, Discrete Math. 341, 380-387; Corrigendum ibid. 342 (2019), 2161 (recorded in [OEIS A003432](https://oeis.org/A003432) by R. P. Brent, 2021).
- **[BEHO2021]** P. Browne, R. Egan, F. Hegarty, and P. O Cathain (2021), A survey of the Hadamard maximal determinant problem, Electron. J. Combin. 28(4), P4.41 ([arXiv:2104.06756](https://arxiv.org/abs/2104.06756)).
- **[BG2026]** G. Butbaia, J. Tan, P. Vipulanandan, X. Huang, T. Saunders-A'Court, L. Fagan, D. Passaro, M. Tarquini, and S. Gukov (Oct 2026), Reinforcement learning for the Hadamard maximal determinant problem, [arXiv:2608.22518](https://arxiv.org/abs/2608.22518) · [Zenodo 22024187](https://zenodo.org/records/22024187) · [Math-AI-Caltech/hadamard-maxdet](https://github.com/Math-AI-Caltech/hadamard-maxdet).
- **[WE2026]** we* (@mizopl, Nov 2026), New constructive Hadamard maximal determinant records via asymmetric Goethals-Seidel cores, Closed-Four Hadamard switching, and Orrick-Solomon/Schur-4 seam ascent.

#### Related Recent Literature Beyond Fixed `1 <= n <= 120` Real Sign Matrices (2022–2026)
- **[MSX2022]** K. Momihara, S. Suda, and Q. Xiang (2022), *Square $(1,-1)$-matrices with large determinants and near-Hadamard matrices*, Linear Algebra Appl. 654, 28–55 (asymptotic near-Hadamard families for $n \equiv 2 \pmod 4$ via two-block almost supplementary difference sets).
- **[LMX2026]** S.-X. Li, K. Momihara, and Q. Xiang (Aug 2026), *New constructions of square $(\pm 1)$-matrices with large determinants*, Linear Algebra Appl. 752, 66–106.
- **[GS2026]** S. D. Georgiou and S. Stylianou (Oct 2026), *An exact compression-and-lifting procedure and $D$-optimal designs of orders 198 and 170*, [arXiv:2610.33418](https://arxiv.org/abs/2610.33418) (resolves $n = 170, 198$ above this catalog's $n \le 120$ window).
- **[TE2025]** Y. Wang et al. (Nov 2025), *ThetaEvolve: Test-time learning on open problems*, [arXiv:2511.23473](https://arxiv.org/abs/2511.23473) (evaluates test-time RL on open order $n = 29$).
- **Scope note on complex generalizations:** Recent papers on roots-of-unity determinants (e.g., G. Nuñez Ponasso, [arXiv:2503.11114](https://arxiv.org/abs/2503.11114), 2025; X. Dong, [arXiv:2610.40679](https://arxiv.org/abs/2610.40679), 2026) study entries in $\mu_m \subset \mathbb{C}$ ($m \ge 3$) rather than real $\{-1, +1\}$ sign matrices.
</details>

---

<sub>Engineered and verified by **we\*** (`PROOFER v20` powered by Gemini + `@mizopl`) · Released under the [MIT License](LICENSE).</sub>
