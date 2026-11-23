# Negative answers to Ta's Problems 11.6–11.7 and a proof of Conjecture 5.13

Posing locators: Lực Ta, *Good involutions of conjugation subquandles*, arXiv:2505.08090v5, Section 11, Problems 11.6–11.7, page 27; and *Good involutions of twisted conjugation subquandles and Alexander quandles*, arXiv:2508.16772v1, Section 5.2.2, Conjecture 5.13, page 17. Problems 11.6–11.7 first occur in version 3, July 15, 2025.

Haobo Ma (ChronoAI Pte Ltd, Singapore), then Wenlin Zhang (National University of Singapore; ChronoAI Pte Ltd, Singapore). CC BY 4.0.

The two semidirect products `C16 ⋊ C4` with actions 5 and 7 have invariant triple `(64,4,22)` but respectively `8847360` and `2160000` good involutions. The proper quaternion subquandle `{±i,±j}` has four good involutions and eight automorphisms, attaining `min(8,2^2)=4`. Thus both posed questions have negative answers. Maps are counted individually, not up to symmetric-quandle isomorphism.

The note extends Ta's central-label characterization to a weighted counting formula for finite conjugation subquandles. Weighted involutions of an orbit set and component factorization already appear for a linear Alexander family in Ta (arXiv:2508.16772v1, August 22, 2025), Lemma 5.10 and Theorem 5.11, pages 15–16; the formula here permits general central-twist weights. It proves the exact equality criterion for the central-label bound and the sharper whole-group bound `I(q) q^(k-q)`, where `I(q)` counts involutions on the center. It also gives the three symmetric isomorphism classes of the quaternion example and corrects the distinction between inversion-closure and classwise inverse stability in Remark 2.14.

For the cross-referenced linear family `Lambda_(4n,2n-1)`, a complete written proof gives ten good involutions for every odd positive `n` and four for every even positive `n`. The odd case proves Ta's Conjecture 5.13 (arXiv:2508.16772v1, Section 5.2.2, page 17), with `n=1` an additional case not posed there. The even case reproves his Remark 5.14. For odd `n`, all good maps are exactly `(u,v) -> (sigma(u),v)` on `Z4 × Zn`, where `sigma` is an involution of four points. For even `n`, the two parity classes independently carry labels `0` or `2n`.

The quaternion quandle is `R4`: its four good involutions reprove the corresponding case of Kamada–Oshiro (arXiv:0902.4277v1, February 25, 2009; published 2010, Theorem 3.2). The central-component decomposition is credited to Nakajima–Yamaguchi (arXiv:2608.19951v1, August 20, 2026, Theorems 1.1 and 3.3). Their separate counterexamples to Conjectures 11.8–11.9 are cited, not presented as new results.

## Files

- `note.tex` and `note.pdf`: statements, complete written proofs, class tables, source coverage and references.
- `verify.py` and `verification.json`: complete order-64 class-block counts.
- `verify_siblings.py` and `sibling-verification.json`: the quaternion examples, direct maps and isomorphism classes, all 55 nonempty quaternion conjugation subquandles, the two distinct inverse classes in `G5`, and exhaustive allowed-label counts for the linear family with `n=1,...,31`.
- `metadata.json`: authorship, license, problem identities and mathematical coverage.
- `SHA256SUMS`: member hashes in the source archive.

Both programs use only the Python 3 standard library.

## Build

With Tectonic 0.17.0:

```sh
tectonic note.tex
```

Alternatively run `pdflatex note.tex` twice. The bibliography is embedded; BibTeX is unnecessary.

## Reproduce the finite calculations

```sh
python3 verify.py --output computed.json
python3 verify_siblings.py --output computed-siblings.json
python3 -c 'import json; assert json.load(open("computed.json")) == json.load(open("verification.json")); assert json.load(open("computed-siblings.json")) == json.load(open("sibling-verification.json")); print("Exact counts, class blocks and quaternion calculations agree.")'
```

The first checker constructs both groups directly, checks the group laws, computes centers and conjugacy classes, exhausts every compatible label assignment within each central-twist block, and separately counts weighted matchings. Every accepted local assignment satisfies the defining good-involution identities. The second examines all 256 functions of the four-point quaternion subquandle and all 24 permutations for its automorphisms and antiautomorphisms. For each of the 55 nonempty conjugation-closed quaternion subsets, it determines the generated subgroup and its center, compares label enumeration with weighted matching, and checks the equality criterion. These finite examples accompany, rather than replace, the general written proofs.

## Pinned Lean declarations

All links below use revision `d9880e9f2c329bb3db982f33791419d5d90126b9`.

- [TaThreeInvariantRefutation.lean](https://github.com/the-omega-institute/trureturing/blob/d9880e9f2c329bb3db982f33791419d5d90126b9/D5/S3/FiniteGroups/Quandles/TaThreeInvariantRefutation.lean): `good_counts` states both complete numerical counts; `result` negates determination by the three invariants for finite nonabelian groups. `Gr false` is the action by 5 and `Gr true` the action by 7.
- [TaGoodInvolutionBoundRefutation.lean](https://github.com/the-omega-institute/trureturing/blob/d9880e9f2c329bb3db982f33791419d5d90126b9/D5/S3/FiniteGroups/Quandles/TaGoodInvolutionBoundRefutation.lean): `result` negates universal strictness with the four-point quaternion example, center order two, two orbits, four good maps and eight automorphisms.
- [ConjugationSubquandle.lean](https://github.com/the-omega-institute/trureturing/blob/d9880e9f2c329bb3db982f33791419d5d90126b9/D5/S3/FiniteGroups/Quandles/ConjugationSubquandle.lean): definitions `Good`, `Aut` and `orbitCount`.
- [CentralLabels.lean](https://github.com/the-omega-institute/trureturing/blob/d9880e9f2c329bb3db982f33791419d5d90126b9/D5/S3/FiniteGroups/Quandles/CentralLabels.lean): `goodLabelEquiv` and `wholeGoodEquiv` identify all whole-group good maps with central labels and the shared definition.
- [BlockCounting.lean](https://github.com/the-omega-institute/trureturing/blob/d9880e9f2c329bb3db982f33791419d5d90126b9/D5/S3/FiniteGroups/Quandles/BlockCounting.lean): `labelBlockEquiv`, `whole_good_card` and `classEquiv` give the complete whole-group block factorization and its conjugacy-class quotient.

The general subquandle formula and equality criterion, sharper center bound, three symmetric isomorphism classes, correction to Remark 2.14, and all-parameter linear-family classification are established by the written arguments. The programs provide the stated finite checks; for the linear family they exhaust the allowed labels and check the defining identities for `n=1,...,31`. The pinned Lean declarations above cover the group counts and strict-bound witness, not the linear-family classification.

## License

Copyright 2026 Haobo Ma and Wenlin Zhang. The note and accompanying materials are licensed under Creative Commons Attribution 4.0 International: https://creativecommons.org/licenses/by/4.0/.
