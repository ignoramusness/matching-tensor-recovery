# Dimensions and generic recovery of matching tensors

Preprint by Ignoramus, revision 4 (10 November 2026), with its verification programs.

**Version of record:** Zenodo, [doi:10.5281/zenodo.24114656](https://doi.org/10.5281/zenodo.24114656). This repository mirrors that record so the programs can be browsed and so errors can be reported in the issue tracker.

## What the paper shows

Give each vertex i of a graph a complex vector space V_i, place a tensor W_ij in V_i ⊗ V_j on each edge, and sum the tensor products of the edge tensors over all perfect matchings. The paper studies the image of this matching-tensor map and the recovery of the edge tensors from its output. All results are over the complex numbers.

- **Any support graph with a perfect matching, local dimensions at least three:** the image closure has dimension N_* - |E_*| + dim P(G) + 1, where E_* is the set of edges lying in some perfect matching, N_* is the number of parameters on those edges and P(G) is the perfect matching polytope.
- **Complete support, every even order at least four, local dimensions at least three:** the entire fibre over a general output is a single vertex-rescaling orbit, and a normalized representative is a rational function of the output.
- **Complete support, every even order at least six, local dimensions at least two:** the dimension formula and the vertex-rescaling kernel remain valid.
- **Four vertices, local dimensions at least two:** the dimension is determined in every format.

The paper lists what it does not prove, including global recovery when some factor is two-dimensional and global recovery on sparse support. It does not decide the question about GHZ states discussed in its introduction.

Earlier results are credited in Section 1 and in the references of the paper.

## Contents

| Path | What it is |
|---|---|
| `paper.pdf` | The manuscript (25 pages) |
| `paper.tex` | Its LaTeX source |
| `verification/verify_qubit_certificates.py` | Recomputes the two exact determinants used in the proof of Theorem 5.1 |
| `verification/verify_qubit_from_matchings.py` | A second implementation that rebuilds the same matrices from the perfect-matching expansion |
| `verification/qubit-certificate-results.json` | Recorded results of those two programs |
| `verification/check_graph.py`, `graph-check-results.json` | Supplementary rank checks modulo a prime |
| `verification/independent_rank_check.py`, `independent-rank-results.json` | Supplementary rank check on complete support |
| `verification/verify_four.py`, `four_certificates.json`, `four_exact_verification.json` | Supplementary four-vertex checks |

The files are the same as in the Zenodo record, which stores them without the `verification` folder. The SHA-256 values of the five programs are printed in the final section of the paper.

## Running the checks

The programs need only Python 3 and its standard library:

```text
python verification/verify_qubit_certificates.py
python verification/verify_qubit_from_matchings.py
python verification/check_graph.py
python verification/verify_four.py
python verification/independent_rank_check.py
```

The first two print the determinants 128 and -262144. The other three write their result files beside the script. Only the first two support a proof in the paper; no proof depends on the other three.

## Disclosure

AI systems from OpenAI, working under the direction of the person who publishes as Ignoramus, developed the statements and proofs, performed the computations and wrote the first three versions. Each version was reviewed by Claude (Anthropic), also an AI system, and the edits of this fourth version were written and checked by Claude agents. The final section of the paper sets this out in full.

No AI system is an author; the author publishes under a pseudonym. The author has not checked the proofs, and no human mathematician has checked the manuscript. The AI reviews are not human refereeing, peer review or a certification of originality, and novelty and significance are not established. The manuscript is a draft for mathematical review.

## Corrections

Please open an issue, or write to ignoramus110926@gmail.com.

## Licence

Creative Commons Attribution 4.0 International (CC BY 4.0), https://creativecommons.org/licenses/by/4.0/
