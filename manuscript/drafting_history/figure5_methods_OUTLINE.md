# Figure 5 + S10 Methods — outline for review

## Placement: split, not one or the other

The DE half genuinely belongs in the existing DE section; the rest does not.

**(A) Fold into `Cell type-specific differential gene expression`** — ~130 words.
That section already describes two aggregation levels (subclass, supertype) with
one pipeline. Depletion-group DE is the *same* pipeline at a third level, so
describing it anywhere else means repeating TMM → voom → limma → metafor.

**(B) New section, after Patch-seq: `Gene-set enrichment and depletion-group
interaction`** — ~550 words. GSEA, the interaction model and the cell-count
control are not differential expression, and folding ~600 words into a 399-word
section would bury the Figure 2 method. Note **GSEA has no Methods home anywhere
in the paper at present**, so it needs this section regardless of Figure 5.

Dependency to handle: the groups are defined from the abundance analysis (Fig 3a),
which is described *after* the DE section. Simplest fix — define the groups in (B)
and have (A) point forward.

---

## (A) Insert into the DE section

One paragraph, after the existing supertype sentence:

- Third aggregation level: nuclei pooled across the supertypes of each depletion
  group, per donor per dataset (groups defined in section B).
- Same retention rules: donors ≥ 10 nuclei of the group; genes ≥ 1 count in ≥ 80%
  of retained donors.
- Same model (1), same metafor REML meta-analysis, k ≥ 5, BH within group.
- Gene identity: symbols; the two Ensembl-native cohorts (HBCC, MSSM2) mapped via
  GENCODE v44, duplicate symbols summed at the count level, unmapped dropped.
- Coverage: 395 / 411 / 342 donors and 44,744 / 42,814 / 13,157 nuclei for the
  depleted / intermediate / non-depleted groups.
- The "All Sst" reference column = same procedure pooled over all 16 supertypes.

---

## (B) New section

**1. Depletion groups** (~90 w)
- From the crumblr 7-cohort compositional meta-analysis (Fig. 3a): depleted =
  β < 0 and FDR < 0.20 (Sst_2, 3, 20, 22, 25); intermediate = β < 0, not
  significant (n = 6); non-depleted = β ≥ 0 (n = 5).
- State *why* FDR < 0.20: admits Sst_20 (FDR = 0.19), which follows the same
  laminar and genetic pattern as the other four. DECISION: do we also state that
  the grouping was fixed before the expression analysis?
- Median cortical depth per supertype from Xenium; Spearman ρ = 0.74 with β.

**2. Gene-set enrichment** (~130 w)
- fgsea, preranked on z = estimate/SE from the meta-analysis.
- MSigDB via msigdbr: GO BP/CC/MF + Reactome. **Hallmark deliberately excluded** —
  say so and why (coarse meta-signatures spanning pathways; hard to interpret
  beside named pathways).
- Sets restricted to 10–500 genes; 10,000 permutations; BH within group;
  significance at FDR < 0.10.
- Blocks in Fig. 5d/e are representative sets; the module colours in Fig. 5c are
  the leading-edge union of each block in the depleted group.
- DECISION: state the gene-permutation caveat (nulls assume genes exchangeable;
  anti-conservative for co-expressed sets such as the ribosome)? Honest, and a
  reviewer may raise it — but it slightly undercuts the headline set.

**3. Diagnosis × depletion-group interaction** (~150 w)
- Key design point: every donor contributes one pseudobulk **per group**, so the
  contrast is within donor and donor-level factors cancel by construction.
- Repeated measures via `limma::duplicateCorrelation(block = donor)`; voom refit
  with the consensus correlation.
- Model (2): ~ Diagnosis × group + age + sex + PMI, reference = non-depleted →
  coefficient `Diagnosis:group_depleted`.
- Model (3): single-df linear trend, ~ Diagnosis × score, score = 0/1/2.
- metafor REML across datasets, k ≥ 5, BH within coefficient; fgsea on the
  interaction z.
- DECISION: report model (3) at all? It is null for the key sets (ribosomal
  subunit FDR = 0.32), so including it is honest disclosure that the effect is a
  depleted-vs-non contrast rather than a monotonic gradient — but the results text
  makes no trend claim, so it could be omitted.

**4. Cell-count control (Fig. S10)** (~150 w)
- Motivation: groups are *defined* by reduced abundance in SCZ, so cases
  contribute 29.5% / 17.9% / 4.0% fewer nuclei (P = 1.6 × 10⁻⁵ / 0.005 / 0.52).
- Cell-level Sst counts for all 7 cohorts (101,566 nuclei); summing them by
  donor × supertype reproduces the pseudobulks exactly.
- Within each cohort and group: restrict to donors ≥ 10 nuclei, assign each
  control donor a target count from the case distribution at that donor's
  percentile rank, subsample nuclei without replacement, rebuild pseudobulks from
  cell-level counts, rerun DE and GSEA unchanged.
- 5 independent draws. Result: imbalance → −0.4% / +1.8% / +4.9% (P ≥ 0.43).
- State that matching is on cells, not reads.

**5. Software** (~30 w)
R 4.5.1; limma 3.64.3, edgeR 4.6.3, metafor 4.8.0, fgsea 1.34.2, msigdbr 25.1.1.
Python 3.14 with anndata 0.12.10, numpy 2.4.4, scipy 1.17.1, pandas 2.3.3.
DECISION: is there a code-availability statement to point at? The pipeline is
`transcriptomic/scripts/fig5/`.

---

## Open decisions

1. Split as proposed, or force everything into the DE section?
2. State the GSEA gene-permutation caveat?
3. Report the linear-trend model (3)?
4. Note that grouping was fixed before the expression analysis?
5. Code-availability pointer for the Figure 5 pipeline?
