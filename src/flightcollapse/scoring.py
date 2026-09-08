"""Calibrated novelty model: how likely is this novel event to be real?

The problem with "require every novel event to appear at least N times" is that
N is not comparable across loci.  At GAPDH depth (240 k reads on one locus) ten
reads is noise; at a 50-read locus ten reads is the dominant isoform.  A flat
threshold is simultaneously too permissive and too strict, which is a large part
of why general-purpose tools overcall transcriptional diversity.

What is done here instead
-------------------------
1. **Features, not counts.**  Every novel junction / chain gets a feature vector
   combining depth-relative support, molecule and cell diversity (not read
   counts -- reads are PCR copies), sequence evidence (splice motif, RT-switch
   direct repeat, polyA signal) and geometry (distance to annotated sites).

2. **Supervised anchoring, no training data required.**  The annotation itself
   supplies positives: events that match GENCODE.  A decoy set supplies
   negatives: junctions whose motif is the reverse complement of a canonical one
   (a strand-assignment / alignment error under a stranded protocol), plainly
   non-canonical junctions, and events carrying an RT-switch direct repeat.
   A ridge logistic regression on those two anchor sets gives a score.

3. **Two-component mixture -> local FDR.**  The scores of *novel* events are
   modelled as a mixture ``f(s) = pi0 * f0(s) + (1 - pi0) * f1(s)`` where
   ``f0`` is estimated from the decoys.  The local false discovery rate

       lfdr(s) = pi0 * f0(s) / f(s)

   is the posterior probability that an event with score ``s`` is an artefact,
   i.e. ``1 - lfdr`` is the posterior probability it is real.  Thresholding on
   ``lfdr`` gives a *tunable, calibrated* discovery rate rather than an
   arbitrary count.

4. **Stratification.**  lfdr is estimated separately for events with 2 / 1 / 0
   annotated splice sites, because a novel combination of two known sites and a
   junction with two invented sites are different populations, and pooling them
   would let the abundant former mask the latter.

Everything degrades gracefully: with too few anchors, or with ``scoring.enabled
= false``, the model returns ``NaN`` and the pipeline falls back to explicit
thresholds -- which are never silently disabled (see :meth:`Config.validate`).
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np
import pandas as pd

from .config import ScoringParams

EPS = 1e-9


# ---------------------------------------------------------------------- #
# ridge logistic regression (IRLS) -- no sklearn dependency
# ---------------------------------------------------------------------- #
def fit_logistic(
    X: np.ndarray, y: np.ndarray, l2: float = 1.0, max_iter: int = 50, tol: float = 1e-7
) -> np.ndarray:
    """Newton-IRLS with an L2 penalty on all but the intercept.

    ``X`` must already include an intercept column of ones as column 0.
    Returns the coefficient vector.  Falls back to a damped step if the Hessian
    is ill-conditioned (perfectly separating features).
    """
    n, p = X.shape
    w = np.zeros(p)
    pen = np.full(p, l2)
    pen[0] = 0.0
    for _ in range(max_iter):
        eta = np.clip(X @ w, -30, 30)
        mu = 1.0 / (1.0 + np.exp(-eta))
        s = np.clip(mu * (1 - mu), 1e-6, None)
        g = X.T @ (y - mu) - pen * w
        H = (X * s[:, None]).T @ X + np.diag(pen)
        try:
            step = np.linalg.solve(H + 1e-8 * np.eye(p), g)
        except np.linalg.LinAlgError:  # pragma: no cover
            step = np.linalg.lstsq(H + 1e-6 * np.eye(p), g, rcond=None)[0]
        w_new = w + step
        if np.max(np.abs(w_new - w)) < tol:
            w = w_new
            break
        w = w_new
    return w


def predict_logit(X: np.ndarray, w: np.ndarray) -> np.ndarray:
    return np.clip(X @ w, -30, 30)


# ---------------------------------------------------------------------- #
# local FDR from a decoy-anchored two-component mixture
# ---------------------------------------------------------------------- #
def _smooth(v: np.ndarray, k: int) -> np.ndarray:
    if k <= 1:
        return v
    kern = np.ones(k) / k
    num = np.convolve(v, kern, mode="same")
    norm = np.convolve(np.ones_like(v), kern, mode="same")
    return num / np.maximum(norm, EPS)


def _pava_nonincreasing(y: np.ndarray, w: np.ndarray) -> np.ndarray:
    """Weighted isotonic regression to a non-increasing sequence (PAVA).

    Used instead of a running minimum: a running minimum locks in downward
    noise from low-count bins, which systematically *under*-states the local
    FDR in exactly the region where the estimate is least reliable.
    """
    yr = y[::-1].astype(float)
    wr = np.maximum(w[::-1].astype(float), EPS)
    vals: List[float] = []
    wts: List[float] = []
    sizes: List[int] = []
    for v, ww in zip(yr, wr):
        vals.append(float(v))
        wts.append(float(ww))
        sizes.append(1)
        while len(vals) > 1 and vals[-2] > vals[-1]:
            v2, w2, s2 = vals.pop(), wts.pop(), sizes.pop()
            v1, w1, s1 = vals.pop(), wts.pop(), sizes.pop()
            vals.append((v1 * w1 + v2 * w2) / (w1 + w2))
            wts.append(w1 + w2)
            sizes.append(s1 + s2)
    out = np.empty(yr.size, float)
    i = 0
    for v, s in zip(vals, sizes):
        out[i:i + s] = v
        i += s
    return out[::-1]


def local_fdr(
    scores: np.ndarray,
    decoy_scores: np.ndarray,
    n_bins: int = 60,
    smooth_window: int = 5,
    pi0_method: str = "decoy",
    pi0_fixed: float = 0.8,
) -> Tuple[np.ndarray, float]:
    """``lfdr`` per element of ``scores`` plus the estimated ``pi0``.

    ``scores`` are the scores of the *candidate* population (novel events);
    ``decoy_scores`` are drawn from events believed to be artefacts.
    """
    if scores.size == 0:
        return np.zeros(0), float("nan")
    if decoy_scores.size < 10:
        return np.full(scores.size, np.nan), float("nan")

    lo = float(min(scores.min(), decoy_scores.min()))
    hi = float(max(scores.max(), decoy_scores.max()))
    if hi - lo < 1e-9:
        return np.full(scores.size, np.nan), float("nan")
    edges = np.linspace(lo, hi, n_bins + 1)
    f, _ = np.histogram(scores, bins=edges, density=True)
    f0, _ = np.histogram(decoy_scores, bins=edges, density=True)
    f = _smooth(f, smooth_window) + EPS
    f0 = _smooth(f0, smooth_window) + EPS

    if pi0_method == "fixed":
        pi0 = float(np.clip(pi0_fixed, 0.0, 1.0))
    else:
        # Storey-style tail counting, which is far more stable than a ratio of
        # density estimates.  In the low-score region essentially everything is
        # an artefact, so the candidate and decoy cumulative fractions there
        # differ exactly by pi0.  Averaged over several cut-points for
        # robustness.
        est = []
        for q in (0.1, 0.2, 0.3, 0.4):
            cut = float(np.quantile(decoy_scores, q))
            f_cand = float(np.mean(scores <= cut))
            if q > 0:
                est.append(f_cand / q)
        pi0 = float(np.clip(np.median(est), 0.0, 1.0)) if est else 1.0

    lf_bin = np.clip(pi0 * f0 / f, 0.0, 1.0)
    # a higher score can never be more suspicious, so the lfdr must be
    # non-increasing across bins of increasing score
    counts = np.histogram(scores, bins=edges)[0].astype(float)
    lf_bin = np.clip(_pava_nonincreasing(lf_bin, counts), 0.0, 1.0)
    idx = np.clip(np.digitize(scores, edges) - 1, 0, n_bins - 1)
    return lf_bin[idx], pi0


# ---------------------------------------------------------------------- #
# feature construction
# ---------------------------------------------------------------------- #
def _logit(p: np.ndarray) -> np.ndarray:
    p = np.clip(p.astype(float), 1e-4, 1 - 1e-4)
    return np.log(p / (1 - p))


JUNCTION_FEATURES: Tuple[str, ...] = (
    "log_support",
    "log_cells",
    "log_reads_per_mol",
    "logit_donor_share",
    "logit_acceptor_share",
    "log_intron_len",
    "direct_repeat",
    "is_canonical",
    "is_semi_canonical",
    "log_dist_ann_donor",
    "log_dist_ann_acceptor",
)


def junction_design(feat: pd.DataFrame) -> np.ndarray:
    n = len(feat)
    reads = feat["n_reads"].to_numpy(float)
    mols = feat["n_mols"].to_numpy(float)
    support = np.where(mols > 0, mols, reads)
    cols = {
        "log_support": np.log10(1.0 + support),
        "log_cells": np.log10(1.0 + feat["n_cells"].to_numpy(float)),
        "log_reads_per_mol": np.log10(np.where(mols > 0, reads / np.maximum(mols, 1), 1.0)),
        "logit_donor_share": _logit(feat["donor_share"].to_numpy()),
        "logit_acceptor_share": _logit(feat["acceptor_share"].to_numpy()),
        "log_intron_len": np.log10(1.0 + feat["intron_len"].to_numpy(float)),
        "direct_repeat": np.minimum(feat["direct_repeat"].to_numpy(float), 20.0),
        "is_canonical": (feat["motif_class"] == "canonical").to_numpy(float),
        "is_semi_canonical": (feat["motif_class"] == "semi_canonical").to_numpy(float),
        "log_dist_ann_donor": np.log10(1.0 + feat["dist_ann_donor"].to_numpy(float)),
        "log_dist_ann_acceptor": np.log10(1.0 + feat["dist_ann_acceptor"].to_numpy(float)),
    }
    X = np.column_stack([np.ones(n)] + [cols[k] for k in JUNCTION_FEATURES])
    return np.nan_to_num(X, nan=0.0, posinf=0.0, neginf=0.0)


def junction_anchors(feat: pd.DataFrame) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """``(positive, decoy, exempt)`` masks for junctions.

    ``positive`` trains the score.  ``exempt`` is the set that is not a
    discovery at all and so receives ``lfdr = 0`` without being scored.  They
    used to be the same mask; short reads make them different, and conflating
    them would quietly turn the anchor into a filter.

    Annotated junctions are both positive and exempt.  A short-read-supported
    novel junction is positive -- it is a true positive that annotation cannot
    supply, and it is what stops the model from being a pure annotation
    detector -- but it is *not* exempt unless rescue is explicitly enabled, so
    it still has to earn its local FDR from long-read evidence alone.

    Short-read support also removes a junction from the decoy set. A
    non-canonical motif with fifty uniquely-mapping short reads behind it is a
    real minor-class junction, and leaving it in the null would bias the whole
    calibration.
    """
    ann = feat["annotated"].to_numpy()
    cls = feat["motif_class"].to_numpy()
    anchor = _col(feat, "sr_anchor", len(feat))
    supported = _col(feat, "sr_supported", len(feat))
    rescue = _col(feat, "sr_rescue_enabled", len(feat))
    pos = ann | anchor
    # `semi_canonical` belongs here and was missing, which is why the fitted
    # score could not tell the class apart: annotated GC-AG junctions are real
    # and sit in the positive set, so the model learned "semi-canonical is
    # fine" and applied it to the novel ones -- which the short reads confirm at
    # 0.13% against 33.4% for the annotated ones. The two populations share a
    # motif and nothing else. 97.1% of novel semi-canonical junctions were
    # scored at lfdr <= 0.05 before this.
    decoy = (
        ~ann
        & ~supported
        & (
            (cls == "antisense")
            | (cls == "non_canonical")
            | (cls == "semi_canonical")
            | (feat["direct_repeat"].to_numpy() >= 10)
        )
    )
    exempt = ann | (supported & rescue)
    return pos, decoy, exempt


def _col(feat: pd.DataFrame, name: str, n: int) -> np.ndarray:
    return (
        feat[name].to_numpy().astype(bool) if name in feat.columns
        else np.zeros(n, bool)
    )


#: Chain-level features.
#:
#: ``min_junction_posterior`` and ``frac_junctions_annotated`` are deliberately
#: absent: half the decoy set is "carries a junction that failed curation",
#: which those two separate perfectly by construction, so including them would
#: make the model rediscover its own null.  ``log_support`` is kept, because
#: excluding the single most informative feature makes the score blind to the
#: difference between a 2,000-read novel isoform and a singleton -- but note
#: that support also enters the decoy definition, so the chain-level score is
#: mildly circular and is therefore **advisory by default**
#: (``chains.max_chain_lfdr = None``).  See docs/ALGORITHM.md section 3.
CHAIN_FEATURES: Tuple[str, ...] = (
    "log_support",
    "log_cells",
    "log_reads_per_mol",
    "logit_gene_share",
    "log_n_introns",
    "log_dist_ann_tts",
    "log_dist_ann_tss",
    "polya_evidence",
)


def chain_design(feat: pd.DataFrame) -> np.ndarray:
    n = len(feat)
    reads = feat["n_reads"].to_numpy(float)
    mols = feat["n_mols"].to_numpy(float)
    support = np.where(mols > 0, mols, reads)
    cols = {
        "log_support": np.log10(1.0 + support),
        "log_cells": np.log10(1.0 + feat["n_cells"].to_numpy(float)),
        "log_reads_per_mol": np.log10(np.where(mols > 0, reads / np.maximum(mols, 1), 1.0)),
        "logit_gene_share": _logit(feat["gene_share"].to_numpy()),
        "log_n_introns": np.log10(1.0 + feat["n_introns"].to_numpy(float)),
        "log_dist_ann_tts": np.log10(1.0 + np.abs(feat["dist_ann_tts"].to_numpy(float))),
        "log_dist_ann_tss": np.log10(1.0 + np.abs(feat["dist_ann_tss"].to_numpy(float))),
        "polya_evidence": feat["polya_evidence"].to_numpy(float),
    }
    X = np.column_stack([np.ones(n)] + [cols[k] for k in CHAIN_FEATURES])
    return np.nan_to_num(X, nan=0.0, posinf=0.0, neginf=0.0)


#: novel chains at or below this support join the decoy set
CHAIN_DECOY_MAX_SUPPORT = 2


def chain_anchors(feat: pd.DataFrame) -> Tuple[np.ndarray, np.ndarray]:
    """Positives: chains matching an annotated intron chain.

    Decoys: novel chains that either carry a junction which failed curation, or
    sit at singleton/doubleton support.  52% of distinct chains are singletons
    and together they carry 1.4% of read mass, so that population is
    overwhelmingly artefactual and is a usable stand-in for the null.
    """
    pos = feat["chain_annotated"].to_numpy()
    novel = ~pos
    reads = feat["n_reads"].to_numpy(float)
    mols = feat["n_mols"].to_numpy(float)
    support = np.where(mols > 0, mols, reads)
    decoy = novel & (
        (feat["n_decoy_junctions"].to_numpy() > 0)
        | (support <= CHAIN_DECOY_MAX_SUPPORT)
    )
    return pos, decoy


# ---------------------------------------------------------------------- #
@dataclass
class NoveltyModel:
    """Fitted score + calibrated local FDR for one level (junction or chain)."""

    name: str
    coef: Optional[np.ndarray] = None
    feature_names: Tuple[str, ...] = ()
    pi0: Dict[int, float] = field(default_factory=dict)
    n_pos: int = 0
    n_decoy: int = 0
    n_exempt: int = 0
    fitted: bool = False
    note: str = ""
    #: out-of-model validation from the held-out short-read junctions
    holdout: Dict[str, object] = field(default_factory=dict)

    # ------------------------------------------------------------------ #
    def fit_predict(
        self,
        feat: pd.DataFrame,
        design_fn,
        anchor_fn,
        feature_names: Tuple[str, ...],
        params: ScoringParams,
        stratum: Optional[np.ndarray] = None,
    ) -> pd.DataFrame:
        """Return a frame with ``score``, ``lfdr`` and ``posterior_real``."""
        self.feature_names = feature_names
        out = pd.DataFrame(index=feat.index)
        out["score"] = np.nan
        out["lfdr"] = np.nan
        out["posterior_real"] = np.nan
        if not params.enabled or len(feat) == 0:
            self.note = "scoring disabled"
            return out

        X = design_fn(feat)
        anchors = anchor_fn(feat)
        pos, decoy = anchors[0], anchors[1]
        exempt = anchors[2] if len(anchors) > 2 else pos
        self.n_pos = int(pos.sum())
        self.n_decoy = int(decoy.sum())
        self.n_exempt = int(exempt.sum())
        if self.n_pos < params.min_anchor_size or self.n_decoy < params.min_anchor_size:
            self.note = (
                f"too few anchors (pos={self.n_pos}, decoy={self.n_decoy}, "
                f"need {params.min_anchor_size}); falling back to thresholds"
            )
            return out

        anchor = pos | decoy
        y = pos[anchor].astype(float)
        Xa = X[anchor]
        mu = Xa.mean(axis=0)
        sd = Xa.std(axis=0)
        sd[sd < 1e-9] = 1.0
        mu[0], sd[0] = 0.0, 1.0  # keep the intercept
        self.coef = fit_logistic((Xa - mu) / sd, y, l2=params.l2, max_iter=params.max_iter)
        self._mu, self._sd = mu, sd
        score = predict_logit((X - mu) / sd, self.coef)
        out["score"] = score
        self.fitted = True

        strat = np.zeros(len(feat), int) if stratum is None else np.asarray(stratum, int)
        lf = np.full(len(feat), np.nan)
        for s in np.unique(strat):
            m = strat == s
            # candidates are everything not exempt -- which includes the
            # short-read-anchored positives, so they are scored like anything
            # else and the held-out half is directly comparable to them
            cand = m & ~exempt
            dec = m & decoy
            if cand.sum() == 0:
                continue
            if dec.sum() < 10:
                # borrow the decoy distribution from the whole contig
                dec = decoy
            vals, pi0 = local_fdr(
                score[cand],
                score[dec],
                n_bins=params.n_score_bins,
                smooth_window=params.smooth_window,
                pi0_method=params.pi0_method,
                pi0_fixed=params.pi0_fixed,
            )
            lf[np.flatnonzero(cand)] = vals
            self.pi0[int(s)] = pi0
        lf[exempt] = 0.0  # annotated events are not discoveries
        out["lfdr"] = lf
        out["posterior_real"] = 1.0 - lf
        self._holdout_report(feat, lf)
        return out

    # ------------------------------------------------------------------ #
    def _holdout_report(self, feat: pd.DataFrame, lf: np.ndarray) -> None:
        """Recall on junctions an orthogonal library confirmed but the model never saw.

        This is the one number in the whole package that is not, in some way,
        the long-read data grading itself.  A held-out short-read-supported
        novel junction is known-real by evidence the model had no access to, so
        the fraction of them that clear an lFDR threshold is an honest recall
        estimate -- and the fraction of *unsupported* novel junctions clearing
        the same threshold bounds how much of what is kept could be noise.
        """
        if "sr_holdout" not in feat.columns:
            return
        hold = feat["sr_holdout"].to_numpy().astype(bool)
        if hold.sum() < 20:
            return
        ann = feat["annotated"].to_numpy().astype(bool)
        sup = feat["sr_supported"].to_numpy().astype(bool)
        unsup = ~ann & ~sup
        rep: Dict[str, object] = {"n_holdout": int(hold.sum()),
                                  "n_novel_unsupported": int(unsup.sum())}
        for t in (0.01, 0.05, 0.1):
            h = lf[hold]
            u = lf[unsup]
            h = h[~np.isnan(h)]
            u = u[~np.isnan(u)]
            rep[f"recall_at_lfdr_{t}"] = (
                round(float(np.mean(h <= t)), 4) if h.size else None
            )
            rep[f"unsupported_novel_kept_at_lfdr_{t}"] = (
                round(float(np.mean(u <= t)), 4) if u.size else None
            )
        self.holdout = rep

    # ------------------------------------------------------------------ #
    def report(self) -> Dict[str, object]:
        d: Dict[str, object] = {
            "level": self.name,
            "fitted": self.fitted,
            "n_positive_anchors": self.n_pos,
            "n_decoy_anchors": self.n_decoy,
            "pi0_by_stratum": {str(k): round(v, 4) for k, v in self.pi0.items()},
            "note": self.note,
        }
        if self.holdout:
            d["short_read_holdout"] = self.holdout
        if self.coef is not None:
            d["coefficients"] = {
                n: round(float(c), 4)
                for n, c in zip(("intercept",) + tuple(self.feature_names), self.coef)
            }
        return d


def score_junctions(feat: pd.DataFrame, params: ScoringParams) -> Tuple[pd.DataFrame, NoveltyModel]:
    m = NoveltyModel("junction")
    res = m.fit_predict(
        feat,
        junction_design,
        junction_anchors,
        JUNCTION_FEATURES,
        params,
        stratum=feat["n_sites_annotated"].to_numpy(int),
    )
    return res, m


def score_chains(feat: pd.DataFrame, params: ScoringParams) -> Tuple[pd.DataFrame, NoveltyModel]:
    m = NoveltyModel("chain")
    strat = np.where(
        feat["n_novel_junctions"].to_numpy() == 0, 0,
        np.where(feat["n_novel_junctions"].to_numpy() == 1, 1, 2),
    )
    res = m.fit_predict(feat, chain_design, chain_anchors, CHAIN_FEATURES, params, stratum=strat)
    return res, m
