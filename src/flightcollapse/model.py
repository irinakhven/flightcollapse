"""The transcript model produced by the collapse, and 3'-end evidence helpers."""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np

from .genome import Genome, PolyASiteAtlas
from .intervals import Chain, Exon, exonic_length


@dataclass
class TranscriptModel:
    model_id: str
    contig: str
    strand: str
    exons: List[Exon]
    chain: Chain
    category: str
    reads: np.ndarray = field(default_factory=lambda: np.zeros(0, np.int64))
    n_reads: int = 0
    n_mols: int = 0
    n_cells: int = 0
    #: transcript whose intron chain the model matches **exactly**; None for
    #: everything novel.  This is a lookup, not a similarity search.
    parent_tx: Optional[str] = None
    #: transcript of the same gene the model is structurally *closest* to, which
    #: for a novel model is the only reference handle it has.  Annotation only:
    #: it never feeds back into merging, thresholding or the novelty
    #: calibration, because that calibration is fitted against annotation and
    #: the loop would be circular.
    associated_tx: Optional[str] = None
    structural_diff: str = ""
    n_equally_close: int = 0
    #: whether the difference from the reference lies entirely outside the CDS.
    #: None where the question has no answer (non-coding reference transcript).
    is_utr_only: Optional[bool] = None
    gene_id: Optional[str] = None
    gene_name: str = ""
    posterior_real: float = float("nan")
    lfdr: float = float("nan")
    flags: List[str] = field(default_factory=list)
    evidence: Dict[str, object] = field(default_factory=dict)
    pb_id: str = ""

    @property
    def start(self) -> int:
        return self.exons[0][0]

    @property
    def end(self) -> int:
        return self.exons[-1][1]

    @property
    def n_exons(self) -> int:
        return len(self.exons)

    @property
    def tx_len(self) -> int:
        return exonic_length(self.exons)

    @property
    def span(self) -> int:
        return self.end - self.start

    @property
    def support(self) -> int:
        return self.n_mols if self.n_mols > 0 else self.n_reads

    def sort_key(self) -> Tuple[str, int, int]:
        return (self.contig, self.start, self.end)


# ---------------------------------------------------------------------- #
def polya_evidence(
    contig: str,
    pos: int,
    strand: str,
    genome: Optional[Genome],
    atlas: Optional[PolyASiteAtlas],
    motif_window: int = 50,
    perc_a_window: int = 20,
    max_perc_a: float = 60.0,
    atlas_window: int = 50,
    motif_min_dist: int = 0,
    motif_max_dist: int = 10 ** 6,
    motif_strict: bool = False,
    require_atlas_when_available: bool = False,
) -> Dict[str, object]:
    """Is this 3' end a genuine cleavage site, or an internal-priming event?

    For an oligo-dT protocol an internal 3' end should not exist: priming is at
    the polyA tail, so a 3' end is either a real cleavage site or an internal
    priming artefact.  Three independent readouts are combined:

    * a polyA signal hexamer 10-50 nt upstream (sequence evidence for cleavage)
    * proximity to a catalogued polyA site, when an atlas BED is supplied
    * genomic A-content immediately downstream, which is *high* for internal
      priming.  Note the 0-100 percent scale: the SQANTI cutoff is 60, not 0.6.
    """
    out: Dict[str, object] = {
        "polya_motif_found": False,
        "polya_motif": "",
        "polya_motif_dist": -1,
        "dist_to_polya_site": None,
        "perc_a_downstream": float("nan"),
        "internal_priming": False,
        "polya_supported": False,
    }
    if genome is not None:
        found, motif, dist = genome.polya_motif(
            contig, pos, strand, motif_window,
            min_dist=motif_min_dist, max_dist=motif_max_dist, strict=motif_strict,
        )
        out["polya_motif_found"] = bool(found)
        out["polya_motif"] = motif
        out["polya_motif_dist"] = int(dist)
        pa = genome.perc_a_downstream(contig, pos, strand, perc_a_window)
        out["perc_a_downstream"] = float(pa)
        out["internal_priming"] = bool(pa == pa and pa >= max_perc_a)
    if atlas is not None:
        d = atlas.distance(contig, pos, strand)
        out["dist_to_polya_site"] = int(d)
    has_atlas = out["dist_to_polya_site"] is not None
    near_site = has_atlas and out["dist_to_polya_site"] <= atlas_window
    if has_atlas and require_atlas_when_available:
        # A catalogued cleavage site is direct evidence; a hexamer is a guess,
        # and 3'UTRs are AT-rich enough to supply one by chance. When the atlas
        # is there, it decides.
        supported = near_site
    else:
        supported = near_site or bool(out["polya_motif_found"])
    out["polya_supported"] = bool(supported and not out["internal_priming"])
    return out
