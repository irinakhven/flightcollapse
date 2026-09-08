import numpy as np
import os
import pytest

from flightcollapse.molecules import (
    MoleculeIndex,
    collapse_umis,
    decode_umi,
    encode_umi,
    hamming1_neighbours,
    n_cells,
    n_molecules,
)


def test_umi_encoding_roundtrip():
    for u in ("ACGT", "AAAAAAAAAA", "GATTACAGATTACA"):
        assert decode_umi(encode_umi(u)) == u


def test_umi_encoding_handles_non_acgt():
    code = encode_umi("ACGN")
    assert decode_umi(code) is None          # hashed fallback, not decodable
    assert code != encode_umi("ACGT")


def test_hamming_neighbours_are_distance_one():
    u = "ACGTAC"
    nb = set(hamming1_neighbours(encode_umi(u)))
    assert len(nb) == 3 * len(u)
    for code in nb:
        v = decode_umi(code)
        assert v is not None and sum(a != b for a, b in zip(u, v)) == 1


def test_directional_collapse_absorbs_sequencing_errors():
    pairs = [(1, encode_umi("ACGTAC"))] * 20 + [(1, encode_umi("ACGTAG"))]
    assert n_molecules(pairs, hamming=1) == 1
    assert n_molecules(pairs, hamming=0) == 2


def test_directional_collapse_keeps_comparable_neighbours():
    """Two real molecules that happen to be neighbours must not be merged."""
    pairs = [(1, encode_umi("ACGTAC"))] * 10 + [(1, encode_umi("ACGTAG"))] * 9
    assert n_molecules(pairs, hamming=1) == 2


def test_molecules_are_per_cell():
    u = encode_umi("ACGTAC")
    pairs = [(1, u), (2, u), (3, u)]
    assert n_molecules(pairs) == 3
    assert n_cells(pairs) == 3


def test_pcr_duplicates_collapse_to_one_molecule():
    """40 reads of one molecule must not look like support from 40 molecules."""
    pairs = [(7, encode_umi("TTTTTTTTTT"))] * 40
    assert len(pairs) == 40
    assert n_molecules(pairs) == 1


def test_index_build_and_lookup(tmp_path):
    tsv = tmp_path / "bc.tsv"
    names = [f"m/{i}/ccs" for i in range(500)]
    with open(tsv, "w") as fh:
        fh.write("cell_barcode\tumi\tread_name\n")
        for i, n in enumerate(names):
            fh.write(f"CB{i % 7}\tACGTACGTAC\t{n}\n")
    idx = MoleculeIndex.build(str(tsv), verbose=False)
    assert len(idx) == 500
    cid, umi = idx.lookup(names[:10] + ["not/a/read"])
    assert (cid[:10] == np.array([i % 7 for i in range(10)], np.uint32)).all()
    assert cid[10] == np.uint32(0xFFFFFFFF)
    assert decode_umi(int(umi[0])) == "ACGTACGTAC"


def test_read_name_normalisation_modes():
    from flightcollapse.molecules import normalise_name

    n = "m84151_260724_132249_s1/100008021/ccs/2645_3379"
    assert normalise_name(n, "none") == n
    assert normalise_name(n, "strip_segment") == "m84151_260724_132249_s1/100008021/ccs"
    assert normalise_name(n, "strip_ccs") == "m84151_260724_132249_s1/100008021"
    assert normalise_name(n, "zmw") == "m84151_260724_132249_s1/100008021"
    with pytest.raises(ValueError):
        normalise_name(n, "nope")


def test_index_matches_across_a_segmentation_rename(tmp_path):
    """BAM has segment spans, the table does not -- a silent mismatch here
    would zero every molecule count."""
    tsv = tmp_path / "bc.tsv"
    with open(tsv, "w") as fh:
        fh.write("cell_barcode\tumi\tbarcode_umi\tread_name\n")
        for i in range(100):
            fh.write(f"CB{i % 5}\tACGTACGTAC\tx\tmovie/{i}/ccs\n")
    idx = MoleculeIndex.build(str(tsv), normalise="strip_segment", verbose=False)
    bam_names = [f"movie/{i}/ccs/100_900" for i in range(100)]
    cid, _u = idx.lookup(bam_names)
    assert (cid != np.uint32(0xFFFFFFFF)).all()

    plain = MoleculeIndex.build(str(tsv), verbose=False)
    cid2, _ = plain.lookup(bam_names)
    assert (cid2 == np.uint32(0xFFFFFFFF)).all()


def test_index_save_load(tmp_path):
    tsv = tmp_path / "bc.tsv"
    with open(tsv, "w") as fh:
        fh.write("cell_barcode\tumi\tread_name\n")
        for i in range(50):
            fh.write(f"CB{i}\tACGTACGTAC\tr/{i}/ccs\n")
    a = MoleculeIndex.open_or_build(str(tsv), verbose=False)
    b = MoleculeIndex.open_or_build(str(tsv), verbose=False)   # cached path
    assert (a.read_hash == b.read_hash).all()
    assert a.barcodes == b.barcodes


def test_prebuilt_index_is_actually_reused(tmp_path, capsys):
    """`index-umi` then `run` must agree on the cache filename."""
    from flightcollapse.cli import main
    from flightcollapse.molecules import MoleculeIndex, default_index_path

    tsv = tmp_path / "bc.tsv"
    with open(tsv, "w") as fh:
        fh.write("cell_barcode\tumi\tread_name\n")
        for i in range(200):
            fh.write(f"CB{i % 4}\tACGTACGTAC\tm/{i}/ccs\n")

    assert main(["index-umi", str(tsv)]) == 0
    expected = default_index_path(str(tsv), "none")
    assert os.path.exists(expected), capsys.readouterr().out

    capsys.readouterr()
    MoleculeIndex.open_or_build(str(tsv), verbose=True)
    assert "loading cached index" in capsys.readouterr().out
