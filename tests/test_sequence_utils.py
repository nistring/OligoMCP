"""Tests for sequence_utils."""
import numpy as np

from oligoclaude.sequence_utils import one_hot_encode, reverse_complement


def test_reverse_complement_basic():
    assert reverse_complement("ACGT") == "ACGT"
    assert reverse_complement("AAAA") == "TTTT"
    assert reverse_complement("GCAT") == "ATGC"


def test_reverse_complement_roundtrip():
    seq = "TTCCACTTCCTTTCTGTAAGC"
    assert reverse_complement(reverse_complement(seq)) == seq


def test_reverse_complement_lowercase():
    assert reverse_complement("acgt") == "acgt"
    assert reverse_complement("aAcC") == "GgTt"


def test_reverse_complement_preserves_n():
    assert reverse_complement("ANCGT") == "ACGNT"


def test_one_hot_encode_basic():
    arr = one_hot_encode("ACGT")
    assert arr.shape == (4, 4)
    assert arr.dtype == np.float32
    expected = np.eye(4, dtype=np.float32)
    np.testing.assert_array_equal(arr, expected)


def test_one_hot_encode_uniform_for_n():
    arr = one_hot_encode("N")
    np.testing.assert_array_equal(arr, np.full((1, 4), 0.25, dtype=np.float32))


def test_one_hot_encode_mixed():
    arr = one_hot_encode("ANT")
    expected = np.array(
        [[1, 0, 0, 0], [0.25, 0.25, 0.25, 0.25], [0, 0, 0, 1]], dtype=np.float32
    )
    np.testing.assert_array_equal(arr, expected)
