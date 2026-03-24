"""Tests for wifi-ndf."""

import numpy as np
import pytest
from wifi_ndf import compute_ndf


def test_basic_4_nodes():
    result = compute_ndf([(0, 0), (4, 0), (4, 3), (0, 3)], room=(5, 4))
    assert result.ndf > 0
    assert result.num_nodes == 4
    assert result.num_links == 6
    assert result.efficiency > 0.9  # few nodes = high efficiency


def test_14_nodes_reference_room():
    np.random.seed(42)
    nodes = [(np.random.uniform(0.3, 4.7), np.random.uniform(0.3, 3.7)) for _ in range(14)]
    result = compute_ndf(nodes, room=(5.0, 4.0), freq_ghz=2.4)
    assert result.ndf > 0
    assert result.num_links == 91
    assert 0.01 <= result.threshold <= 0.1


def test_collinear_low_ndf():
    nodes = [(i * 0.5, 2.0) for i in range(14)]
    result = compute_ndf(nodes, room=(7.0, 4.0))
    assert result.efficiency < 0.5  # collinear = low efficiency


def test_efficiency_well_spread():
    nodes = [(0.3, 0.3), (4.7, 0.3), (4.7, 3.7), (0.3, 3.7),
             (2.5, 0.3), (2.5, 3.7), (0.3, 2.0), (4.7, 2.0)]
    result = compute_ndf(nodes, room=(5.0, 4.0))
    assert result.efficiency > 0.8


def test_different_frequencies():
    nodes = [(0, 0), (3, 0), (3, 3), (0, 3)]
    r24 = compute_ndf(nodes, freq_ghz=2.4)
    r50 = compute_ndf(nodes, freq_ghz=5.0)
    assert r24.ndf > 0
    assert r50.ndf > 0


def test_auto_room_bounds():
    nodes = [(1, 1), (4, 1), (4, 3), (1, 3)]
    result = compute_ndf(nodes)  # no room= argument
    assert result.ndf > 0


def test_diagnose():
    nodes = [(0, 0), (4, 0), (4, 3), (0, 3)]
    result = compute_ndf(nodes, room=(5, 4))
    assert isinstance(result.diagnose(), str)
    assert len(result.diagnose()) > 0


def test_summary():
    nodes = [(0, 0), (4, 0), (4, 3), (0, 3)]
    result = compute_ndf(nodes, room=(5, 4))
    s = result.summary()
    assert "NDF=" in s
    assert "efficiency=" in s


def test_too_few_nodes():
    with pytest.raises(ValueError):
        compute_ndf([(0, 0)])


def test_wrong_shape():
    with pytest.raises(ValueError):
        compute_ndf([(0, 0, 0), (1, 1, 1)])
