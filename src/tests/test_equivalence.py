# ruff: noqa: D100, D101, D102, D103, D104
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest
from ase import Atoms
from ase.cluster import Octahedron

from graphatoms.system import Cluster, System
from rxnetgen.config import GeometryConfig
from rxnetgen.network import ReactionNetwork
from rxnetgen.network._equivalence import EquivalenceRegistry
from rxnetgen.runner.common import _helper
from rxnetgen.runner.common._base import _unique_clusters_by_itlee
from rxnetgen.runner.common._helper import record_reconstruction


@pytest.fixture(scope="module")
def system() -> System:
    return System.from_ase(
        Octahedron("Pd", 8),
        parse_atoms_is_outer_or_not=True,
    )


@pytest.fixture(scope="module")
def cluster(system: System) -> Cluster:
    return Cluster.from_select(
        system,
        np.array([0]),
        method="hop",
        max_moved_threshold=1,
        env_threshold=2,
    )


class Test_EquivalenceRegistry:
    def test_persistence_and_halving(self, tmp_path: Path) -> None:
        registry = EquivalenceRegistry()
        registry.set_delta("cluster-a", 0.2)
        registry.record_attempt("cluster-a")
        registry.record_failure("cluster-a")
        assert registry.false_reuse_rate("cluster-a") == pytest.approx(0.5)
        new_delta = registry.halve("cluster-a")
        assert new_delta == pytest.approx(0.1)

        path = tmp_path / "equivalence.json"
        registry.write_json(path)
        loaded = EquivalenceRegistry.read_json(path)
        assert loaded.get_delta("cluster-a", 0.0) == pytest.approx(0.1)
        assert loaded.clusters["cluster-a"].attempts == 2
        assert loaded.clusters["cluster-a"].failures == 1

    def test_missing_file(self, tmp_path: Path) -> None:
        registry = EquivalenceRegistry.read_json(tmp_path / "missing.json")
        assert registry.clusters == {}
        assert registry.events == {}

    def test_network_persistence(self, tmp_path: Path) -> None:
        path = tmp_path / "event"
        network = ReactionNetwork(path, restart=False)
        network.equivalence.set_delta("cluster-b", 0.2)
        network.persistence()
        loaded = ReactionNetwork(path, restart=True)
        assert loaded.equivalence.get_delta("cluster-b", 0.0) == pytest.approx(
            0.2
        )


class Test_ClusterDedup:
    def test_merge_same_geometry(self, cluster: Cluster) -> None:
        dct = cluster.to_dict()
        dct.pop("hashes", None)
        dct["coordination"] = (cluster.CN_MATRIX + 1).astype("uint8")
        other = Cluster.model_validate(dct)
        values = [cluster, other]
        keys = [
            (True, 1, cluster.hash),
            (True, 1, other.hash),
        ]
        registry = EquivalenceRegistry()
        selected = _unique_clusters_by_itlee(
            values=values,
            keys=keys,
            candidates=np.array([0, 1]),
            tolerance=0.05,
            equivalence=registry,
        )
        assert len(selected) == 1
        assert int(selected[0]) == 0

    def test_role_protection(self, cluster: Cluster) -> None:
        dct = cluster.to_dict()
        dct.pop("hashes", None)
        dct["coordination"] = (cluster.CN_MATRIX + 1).astype("uint8")
        dct["is_core"] = np.zeros_like(cluster.is_core, dtype=bool)
        other = Cluster.model_validate(dct)
        values = [cluster, other]
        keys = [
            (True, 1, cluster.hash),
            (True, 1, other.hash),
        ]
        registry = EquivalenceRegistry()
        selected = _unique_clusters_by_itlee(
            values=values,
            keys=keys,
            candidates=np.array([0, 1]),
            tolerance=0.05,
            equivalence=registry,
        )
        assert len(selected) == 2


class _StubNetwork:
    def __init__(self) -> None:
        self.equivalence = EquivalenceRegistry()
        self.persisted = 0

    def persistence_equivalence(self) -> None:
        self.persisted += 1


def test_record_reconstruction() -> None:
    class _Config:
        geometry = GeometryConfig(
            tolerance=0.1,
            reconstruction_check=True,
            delta_halving=True,
        )

    network = _StubNetwork()
    ok, delta, limit = record_reconstruction(
        config=_Config(),
        network=network,
        for_cluster="c1",
        rmsd=0.001,
        n_atoms=100,
    )
    assert ok
    assert delta == pytest.approx(0.1)
    assert limit == pytest.approx(0.01)
    assert network.equivalence.clusters["c1"].attempts == 1
    assert network.equivalence.clusters["c1"].failures == 0

    ok, _, _ = record_reconstruction(
        config=_Config(),
        network=network,
        for_cluster="c1",
        rmsd=0.02,
        n_atoms=100,
    )
    assert not ok
    assert network.equivalence.clusters["c1"].failures == 1
    assert network.equivalence.get_delta("c1", 0.0) == pytest.approx(0.05)
    assert network.persisted == 2


def test_geometry_config_defaults() -> None:
    config = GeometryConfig()
    assert config.tolerance == 0.0
    assert config.allow_reflection is False
    assert config.max_branches > 0


class _FakeInfo:
    for_cluster = "cluster-test"
    key_r = "r"
    key_p = "p"


class _FakeRxn:
    R = SimpleNamespace()

    def __init__(self) -> None:
        self.calls: list[np.ndarray | None] = []

    def apply(self, system, matched_indxs=None):
        self.calls.append(matched_indxs)
        return Atoms("H2", positions=np.zeros((2, 3))), 0.01


class _FakeSystem:
    def __init__(self, mapping=None) -> None:
        self.mapping = mapping
        self.kwargs = None

    def get_match_mode(self, pattern, **kwargs):
        self.kwargs = kwargs
        return self.mapping


class _FakeNetwork:
    def __init__(self, rxn) -> None:
        self.rxn = rxn
        self.equivalence = EquivalenceRegistry()

    def read_event_info(self, key):
        return _FakeInfo()

    def read_event(self, key):
        return _FakeInfo(), self.rxn

    def read_sysgraph(self, key, type, return_sysgraph_for_minima_or_ts):
        return SimpleNamespace()


def test_helper_apply_tolerant_match(monkeypatch) -> None:
    rxn = _FakeRxn()
    network = _FakeNetwork(rxn)
    mapping = np.array([[0, 1]])
    system = _FakeSystem(mapping)
    monkeypatch.setattr(
        _helper.ReactionNetwork,
        "from_config",
        lambda *args, **kwargs: network,
    )
    monkeypatch.setattr(_helper, "System", _FakeSystem)
    config = SimpleNamespace(geometry=GeometryConfig(tolerance=0.1))

    result = _helper.helper_apply(
        rxn_key="rxn-1",
        config=config,
        system=system,
        match_mode=None,
        raise_when_fail=True,
    )
    assert result[0] == "rxn-1"
    assert isinstance(result[2], Atoms)
    assert rxn.calls == [mapping]
    assert system.kwargs is not None
    assert system.kwargs["tolerance"] == pytest.approx(0.1)


def test_helper_apply_no_tolerant_match(monkeypatch) -> None:
    rxn = _FakeRxn()
    network = _FakeNetwork(rxn)
    system = _FakeSystem(None)
    monkeypatch.setattr(
        _helper.ReactionNetwork,
        "from_config",
        lambda *args, **kwargs: network,
    )
    monkeypatch.setattr(_helper, "System", _FakeSystem)
    config = SimpleNamespace(geometry=GeometryConfig(tolerance=0.1))

    result = _helper.helper_apply(
        rxn_key="rxn-1",
        config=config,
        system=system,
        match_mode=None,
        raise_when_fail=True,
    )
    assert isinstance(result[2], str)
    assert rxn.calls == []


def test_helper_apply_exact_fallback(monkeypatch) -> None:
    rxn = _FakeRxn()
    network = _FakeNetwork(rxn)
    system = _FakeSystem(None)
    monkeypatch.setattr(
        _helper.ReactionNetwork,
        "from_config",
        lambda *args, **kwargs: network,
    )
    monkeypatch.setattr(_helper, "System", _FakeSystem)
    config = SimpleNamespace(geometry=GeometryConfig(tolerance=0.0))

    result = _helper.helper_apply(
        rxn_key="rxn-1",
        config=config,
        system=system,
        match_mode=None,
        raise_when_fail=True,
    )
    assert isinstance(result[2], Atoms)
    assert rxn.calls == [None]


def test_helper_match_tolerance(monkeypatch) -> None:
    rxn = _FakeRxn()
    network = _FakeNetwork(rxn)
    network.equivalence.set_delta("cluster-test", 0.05)
    mapping = np.array([[0, 1]])
    system = _FakeSystem(mapping)
    monkeypatch.setattr(
        _helper.ReactionNetwork,
        "from_config",
        lambda *args, **kwargs: network,
    )
    config = SimpleNamespace(geometry=GeometryConfig(tolerance=0.1))

    result = _helper.helper_match(
        config=config,
        system=system,
        rxn_is_forward=True,
        rxn_key="rxn-1",
    )
    assert result[0] == "rxn-1"
    assert result[1] is True
    assert np.array_equal(result[2], mapping)
    assert system.kwargs is not None
    assert system.kwargs["tolerance"] == pytest.approx(0.05)
