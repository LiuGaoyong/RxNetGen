from typing import override

from ase import Atoms

from rxnetgen.runner.general.rxngen import (
    ReactionNetworkGenerator as _Rxngen,
)
from rxnetgen.runner.v2.general._expl import ExplorationBase
from graphatoms.system.system import System


class ReactionNetworkGenerator(ExplorationBase, _Rxngen):
    @override
    def explore(
        self,
        system: System | Atoms | None = None,
    ) -> tuple[bool, float]:
        return ExplorationBase.explore(self, system)
