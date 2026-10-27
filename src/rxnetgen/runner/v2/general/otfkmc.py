from typing import override

from ase import Atoms

from rxnetgen.runner.general.otfkmc import OTFKMC as _OTFKMC
from rxnetgen.runner.v2.general._expl import ExplorationBase
from graphatoms.system.system import System


class OTFKMC(ExplorationBase, _OTFKMC):
    @override
    def explore(
        self,
        system: System | Atoms | None = None,
    ) -> tuple[bool, float]:
        return ExplorationBase.explore(self, system)
