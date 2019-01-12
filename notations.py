from algebras.operations import Steenrod, DyerLashof, AR, DualSteenrod
from algebras.polynomials import PolySingZ, PolyAnyVarMod2, PolyAnyVarModP, PolyAnyVarZ, poincare_series
from algebras.constructions import SubRing, QuoRing, FreeModule, FreeModuleMod2
from algebras.mymath import choose_mod2, binom_mod2, multi_nom_mod2
from algebras.linalg import VectorSpaceMod2, GradedVectorSpaceMod2, LinearMapMod2, LinearMapKernelMod2
from algebras.homology import HBOZ
from algebras.my_dyerlashof import MyDyerLashof, simplify_sQ


Q = DyerLashof.gen
Sq = Steenrod.gen
conj_Sq = Steenrod.conj_gen
psi_Sq = Steenrod.coprod_gen
xi = DualSteenrod.gen
psi_xi = DualSteenrod.coprod_gen
xto = PolySingZ.gen


pre_defined = \
    [choose_mod2, binom_mod2, multi_nom_mod2,
     Steenrod, DyerLashof, AR, DualSteenrod,
     PolySingZ, PolyAnyVarMod2, PolyAnyVarModP, PolyAnyVarZ, poincare_series,
     SubRing, QuoRing, FreeModule, FreeModuleMod2,
     VectorSpaceMod2, GradedVectorSpaceMod2, LinearMapMod2, LinearMapKernelMod2,
     Q, Sq, conj_Sq, psi_Sq, xi, psi_xi, xto,
     HBOZ,
     MyDyerLashof, simplify_sQ
     ]
