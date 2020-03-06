import unittest
from typing import Union


class OperationsTestCase(unittest.TestCase):
    def setUp(self):
        pass

    def test_Sq(self):
        a = Sq(4) * Sq(4) * Sq(4)
        b = Sq(8) * Sq(3) * Sq(1) + Sq(9) * Sq(2) * Sq(1) + Sq(10) * Sq(2) + Sq(11) * Sq(1)
        self.assertTrue(a)
        self.assertEqual(a, b)

        c = Sq(8) * Sq(1) + Sq(6) * Sq(2) * Sq(1)
        self.assertEqual(conj_Sq(9), c)

    def test_Q(self):
        a = Q(20) * Q(8)
        b = Q(18) * Q(10) + Q(17) * Q(11)
        self.assertEqual(a, b)

    def test_L(self):
        k_max = 20
        a = AR.sQ(3, k_max) * AR.sQ(12, k_max)
        a = a.truncate(0)
        b = Q(9) * Q(6) + Q(8) * Q(7)
        a = str(a)
        b = str(b)
        self.assertEqual(a, b)

    def test_xi(self):
        self.assertTrue(len(str(psi_xi(4))) == 92)

    def test_dyerlashofx(self):
        from algebras.operations import DyerLashofX
        x = DyerLashofX.gen()
        DyerLashofX.set_x_deg(2)
        y5 = Q(3) * x
        y7 = Q(5) * x
        y9 = Q(7) * x
        y13 = Q(11) * x
        y8 = Q(6) * x + x ** 4
        y10 = Q(8) * x + x ** 2 * (Q(4) * x)
        y12 = Q(10) * x + (Q(4) * x) ** 2
        print(y12)
        z = Q(20) * y10 + Q(18) * y12 + Q(17) * y13 + x ** 4 * (Q(12) * y10) + \
            y9 ** 2 * (Q(4) * x) ** 2 + y7 ** 2 * (Q(9) * Q(5) * x) + y8 ** 2 * (Q(8) * Q(4) * x) + \
            Q(9) * y9 * (Q(4) * x) ** 2 + Q(10) * y8 * (Q(4) * x) ** 2 + \
            y5 ** 2 * (Q(11) * Q(7) * x + Q(10) * Q(8) * x + x ** 4 * (Q(6) * Q(4) * x))
        b = DyerLashofX.zero()
        self.assertEqual(z, b)

    def test_pair(self):
        milnor1 = Sq(1)
        milnor2 = Sq(2) * milnor1 - milnor1 * Sq(2)
        milnor3 = Sq(4) * milnor2 - milnor2 * Sq(4)
        milnor4 = Sq(8) * milnor3 - milnor3 * Sq(8)
        is_first = True
        for poly_xi in DualSteenrod.basis(15):
            a = milnor4.pair(poly_xi)
            b = poly_xi.pair(milnor4)
            if is_first:
                self.assertEqual(a, 1)
                is_first = False
            else:
                self.assertNotEqual(a, 1)
                self.assertEqual(a, b)

    def test_L_work(self):
        from algebras.operations import AR
        from algebras.mymath import choose_mod2
        k_max = 20
        r, s = 3, 8
        lhs = Sq(r) * AR.sQ(s, k_max)
        rhs = sum((AR.sQ(s - j, k_max) * Sq(r - j) for j in range(r + 1) if choose_mod2(s - j, j)), AR.zero())
        self.assertEqual(lhs.truncate(k_max), rhs.truncate(k_max))

    def test_monomial_xi(self):
        for mon in DualSteenrod.basis(0):
            print(mon)
        self.assertTrue(True)

    def test_T2(self):
        pass


class BUTestCase(unittest.TestCase):
    def test_bu(self):
        from algebras.BU import b_rd, chi_b
        for n in range(5, 8):
            a = b_rd(1) % b_rd(n)
            b = b_rd(1) % chi_b[n]
            self.assertEqual(a, b)


class SpecSeqTestCase(unittest.TestCase):
    def setUp(self):
        from GUI.specseq import SpecSeq
        self.SpecSeq = SpecSeq

    def test_serre(self):
        spec = self.SpecSeq(10, 10)
        x = spec.add_gen("x", (1, 0))
        y = spec.add_gen("y", (1, 0))
        w = spec.add_gen("w", (0, 1))
        spec.init_diff()
        spec.add_diff(w, x * x + x * y + y * y)
        spec.draw()
        self.assertTrue(True)

    def test_adams(self):
        spec = self.SpecSeq(13, 13)
        spec.add_gen("h_{10}", (1, 1))
        spec.add_gen("h_{11}", (1, 2))
        spec.add_gen("h_{12}", (1, 4))
        spec.add_gen("h_{13}", (1, 8))
        spec.add_gen("h_{20}", (1, 3))
        spec.add_gen("h_{21}", (1, 6))
        spec.add_gen("h_{22}", (1, 12))
        spec.add_gen("h_{30}", (1, 7))
        spec.init_diff()
        spec.draw()
        self.assertTrue(True)


class PolyTestCase(unittest.TestCase):
    def setUp(self):
        from algebras import polynomials
        self.polynomials = polynomials

    def test_poincare_series(self):
        poincare_series = self.polynomials.poincare_series
        a = [1, 3]
        b = [1, 2, 4, 8, 3, 6]
        f = poincare_series(a, [], 10)
        g = poincare_series([], b, 10)
        self.assertEqual(f, g)

    def test_polymod2(self):
        PolyAnyVarMod2 = self.polynomials.PolyAnyVarMod2
        x = PolyAnyVarMod2.gen("x")
        y = PolyAnyVarMod2.gen("y")
        p = (x - y) ** 3
        self.assertEqual("x^3 + x^2y + xy^2 + y^3", str(p))
        self.assertTrue(True)

    def test_polymodp(self):
        PolyAnyVarModP = self.polynomials.PolyAnyVarModP
        PolyAnyVarModP.set_prime(3)
        x = PolyAnyVarModP.gen("x")
        y = PolyAnyVarModP.gen("y")
        p = (x - y) ** 9
        self.assertEqual(x*(-1), x*2)
        self.assertEqual(p, x ** 9 - y ** 9)
        self.assertEqual(str(p), "x^9+2y^9")
        self.assertTrue(True)

    def test_polyZ(self):
        PolyAnyVarZ = self.polynomials.PolyAnyVarZ
        x = PolyAnyVarZ.gen("x")
        y = PolyAnyVarZ.gen("y")
        p = (x - y) ** 3
        self.assertEqual(x**3 - x**2 * y * 3 + x * y**2 * 3 - y**3, p)
        self.assertEqual("x^3-3x^2y+3xy^2-y^3", str(p))
        self.assertTrue(True)

    def test_polySingZ(self):
        PolySingZ = self.polynomials.PolySingZ
        x = PolySingZ.gen(1)
        f = x + x ** 2
        g = f.inverse_composition(10)
        h = g.composition(f, 10)
        i = f.composition(g, 10)
        self.assertEqual(str(h), "x")
        self.assertEqual(str(i), "x")
        self.assertEqual(f.pow(5, 7), x**5 + 5 * x**6+10 * x**7)
        self.assertEqual(f.evaluation(2), 6)


class ResolutionTestCase(unittest.TestCase):
    def setUp(self):
        from algebras.resolution import Ext
        self.Ext = Ext

    def test_ext(self):
        from algebras.operations import Steenrod
        ext = self.Ext(11, 30)
        ext.compute_minimal(Steenrod)
        self.assertTrue(True)

    def test_ext1(self):
        ext = self.Ext(15, 47)
        ext.compute_minimal(Steenrod)
        import algebras.BaseAlgebras as BC
        BC.Monitor.present()
        self.assertTrue(True)

    def test_subring_resolution(self):
        # from work.work_Lrs import MyPoly
        # x = MyPoly.gen
        ring = SubRing(30)
        ring.generate_non_comm([Sq(1), Sq(2)])
        ext = self.Ext(15, 30)
        ext.compute_minimal(ring)
        spec = ext.get_spec(20, 15)
        spec.draw()
        self.assertTrue(True)

    def test_amod(self):
        sq = Steenrod.gen
        FreeModule.set_ring(Steenrod)
        FreeModuleMod2.set_ring(Steenrod)
        x = FreeModule.gen('x')
        y = FreeModuleMod2.gen('y')
        # print(sq(4) * sq(4) * y)
        self.assertEqual(sq(4) * sq(4) * x, sq(4) * (sq(4) * x))
        self.assertEqual(sq(4) * sq(4) * y, sq(4) * (sq(4) * y))


class LinAlgTestCase(unittest.TestCase):
    def test_vector_space_mod2(self):
        vs = VectorSpaceMod2([Sq(1), Sq(1) + Sq(2)])
        a = vs.res(Sq(1) + Sq(2) + Sq(3))
        self.assertEqual(a, Sq(3))
        sub_vs = VectorSpaceMod2((Sq(2),))
        for r in (vs / sub_vs).basis(Steenrod):
            print(r)

    def test_graded_vector_space_mod2(self):
        vs = GradedVectorSpaceMod2()
        vs.add_vectors([Sq(10), Sq(9) * Sq(1)])
        a = vs.res(Sq(10) + Sq(8) * Sq(2))
        self.assertEqual(a, Sq(8) * Sq(2))

    def test_linear_map_mod2(self):
        lin_map = LinearMapMod2()
        lin_map.add_maps([(Sq(1), Sq(1)), (Sq(2), Sq(2)), (Sq(3), Sq(1)+Sq(2))])
        print("\n")
        for r in lin_map.kernel.basis(Steenrod):
            print(r)
        print(lin_map.f(Sq(1) + Sq(3)))
        print(lin_map.g(Sq(1) + Sq(2)))
        self.assertTrue(True)

    def test_lin_map_kernel_mod2(self):
        lin_map = LinearMapKMod2()
        lin_map.add_maps([(Sq(1), Sq(1)), (Sq(2), Sq(2)), (Sq(3), Sq(1)+Sq(2))])
        print("\nkernel:\n")
        for r in lin_map.kernel.basis(Steenrod):
            print(r)
        print("\nimage:\n")
        for r in lin_map.image.basis(Steenrod):
            print(r)
        self.assertTrue(True)


class MyDyerLashofTestCase(unittest.TestCase):
    def test_sq_benchmark(self):
        sQ = MyDyerLashof.gen
        prod = sQ(16) * sQ(34) * sQ(65) * sQ(66)
        print(prod.simplify())
        self.assertTrue(True)

    def test_simplify_sQ_benchmark(self):
        result = simplify_sQ((16, 34, 65, 66))
        print(result)
        self.assertTrue(True)

    def test(self):
        sQ = MyDyerLashof.gen
        prod = sQ(16) * sQ(34) * sQ(65) * sQ(66)
        prod.simplify()
        result = simplify_sQ((16, 34, 65, 66))
        self.assertEqual(prod.data, result)


class MymathTestCase(unittest.TestCase):
    def setUp(self):
        import algebras.mymath
        self.mymath = algebras.mymath

    def test_orderedpartition(self):
        ls = list(self.mymath.orderedpartition(4, 30))
        answer = 3654
        self.assertEqual(answer, len(ls))


class Benchmark(unittest.TestCase):
    def test_alg_B(self):
        pass


class GroebnerTestCase(unittest.TestCase):
    def setUp(self):
        from algebras.groebner import GbAlgMod2
        from itertools import permutations
        self.GbAlgMod2 = GbAlgMod2
        self.permutations = permutations

    def alg_B(self, n_max):
        gens = []
        for i in range(n_max):
            for j in range(i + 1, n_max + 1):
                gens.append((f"R_{{{i}{j}}}", 2 ** j - 2 ** i, j - i))
        gens.sort(key=lambda _x: _x[2])

        E1 = self.GbAlgMod2.new_alg(key=lambda _m: [-_i for _i in _m])
        E1.add_gens(gens)

        def R(S: Union[int, tuple], T: Union[int, tuple]):
            if type(S) is int:
                return E1.gen(f"R_{{{S}{T}}}")
            assert len(S) == len(T)
            S, T = sorted(S), sorted(T)
            result = E1.zero()
            for T1 in self.permutations(T):
                if all(t - s > 0 for s, t in zip(S, T1)):
                    pro = E1.unit()
                    for x in map(R, S, T1):
                        pro *= x
                    result += pro
            return result

        rels = []
        for d in range(2, n_max + 1):
            for i in range(n_max + 1 - d):
                j = i + d
                rel = sum((R(i, k) * R(k, j) for k in range(i + 1, j)), E1.zero())
                rels.append(rel)
        rels.sort(key=lambda x: x.deg())
        E1.add_rels(rels, sorted_=True, clear_cache=True)
        return E1, R

    def test_GbAlgMod2(self):
        E1, _ = self.alg_B(7)
        self.assertEqual(65, len(E1._rels))

    def test_subalgebra(self):
        n_max = 4
        E1, R = self.alg_B(n_max)
        ele_names = []
        for i in range(n_max):
            ele_names.append((R(i, i + 1), f'h_{i}'))
        for i in range(n_max - 2):
            ele_names.append((R((i, i + 1), (i + 2, i + 3)), f'h_{i}(1)'))
        for i in range(n_max - 4):
            ele_names.append((R((i, i + 1, i + 3), (i + 2, i + 4, i + 5)), f'h_{i}(1, 3)'))
        for i in range(n_max - 4):
            ele_names.append((R((i, i + 1, i + 2), (i + 3, i + 4, i + 5)), f'h_{i}(1, 2)'))
        for d in range(2, n_max + 1):
            for i in range(n_max + 1 - d):
                j = i + d
                ele_names.append((R(i, j) * R(i, j), f'b_{{{i}{j}}}'))
        HX = E1.subalgebra(ele_names, key=E1._key)
        self.assertEqual(15, len(HX._rels))

    @classmethod
    def tearDownClass(cls):
        from algebras.BaseAlgebras import Monitor
        Monitor.present()


if __name__ == '__main__':
    unittest.main()
