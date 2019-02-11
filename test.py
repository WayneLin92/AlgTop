import unittest
from spec.specseq import SpecSeq, MyPoly
from spec.resolution import Ext
from notations import *


class OperationsTestCase(unittest.TestCase):
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
    def test_serre(self):
        MyPoly.set_prime(2)
        spec = SpecSeq(20, 16)
        x = spec.add_free_gen("x", (1, 0))
        y = spec.add_free_gen("y", (1, 0))
        w = spec.add_free_gen("w", (0, 1))
        spec.deduce_init()
        poly = spec.add_diff(w, x * x + x * y + y * y)
        spec.deduce_new_diff(poly)
        spec.draw()
        self.assertTrue(True)

    def test_adams(self):
        MyPoly.set_prime(2)
        spec = SpecSeq(13, 13, "Adams")
        spec.add_free_gen("h_{10}", (1, 1))
        spec.add_free_gen("h_{11}", (1, 2))
        spec.add_free_gen("h_{12}", (1, 4))
        spec.add_free_gen("h_{13}", (1, 8))
        spec.add_free_gen("h_{20}", (1, 3))
        spec.add_free_gen("h_{21}", (1, 6))
        spec.add_free_gen("h_{22}", (1, 12))
        spec.add_free_gen("h_{30}", (1, 7))
        spec.deduce_init()
        spec.draw()
        self.assertTrue(True)


class PolyTestCase(unittest.TestCase):
    def test_poincare_series(self):
        a = [1, 3]
        b = [1, 2, 4, 8, 3, 6]
        f = poincare_series(a, [], 10)
        g = poincare_series([], b, 10)
        self.assertEqual(f, g)

    def test_polymod2(self):
        x = PolyAnyVarMod2.gen("x")
        y = PolyAnyVarMod2.gen("y")
        p = (x - y) ** 3
        self.assertEqual(str(p), "x^3+x^2y+xy^2+y^3")
        self.assertTrue(True)

    def test_polymodp(self):
        from algebras.polynomials import PolyAnyVarModP
        PolyAnyVarModP.set_prime(3)
        x = PolyAnyVarModP.gen("x")
        y = PolyAnyVarModP.gen("y")
        p = (x - y) ** 9
        self.assertEqual(x*(-1), x*2)
        self.assertEqual(p, x ** 9 - y ** 9)
        self.assertEqual(str(p), "x^9+2y^9")
        self.assertTrue(True)

    def test_polyZ(self):
        x = PolyAnyVarZ.gen("x")
        y = PolyAnyVarZ.gen("y")
        p = (x - y) ** 3
        self.assertEqual(p, x**3 - x**2 * y * 3 + x * y**2 * 3 - y**3)
        self.assertEqual(str(p), "x^3-3x^2y+3xy^2-y^3")
        self.assertTrue(True)

    def test_polySingZ(self):
        x = xto(1)  # type: PolySingZ
        f = x + x ** 2
        g = f.inv_function(10)
        h = g.composition(f, 10)
        i = f.composition(g, 10)
        self.assertEqual(str(h), "x")
        self.assertEqual(str(i), "x")
        self.assertEqual(f.pow(5, 7), x**5 + 5 * x**6+10 * x**7)
        self.assertEqual(f.eval(2), 6)


class ResolutionTestCase(unittest.TestCase):
    def test_ext(self):
        ext = Ext(11, 30)
        ext.compute_minimal(Steenrod)
        self.assertTrue(True)

    def test_ext1(self):
        ext = Ext(15, 47)
        ext.compute_minimal(Steenrod)
        import algebras.BaseAlgebras as BC
        BC.Monitor.present()
        self.assertTrue(True)

    def test_subring_resolution(self):
        # from work.work_Lrs import MyPoly
        # x = MyPoly.gen
        ring = SubRing(30)
        ring.generate_non_comm([Sq(1), Sq(2)])
        ext = Ext(15, 30)
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
        vs = VectorSpaceMod2()
        vs.add_vectors([Sq(1), Sq(1) + Sq(2)])
        a = vs.res(Sq(1) + Sq(2) + Sq(3))
        self.assertEqual(a, Sq(3))
        sub_vs = VectorSpaceMod2((Sq(2),))
        for r in vs.quotient(sub_vs).basis(Steenrod):
            print(r)

    def test_graded_vector_space_mod2(self):
        vs = GradedVectorSpaceMod2(10)
        vs.add_vectors([Sq(10), Sq(9) * Sq(1)], 10)
        a = vs.res(Sq(10) + Sq(8) * Sq(2))
        self.assertEqual(a, Sq(8) * Sq(2))

    def test_linear_map_mod2(self):
        lin_map = LinearMapMod2()
        lin_map.add_maps([(Sq(1), Sq(1)), (Sq(2), Sq(2)), (Sq(3), Sq(1)+Sq(2))])
        print("\n")
        lin_map.present_kernel(Steenrod)
        print(lin_map.f(Sq(1) + Sq(3)))
        print(lin_map.g(Sq(1) + Sq(2)))
        self.assertTrue(True)

    def test_lin_map_kernel_mod2(self):
        lin_map = LinearMapKernelMod2()
        lin_map.add_maps([(Sq(1), Sq(1)), (Sq(2), Sq(2)), (Sq(3), Sq(1)+Sq(2))])
        print("\nkernel:\n")
        for r in lin_map.kernel.basis(Steenrod):
            print(r)
        print("\nimage:\n")
        for r in lin_map.image().basis(Steenrod):
            print(r)
        self.assertTrue(True)


class ConstructionsTestCase(unittest.TestCase):
    def test_subring(self):
        self.assertTrue(True)

    # def test_quo_ring(self):
    #     x = MyQuoRing.gen
    #     MyQuoRing.init(10)
    #     rels = [x(i) * x(i) for i in range(1, 6)]
    #     MyQuoRing.add_relations(rels)
    #     for r in MyQuoRing.basis(7):
    #         print(r)
    #     print(x(1) * (x(2) + x(1) * x(1)))
    #     self.assertTrue(True)


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


if __name__ == '__main__':
    unittest.main()
