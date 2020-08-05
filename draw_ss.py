"""Draw spectral sequences"""
import argparse
from algebras import mymath
from MaySS.groebner import GbAlgMod2, GbDga
from GUI.specseq import SpecSeq

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('--pred', nargs=2, type=int, help="maximum stem and s to display")
parser.add_argument('-g', action="store_true", help="display only generators")
parser.add_argument('--dga', action="store_true", help="indicate this is a DGA")
parser.add_argument('--arrows', action="store_true", help="draw arrows for dga")
parser.add_argument('--h0', action="store_true", help="indicate this is a localization")
parser.add_argument('path', help="the path of the pickle file to be loaded")
args = parser.parse_args()
# print(args)
# quit()


def pred_graph(d3d):
    return d3d[1] - d3d[0] <= args.pred[0] and d3d[0] <= args.pred[1]


def pred_500(d3d):
    return d3d[1] - d3d[0] <= 500


def pred_E4_149_86(d3d):
    return d3d[1] - d3d[0] <= 149 and d3d[0] <= 86


def pred_118_59(d3d):
    return d3d[1] - d3d[0] <= 118 and d3d[0] <= 59


def pred_50_25(d3d):
    return d3d[1] - d3d[0] <= 50 and d3d[0] <= 25


def pred_200_100(d3d):
    return d3d[1] - d3d[0] <= 200 and d3d[0] <= 100


def pred_256_128(d3d):
    return d3d[1] - d3d[0] <= 256 and d3d[0] <= 128


def pred_152_76(d3d):
    return d3d[1] - d3d[0] <= 152 and d3d[0] <= 76


def pred_71_36(d3d):
    return d3d[1] - d3d[0] <= 71 and d3d[0] <= 36


def pred_70_35(d3d):
    return d3d[1] - d3d[0] <= 70 and d3d[0] <= 35


def pred_124_62(d3d):
    return d3d[1] - d3d[0] <= 124 and d3d[0] <= 62


def pred_123_61(d3d):
    return d3d[1] - d3d[0] <= 123 and d3d[0] <= 61


def pred_200(d3d):
    return d3d[1] <= 200


if __name__ == "__main__":
    directory = "C:\\Users\\lwnpk\\Documents\\MyProgramData\\Math_AlgTop\\pickle\\"
    if args.dga:
        Alg = GbDga.load_alg(directory + args.path)
    else:
        Alg = GbAlgMod2.load_alg(directory + args.path)
    pred = pred_graph if args.pred else Alg.pred
    spec = SpecSeq()
    if args.g:  # Only draw generators
        for gen in Alg.generators:
            s, t, u = gen.deg3d
            x = t - s
            y = u - s if args.h0 else s
            if args.dga and gen.diff is None:
                spec.add_bullet((x, y), gen.name, [0, 0, 255])
            else:
                spec.add_bullet((x, y), gen.name, [255, 0, 0])
    else:
        basis_alg = Alg.basis_mons(pred=pred)
        for d, basis_d in basis_alg.items():
            s, t, u = d
            x = t - s
            y = u - s if args.h0 else s
            for m in basis_d:
                if sum(e for i, e in m) == -1:
                    if args.dga:
                        g = Alg.get_gen(mymath.get_one_element(m)[0])
                        if g.diff is not None:
                            spec.add_bullet((x, y), Alg.str_mon(m), [255, 0, 0])
                        else:
                            spec.add_bullet((x, y), Alg.str_mon(m), [0, 0, 255])
                    else:
                        spec.add_bullet((x, y), Alg.str_mon(m), [255, 0, 0])
                else:
                    spec.add_bullet((x, y), Alg.str_mon(m), [0, 0, 0])
        for d, basis_d in basis_alg.items():
            for m in basis_d:
                if sum(e for i, e in m) < -1:
                    m1 = ((m[0][0], m[0][1] + 1),) + m[1:] if m[0][1] + 1 else m[1:]
                    spec.add_line(Alg.str_mon(m1), Alg.str_mon(m), label=Alg.get_gen(m[0][0]).name,
                                  color=[100, 100, 100], by_label=True)
                for gen in Alg.generators[:3]:
                    if pred(d + gen.deg3d):
                        m1 = Alg(m) * Alg(((gen.index, -1),))
                        m1 = m1.data
                        if len(m1) == 1:
                            m1 = mymath.get_one_element(m1)
                            if m1 and m1[0][0] != gen.index:
                                spec.add_line(Alg.str_mon(m), Alg.str_mon(m1), label=Alg.get_gen(m[0][0]).name,
                                              color=[100, 100, 100], by_label=True)
                if args.dga and args.arrows and pred(d + Alg.deg_diff) and Alg.is_differentiable_mon(m):
                    m1 = Alg(m).diff()
                    m1 = m1.data
                    if len(m1) == 1:
                        m1 = mymath.get_one_element(m1)
                        spec.add_arrow(Alg.str_mon(m), Alg.str_mon(m1), by_label=True)
    spec.draw(args.path)
