"""doc"""
import argparse
import os
import subprocess
import sqlite3
import textwrap
from collections import defaultdict
from MaySS.groebner import GbAlgMod2

def str_to_poly(s, R):
    if not s:
        return R.zero()
    data = set()
    for str_mon in s.split(";"):
        mon = tuple(map(int, str_mon.split(",")))
        it_zip = zip((it := iter(mon)), it)
        mon = tuple((i, -e) for i, e in it_zip)
        data.add(mon)
    return R(data)

def str_to_mon(s):
    if not s:
        return tuple()
    mon = tuple(map(int, s.split(",")))
    it_zip = zip((it := iter(mon)), it)
    return tuple((i, -e) for i, e in it_zip)

def str_to_array(s):
    if not s:
        return tuple()
    a = tuple(map(int, s.split(",")))
    return a

if __name__ == "__main__":
    # parser
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--edit', action='store_true', help='open the script in vscode')
    parser.add_argument('--table', help="sql table name")
    parser.add_argument('--stem_max', type=int, default=30, help='truncate ss by t<=t_max')
    args = parser.parse_args()
    if args.edit:
        subprocess.Popen(f"code {__file__}", shell=True)
        os.sys.exit()

    # actions
    data_dir = R"C:\Users\lwnpk\Documents\MyProgramData\Math_AlgTop\database" "\\"
    html_file_path = R"C:\Users\lwnpk\Dropbox\HTML\Weinan - Website\programs\sstable\MaySS.tex"
    kLevelMax = 10000
    conn = sqlite3.connect(data_dir + 'tmp.db')
    with conn:
        c = conn.cursor()

        # Load E4
        E4 = GbAlgMod2.new_alg()
        for row in c.execute("SELECT gen_name, s, t, v FROM E4_generators ORDER BY gen_id"):
            E4.add_gen(row[0], row[2], deg3d=(row[1], row[2], row[3]))

        # Load degrees
        print("Loading degrees...")
        degs = []
        for row in c.execute("SELECT DISTINCT s, t, v FROM E4_basis WHERE t<=? ORDER BY mon_id", (args.stem_max * 1.5,)):
            degs.append((row[0], row[1], row[2]))
        f = lambda d: (d[1] - d[0], d[0], d[2])
        degs.sort(key=f)

        # Load basis
        print("Loading basis...")
        basis = defaultdict(list)
        for deg in degs:
            print(f(deg), end="      \r")
            for row in c.execute(f"SELECT mon FROM E4_basis WHERE s=? and t=? and v=? ORDER BY mon_id", deg):
                basis[deg].append(str_to_mon(row[0]))
        
        degs1 = [d for d in degs if d[1] - d[0] <= args.stem_max and d[0] <= (d[1] - d[0]) // 2 + 1]
        # write to file
        print("Write to TEX...")
        with open(html_file_path, "w") as file:
            file.write(textwrap.dedent("""\
                \\documentclass{article}
                \\usepackage[utf8]{inputenc}
                \\usepackage[margin=2cm]{geometry}
                \\usepackage{array}
                \\usepackage{longtable}

                \\begin{document}
                \\begin{longtable}{|l|l|>{\\raggedright\\arraybackslash}p{6cm}|>{\\raggedright\\arraybackslash}p{6cm}|}\\hline
                """))
            num_row = 0
            for deg in degs1:
                print(f(deg), end="      \r")
                s, t, v = deg
                for row in c.execute("SELECT base, diff, level FROM E4_ss WHERE s=? and t=? and v=? ORDER BY base_id", deg):
                    base = str_to_array(row[0])
                    diff = str_to_array(row[1]) if row[1] is not None else None
                    level = row[2]
                    num_row += 1
                    file.write(f'{num_row} & ')
                    file.write(f'{t - s},{s};{v} & ')
                    if (level < 5000):
                        r = level
                        base1 = E4({basis[deg][i] for i in base})
                        diff1 = E4({basis[(deg[0] - 1, deg[1], deg[2] + r)][i] for i in diff}) if diff is not None else "?"
                        file.write(f'${base1}$ & ')
                        file.write(f'$d_{{{r}}}^{{-1}}={diff1}$')
                    elif (level<7500):
                        base1 = E4({basis[deg][i] for i in base})
                        file.write(f'${base1}$ & ')
                        file.write(f'Permanent cycle')
                    else:
                        r = 10000 - level
                        base1 = E4({basis[deg][i] for i in base})
                        diff1 = E4({basis[(deg[0] + 1, deg[1], deg[2] - r)][i] for i in diff}) if diff is not None else "?"
                        file.write(f'${base1}$ &')
                        file.write(f'$d_{{{r}}}={diff1}$')
                    file.write("\\\\\n")
                file.write("\\hline\n")
            file.write(textwrap.dedent("""\
                \end{longtable}
                \end{document}
                """))
        c.close()
    conn.close()