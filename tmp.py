def aij_matrix(n, m):
    result = ""
    for i in range(1, n + 1):
        result += " & ".join(f"a_{{{i}{j}}}" for j in range(1, m + 1))
        result += "\\\\\n"
    print(result)


def texify_matrix(text: str):
    result = "\\begin{bmatrix}\n"
    rows = text.split(sep="\n")
    for row in rows:
        r = row.split()
        if r:
            result += " & ".join(row.split())
            result += "\\\\\n"
    result += "\\end{bmatrix}"
    print(result)


if __name__ == "__main__":
    A = """
1 2 0 0 3 2
0 0 1 0 −1 4
0 0 0 1 −2 3
0 0 0 0 0 0
    """
    texify_matrix(A)