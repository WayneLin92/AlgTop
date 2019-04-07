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


def texify_aug_matrix(text: str):
    result = "\\augmatrix{"
    rows = text.split(sep="\n")
    first_row = True
    for row in rows:
        r = row.split()
        if r and first_row:
            result += "c" * (len(r) - 1) + "|c}{\n"
            first_row = False
        if r:
            result += " & ".join(row.split())
            result += "\\\\\n"
    result += "}"
    print(result)


if __name__ == "__main__":
    A = """
1 0 0 \\cdots 0
0 1 0 \\cdots 0
0 0 1 \\cdots 0
\\cdots \\cdots \\cdots \\cdots \\cdots
0 0 0 \\cdots 1
    """
    texify_matrix(A)
