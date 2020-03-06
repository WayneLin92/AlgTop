"""Class of spectral sequences."""
from itertools import groupby
from typing import List, Tuple, NamedTuple

TYPE_IMAGE = 0
TYPE_KERNEL = 1
TYPE_DIFF = 2
TYPE_TBD = 3


class BulletProperties(NamedTuple):
    id: int
    deg: Tuple[int, int]
    label: str
    color: str


class LineProperties(NamedTuple):
    id: int
    src: int
    tgt: int
    label: str
    color: str


class ArrowProperties(NamedTuple):
    id: int
    src: int
    tgt: int


class SpecSeq:
    """This class is the database for the graph of spectral sequence."""
    def __init__(self):
        self.bullets = []  # type: List[BulletProperties]
        self.lines = []  # type: List[LineProperties]
        self.arrows = []  # type: List[ArrowProperties]
        self.bullet_id = 0
        self.line_id = 0
        self.arrow_id = 0

    # getters -------------------------
    def get_bullet_by_id(self, id_):
        for b in self.bullets:
            if b.id == id_:
                return b

    def get_bullets_and_arrows(self):
        bullets = {}
        self.bullets.sort(key=lambda p: p.deg)
        id2index = {}
        for deg, group in groupby(self.bullets, key=lambda p: p.deg):
            bullets[deg] = []
            for i, b in enumerate(group):
                bullets[deg].append((b.id, TYPE_TBD))
                id2index[b.id] = i
        arrows = []
        for a in self.arrows:
            src = self.get_bullet_by_id(a.src)
            tgt = self.get_bullet_by_id(a.tgt)
            arrows.append(((src.deg, id2index[src.id]), (tgt.deg, id2index[tgt.id])))
        return bullets, arrows

    # setters --------------------------
    def add_bullet(self, deg, label="", color="#000000"):
        bullet = BulletProperties(self.bullet_id, deg, label, color)
        self.bullets.append(bullet)
        self.bullet_id += 1

    def add_arrow(self, src, tgt):
        arrow = ArrowProperties(self.arrow_id, src, tgt)
        self.arrows.append(arrow)
        self.arrow_id += 1

    # methods ----------------------------
    def add_diff(self, id1, id2):
        self.add_arrow(id1, id2)

    def draw(self):
        import GUI.draw
        GUI.draw.draw_ss(self)


def test():
    spec = SpecSeq()
    spec.add_bullet((1, 1))
    spec.add_bullet((2, 1))
    spec.add_bullet((1, 2))
    spec.draw()
    return spec


if __name__ == "__main__":
    test()

# 623, 193, 184
