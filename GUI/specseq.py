"""Class of spectral sequences.

This class mainly handles the logical data of the spectral sequence."""
from typing import List, Tuple, NamedTuple


class BulletProperties(NamedTuple):
    id: int
    deg: Tuple[int, int]
    label: str
    color: tuple


class LineProperties(NamedTuple):
    id: int
    src_id: int
    tgt_id: int
    label: str
    color: tuple


class ArrowProperties(NamedTuple):
    id: int
    src_id: int
    tgt_id: int


class TextProperties(NamedTuple):
    id: int
    text: str
    deg: Tuple[int, int]
    offset: Tuple[float, float]
    color: tuple


class SpecSeq:
    """This class is the database for the graph of spectral sequence."""
    def __init__(self):
        self.bullets = []  # type: List[BulletProperties]
        self.lines = []  # type: List[LineProperties]
        self.arrows = []  # type: List[ArrowProperties]
        self.texts = []  # type: List[TextProperties]
        self.bullet_id = 0
        self.line_id = 0
        self.arrow_id = 0
        self.text_id = 0

    # getters -------------------------
    def get_bullet_by_id(self, id_):
        for b in self.bullets:
            if b.id == id_:
                return b
        else:
            raise ValueError(f"Id {id_} not found.")

    def get_bullet_by_label(self, label):
        """Warning: labels might be non-unique."""
        for b in self.bullets:
            if b.label == label:
                return b
        else:
            raise ValueError(f"label {label} not found.")

    # setters --------------------------
    def add_bullet(self, deg, label=None, color=None):
        bullet = BulletProperties(self.bullet_id, deg, label or f"id: {self.bullet_id}", color)
        self.bullets.append(bullet)
        self.bullet_id += 1
        return self.bullet_id - 1

    def add_line(self, src, tgt, label=None, color=None, *, by_label=False):
        if by_label:  # Warning: labels might be non-unique
            src = self.get_bullet_by_label(src).id
            tgt = self.get_bullet_by_label(tgt).id
        line = LineProperties(self.line_id, src, tgt, label or f"id: {self.line_id}", color)
        self.lines.append(line)
        self.line_id += 1

    def add_arrow(self, src, tgt, *, by_label=False):
        if by_label:  # Warning: labels might be non-unique
            src = self.get_bullet_by_label(src).id
            tgt = self.get_bullet_by_label(tgt).id
        arrow = ArrowProperties(self.arrow_id, src, tgt)
        self.arrows.append(arrow)
        self.arrow_id += 1

    def add_remove_arrow(self, src, tgt, *, by_label=False):
        if by_label:  # Warning: labels might be non-unique
            src = self.get_bullet_by_label(src).id
            tgt = self.get_bullet_by_label(tgt).id
        for i, arrow in enumerate(self.arrows):
            if arrow.src_id == src and arrow.tgt_id == tgt:
                break
        else:
            arrow = ArrowProperties(self.arrow_id, src, tgt)
            self.arrows.append(arrow)
            self.arrow_id += 1
            return
        del self.arrows[i]

    def add_text(self, text: str, deg: Tuple[int, int], offset: Tuple[float, float] = (0.0, 0.0), color=None):
        text_p = TextProperties(self.text_id, text, deg, offset, color)
        self.texts.append(text_p)
        self.text_id += 1

    # methods ----------------------------
    def draw(self, caption=None):
        import GUI.event_loop
        GUI.event_loop.draw_ss(self, caption)


def test():
    spec = SpecSeq()
    b1 = spec.add_bullet((1, 1), color=(255, 255, 0))
    b2 = spec.add_bullet((2, 1))
    b3 = spec.add_bullet((1, 2))
    spec.add_arrow(b1, b2)
    spec.add_line(b1, b3, color=(0, 255, 0))
    spec.add_text("12341123412341234", (4, 4), color=(0, 255, 255))

    for i in range(12):
        for j in range(i):
            spec.add_bullet((i, 2), color=(0, 0, 255))
    spec.draw()


if __name__ == "__main__":
    test()

# 623, 193, 184
