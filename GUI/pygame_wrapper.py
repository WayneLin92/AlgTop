"""pygame wrapper"""
import os
import json
import pygame
import pygame.gfxdraw
from algebras.mymath import Vector


with open(f"{os.path.dirname(__file__)}\\pen_ss.json", "r") as file:
    config = json.load(file)
    config["bullet.patterns.offset"] = [{tuple(pos): i for i, pos in enumerate(pattern)}
                                        for pattern in config["bullet.patterns"]]
    config["bullet.patterns"] = [[Vector(pos) for pos in pattern] for pattern in config["bullet.patterns"]]


def myroundv(t) -> Vector:
    return Vector(map(round, t))


def c2Vector(z: complex):
    return Vector((z.real, z.imag))


def Vector2c(a):
    return complex(a[0], a[1])


class Paint:
    def __init__(self, surface):
        self.surface = surface
        self.font = pygame.font.SysFont("Arial", 20)

    def clear_screen(self):
        self.surface.fill(config["bg_color"])

    def draw_line(self, color, start_pos, end_pos, width=config["pen_width"]):
        pygame.draw.line(self.surface, color or config["pen_color"], myroundv(start_pos), myroundv(end_pos), width)

    def draw_circle(self, color, pos, radius, b_fill=True):
        if b_fill:
            pygame.gfxdraw.filled_circle(self.surface, *myroundv(pos), int(radius), color or config["pen_color"])
        else:
            pygame.gfxdraw.circle(self.surface, *myroundv(pos), int(radius), color or config["pen_color"])

    def draw_rect(self, color, rect, width=0):
        pygame.draw.rect(self.surface, color or config["pen_color"], rect, width)

    def draw_text(self, text, pos, color=None):
        """The center of the text is placed at pos."""
        text_img = self.font.render(text, True, color or config["pen_color"], config["bg_color"])
        w, h = text_img.get_size()
        self.surface.blit(text_img, (round(pos[0]) - w // 2, round(pos[1]) - h // 2))

    def draw_arrow(self, start_pos, end_pos):
        """ Draw a classic arrow """
        d = Vector2c(end_pos) - Vector2c(start_pos)
        if abs(d) == 0:
            return
        d = d / abs(d)
        d1, d2 = d * (10+5j), d * (10-5j)
        end_pos = Vector(end_pos)
        self.draw_line(config["pen_color"], start_pos, end_pos)
        self.draw_line(config["pen_color"], end_pos - c2Vector(d1), end_pos)
        self.draw_line(config["pen_color"], end_pos - c2Vector(d2), end_pos)
