import os, json
import pygame
import pygame.gfxdraw
from algebras.mymath import Vector


def myroundv(t) -> Vector:
    return Vector(map(round, t))


with open(f"{os.path.dirname(__file__)}\\pen_ss.json", "r") as file:
    config = json.load(file)
    config["bullet.patterns.offset"] = [{tuple(pos): i for i, pos in enumerate(pattern)}
                                        for pattern in config["bullet.patterns"]]
    config["bullet.patterns"] = [[Vector(pos) for pos in pattern] for pattern in config["bullet.patterns"]]


def draw_line(surface, color, start_pos, end_pos, width=config["pen_width"]):
    pygame.draw.line(surface, color, myroundv(start_pos), myroundv(end_pos), width)


def draw_circle(surface, color, pos, radius, b_fill=True):
    if b_fill:
        pygame.gfxdraw.filled_circle(surface, *myroundv(pos), int(radius), color)
    else:
        pygame.gfxdraw.circle(surface, *myroundv(pos), int(radius), color)


def draw_text(surface, text, pos, font):
    text_img = font.render(text, True, config["pen_color"], config["bg_color"])
    w, h = text_img.get_size()
    surface.blit(text_img, (round(pos[0]) - w // 2, round(pos[1]) - h // 2))


def draw_rect(surface, color, rect, width=0):
    pygame.draw.rect(surface, color, rect, width)


def c2Vector(z: complex):
    return Vector((z.real, z.imag))


def Vector2c(a):
    return complex(a[0], a[1])


def draw_arrow(surface, start_pos, end_pos):
    """ Draw a classic arrow """
    d = Vector2c(end_pos) - Vector2c(start_pos)
    if abs(d) == 0:
        return
    d = d / abs(d)
    d1, d2 = d * (10+5j), d * (10-5j)
    end_pos = Vector(end_pos)
    draw_line(surface, config["pen_color"], start_pos, end_pos)
    draw_line(surface, config["pen_color"], end_pos - c2Vector(d1), end_pos)
    draw_line(surface, config["pen_color"], end_pos - c2Vector(d2), end_pos)
