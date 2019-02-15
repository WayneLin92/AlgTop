import pygame
from .pen_ss import Pen
from .constants import *


def draw_ss(spec_seq):
    """ interface for drawing spectral sequences """
    pygame.init()
    surface = pygame.display.set_mode([WIN_WIDTH, WIN_HEIGHT])
    pygame.display.set_caption("Spectral Sequence")

    pen = Pen(spec_seq)
    pen.bg(surface)

    loop_init = 0
    loop_exit = 1
    loop_index = loop_init
    while True: 
        if loop_index == loop_init:
            while True:
                event = pygame.event.wait()
                msg = event.dict
                if event.type == pygame.MOUSEBUTTONDOWN:
                    pen.mouse_down(msg)
                elif event.type == pygame.MOUSEBUTTONUP:
                    pen.mouse_up(msg)
                elif event.type == pygame.MOUSEMOTION:
                    pen.mouse_motion(msg)
                elif event.type == pygame.KEYDOWN:
                    pen.key_down(msg)
                elif event.type == pygame.KEYUP:
                    pen.key_up(msg)
                elif event.type == pygame.QUIT:
                    pygame.quit()
                    return None

                pen.render(surface)
                pygame.display.flip()
                pygame.time.wait(5)

        elif loop_index == loop_exit:
            pygame.display.quit()
            pygame.quit()
            return None


def draw_01matrix(m):
    """ draw a 0-1 (black-white) square matrix """
    pygame.init()
    surface = pygame.display.set_mode([640, 640])
    surface.fill((0, 0, 0))
    d = len(m)
    grid_width = 640 // d
    for i in range(d):
        for j in range(d):
            if m[i][j]:
                rect = pygame.Rect(grid_width * j, grid_width * i, grid_width, grid_width)
                pygame.draw.rect(surface, (255, 255, 255), rect)
    pygame.display.flip()
    while True:
        event = pygame.event.wait()
        if event.type == pygame.QUIT:
            pygame.quit()
            break


def draw_01triangle(m):
    """ draw an 0-1 (red-black) triangle """
    pygame.init()
    surface = pygame.display.set_mode([640, 640])
    surface.fill((250, 250, 250))
    d = len(m)
    top, w = 300j, -0.5 + 3 ** 0.5 * 0.5j
    w1 = 0.5 + 3 ** 0.5 * 0.5j
    origin = 320+320j + (top * w).conjugate()
    grid_width = 900 / 3 ** 0.5 / d
    v1 = grid_width + 0j
    v2 = v1 / w1
    for i in range(d + 1):
        for j in range(d + 1):
            if i + j < d:
                pt = origin + (v1 * i + v2 * j)
                x, y = int(round(pt.real)), int(round(pt.imag))
                if m[i][j]:
                    pygame.draw.circle(surface, (0, 0, 0), (x, y), int(grid_width / 2))
                else:
                    pygame.draw.circle(surface, (255, 128, 128), (x, y), int(grid_width / 2))
    pygame.display.flip()
    while True:
        event = pygame.event.wait()
        if event.type == pygame.QUIT:
            pygame.quit()
            break
