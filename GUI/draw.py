import pygame
from .pen_ss import Pen
from .constants import *


def draw_ss(spec_seq):  # TODO: fix the reopen failing problem.
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
            pygame.quit()
            return None
