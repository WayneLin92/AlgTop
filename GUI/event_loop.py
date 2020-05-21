import pygame
from GUI.app_ss import App, config


def draw_ss(spec_seq):
    """Interface for drawing spectral sequences."""
    pygame.init()
    surface = pygame.display.set_mode([config["win_width"], config["win_height"]])
    pygame.display.set_caption("Spectral Sequence")

    app = App(spec_seq, surface)
    app.render_grid_lines()

    loop_init = 0
    loop_exit = 1
    loop_index = loop_init
    while True: 
        if loop_index == loop_init:
            while True:
                event = pygame.event.wait()
                msg = event.dict
                if event.type == pygame.MOUSEBUTTONDOWN:
                    app.mouse_down(msg)
                elif event.type == pygame.MOUSEBUTTONUP:
                    app.mouse_up(msg)
                elif event.type == pygame.MOUSEMOTION:
                    app.mouse_motion(msg)
                elif event.type == pygame.KEYDOWN:
                    app.key_down(msg)
                elif event.type == pygame.KEYUP:
                    app.key_up(msg)
                if event.type == pygame.QUIT:
                    pygame.quit()
                    return None

                app.render()
                pygame.display.flip()
                pygame.time.wait(5)

        elif loop_index == loop_exit:
            pygame.quit()
            return None
