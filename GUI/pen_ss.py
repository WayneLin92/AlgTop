import pygame
import pygame.gfxdraw
import math
from tkinter import Tk
from numpy import array
from .constants import *

BULLETS_SEP = 10
BULLETS_RADIUS = 3
BULLETS = None
BULLETS_DICTS = None

STATUS_START = 0
STATUS_ON_BULLET = 1
STATUS_ON_CANVAS = 2
STATUS_ON_CANVAS_MOVING = 3
ZOOM_RATE = 1.3

TYPE_IMAGE = 0
TYPE_KERNEL = 1
TYPE_DIFF = 2
TYPE_TBD = 3

font = None


def myround(r):
    return int(round(r))


def const_adjust(grid_width):
    global BULLETS_SEP, BULLETS_RADIUS, BULLETS, BULLETS_DICTS
    BULLETS_SEP = max(min(grid_width // 5, 25), 10)
    BULLETS_RADIUS = max(min(grid_width // 20, 5), 3)
    BULLETS = (
        ((0, 0),),
        ((-BULLETS_SEP // 2, 0), (BULLETS_SEP // 2, 0)),
        ((-BULLETS_SEP, 0), (0, 0), (BULLETS_SEP, 0)),
        ((-BULLETS_SEP // 2, -BULLETS_SEP // 2), (BULLETS_SEP // 2, -BULLETS_SEP // 2),
         (-BULLETS_SEP // 2, BULLETS_SEP // 2), (BULLETS_SEP // 2, BULLETS_SEP // 2)),
        ((0, -BULLETS_SEP), (-BULLETS_SEP, 0), (0, 0), (BULLETS_SEP, 0), (0, BULLETS_SEP)),
        ((-BULLETS_SEP, -BULLETS_SEP // 2), (0, -BULLETS_SEP // 2), (BULLETS_SEP, -BULLETS_SEP // 2),
         (-BULLETS_SEP, BULLETS_SEP // 2), (0, BULLETS_SEP // 2), (BULLETS_SEP, BULLETS_SEP // 2)),
        ((0, -BULLETS_SEP),
         (-BULLETS_SEP, 0), (0, 0), (BULLETS_SEP, 0),
         (-BULLETS_SEP, BULLETS_SEP), (0, BULLETS_SEP), (BULLETS_SEP, BULLETS_SEP)),
        ((-BULLETS_SEP, -BULLETS_SEP), (0, -BULLETS_SEP), (BULLETS_SEP, -BULLETS_SEP),
         (-BULLETS_SEP, 0), (BULLETS_SEP, 0),
         (-BULLETS_SEP, BULLETS_SEP), (0, BULLETS_SEP), (BULLETS_SEP, BULLETS_SEP)),
        ((-BULLETS_SEP, -BULLETS_SEP), (0, -BULLETS_SEP), (BULLETS_SEP, -BULLETS_SEP),
         (-BULLETS_SEP, 0), (0, 0), (BULLETS_SEP, 0),
         (-BULLETS_SEP, BULLETS_SEP), (0, BULLETS_SEP), (BULLETS_SEP, BULLETS_SEP)),
    )
    _bullet_diff = (tuple(((round((BULLETS[i][j][0] - BULLETS[i][0][0]) / BULLETS_SEP),
                            round((BULLETS[i][j][1] - BULLETS[i][0][1]) / BULLETS_SEP)), j)
                          for j in range(len(BULLETS[i]))) for i in range(len(BULLETS)))
    BULLETS_DICTS = [dict(d) for d in _bullet_diff]


def my_round(a):
    return array([int(round(x)) for x in a])


class VirtualScreen:
    def __init__(self):
        self.vir_origin = array([0, 0])
        self.vir_rate = 1.0

    def zoom(self, pos, rate):
        pos = array(pos)
        self.vir_rate *= rate
        self.vir_origin = pos + (self.vir_origin - pos) * rate
        const_adjust(GRID_WIDTH * self.vir_rate)

    def translation(self, rel):
        self.vir_origin += array(rel)

    def vir2surf(self, pos):
        screen = self.vir_origin + array(pos) * self.vir_rate
        return screen

    def surf2vir(self, pos):
        vir = (array(pos) - self.vir_origin) / self.vir_rate
        return vir


vs = VirtualScreen()


def draw_line(surface, color, start_pos, end_pos, width=1):
    pygame.draw.line(surface, color, my_round(start_pos), my_round(end_pos), width)


def draw_circle(surface, color, pos, radius, b_fill=True):
    if b_fill:
        pygame.gfxdraw.filled_circle(surface, *my_round(pos), int(round(radius)), color)
    else:
        pygame.gfxdraw.circle(surface, *my_round(pos), int(round(radius)), color)


def draw_text(surface, text, pos):
    text_img = font.render(text, True, (0, 0, 0), BG_COLOR)
    w, h = text_img.get_size()
    surface.blit(text_img, (myround(pos[0]) - w // 2, myround(pos[1]) - h // 2))


def draw_rect(surface, color, rect, width=0):
    pygame.draw.rect(surface, color, rect, width)


def c2array(z):
    return array((z.real, z.imag))


def array2c(a):
    return complex(a[0], a[1])


def draw_arrow(surface, start_pos, end_pos):
    """ Draw a classic arrow """
    d = array2c(end_pos) - array2c(start_pos)
    if abs(d) == 0:
        return
    d = d / abs(d)
    d1, d2 = d * (10+5j), d * (10-5j)
    end_pos = array(end_pos)
    draw_line(surface, PEN_COLOR, start_pos, end_pos, PEN_WIDTH)
    draw_line(surface, PEN_COLOR, end_pos - c2array(d1), end_pos, PEN_WIDTH)
    draw_line(surface, PEN_COLOR, end_pos - c2array(d2), end_pos, PEN_WIDTH)


def sp2deg(surf_pos):
    """ Return the closest grid point to surf_pos """
    vir_pos = vs.surf2vir(surf_pos)
    x = (vir_pos[0] - MARGIN_LEFT_RIGHT) / GRID_WIDTH
    y = (vir_pos[1] - MARGIN_TOP) / GRID_HEIGHT
    x = int(round(x))
    y = NUM_GRID_Y - int(round(y))
    return x, y


def deg2sp(deg):
    """ Return the surface coordinate of the deg """
    return vs.vir2surf([deg[0] * GRID_WIDTH + MARGIN_LEFT_RIGHT, (NUM_GRID_Y - deg[1]) * GRID_HEIGHT + MARGIN_TOP])


def offset_to_bullet(offset, n):
    """ Need 1 <= n <= 9 """
    diff = array(offset) - array(BULLETS[n-1][0])
    x = int(round(diff[0] / BULLETS_SEP))
    y = int(round(diff[1] / BULLETS_SEP))
    if (x, y) in BULLETS_DICTS[n-1]:
        return BULLETS_DICTS[n-1][(x, y)]
    else:
        return None


def copy_to_clipboard(content):
    """ Copy content to the clipboard """
    r = Tk()
    r.withdraw()
    r.clipboard_clear()
    r.clipboard_append(content)
    r.update()
    r.destroy()


class Pen:
    """
    The class for drawing
    self.expansion is a dictionary where the values are dictionaries
    """
    def __init__(self, ss):
        # for drawing
        self.status = STATUS_START
        self.addr_mouse_down = None
        self.alg_hover_on = None
        self.expansion = {}
        self.pos_current = None
        self.ctrl_down = False
        self.wait_for_update = True

        self.ss = ss
        self.bullets, self.arrows = ss.get_bullets_and_arrows()

        const_adjust(GRID_WIDTH)
        global font
        font = pygame.font.SysFont("Arial", 20)

    def exp_collide_point(self, pos):
        for deg, v in self.expansion.items():
            rect = v['rect']
            if rect.collidepoint(pos):
                return deg
        return None

    def expand(self, deg):
        """ rect is a surface rect """
        n = len(self.bullets[deg])
        num_edge = math.ceil(math.sqrt(n))
        length_edge = BULLETS_SEP * (num_edge + 1)
        size = array([length_edge, length_edge])
        center = deg2sp(deg)
        rect = pygame.Rect((center - size // 2), size)
        self.expansion[deg] = {'n': n, 'num_edge': num_edge, 'rect': rect}

    def exp_regen(self):
        for deg in self.expansion:
            self.expand(deg)

    def exp_close(self, deg):
        del self.expansion[deg]

    def exp_close_all(self):
        self.expansion = {}

    def sp2addr(self, surf_pos):
        """ Return (deg, index) where deg is for the expansion """
        deg = self.exp_collide_point(surf_pos)
        if deg is not None:
            rect = self.expansion[deg]['rect']
            n = self.expansion[deg]['n']
            num_edge = self.expansion[deg]['num_edge']
            diff = array(surf_pos) - rect.topleft
            x = int(round(diff[0] / BULLETS_SEP))
            y = int(round(diff[1] / BULLETS_SEP))
            if 0 < x <= num_edge and 0 < y <= num_edge:
                index = x - 1 + num_edge * (y - 1)
                if index < n:
                    return deg, index
            return None
        deg = sp2deg(surf_pos)
        if deg in self.bullets:
            n = len(self.bullets[deg])
            if 0 < n <= 9:
                offset = array(surf_pos) - array(deg2sp(deg))
                index = offset_to_bullet(offset, n)
                if index is not None:
                    return deg, index
                else:
                    return None
            elif n > 9:
                offset = array(surf_pos) - array(deg2sp(deg))
                index = offset_to_bullet(offset, 9)
                if index is not None:
                    return deg, None
                else:
                    return None
        else:
            return None

    def addr2sp(self, addr):
        deg, index = addr
        if index is None:
            return deg2sp(deg)
        n = len(self.bullets[deg])
        if n <= 9:
            return deg2sp(deg) + array(BULLETS[n - 1][index])
        else:
            if deg in self.expansion:
                rect = self.expansion[deg]['rect']
                num_edge = self.expansion[deg]['num_edge']
                y, x = divmod(index, num_edge)
                pos_bullet = rect.topleft + array([x + 1, y + 1]) * BULLETS_SEP
                return pos_bullet
            else:
                return deg2sp(deg)
    
    def addr2alg(self, addr):
        deg, index = addr
        return self.bullets[deg][index][0]
    
    def addr2color(self, addr):
        deg, index = addr
        type_bullet = self.bullets[deg][index][1]
        return COLOR_BULLET_TBD if type_bullet == TYPE_TBD else COLOR_BULLET_NON_TBD

    def alg2addr(self, alg):
        deg = alg.deg()
        for i, bullet in enumerate(self.bullets[deg]):
            if alg == bullet[0]:
                return deg, i

    def mouse_down(self, msg):
        # todo: context manager
        # todo: add generators
        if msg['button'] == 1:
            if self.status == STATUS_START:
                self.addr_mouse_down = self.sp2addr(msg['pos'])
                if self.addr_mouse_down is None:
                    self.status = STATUS_ON_CANVAS
                elif self.addr_mouse_down[1] is None:
                    if not self.ctrl_down:
                        self.exp_close_all()
                    self.expand(self.addr_mouse_down[0])
                    self.wait_for_update = True
                else:
                    deg, index = self.addr_mouse_down
                    if self.bullets[deg][index][1] == TYPE_TBD:
                        self.status = STATUS_ON_BULLET
                    else:
                        self.addr_mouse_down = None
                        self.status = STATUS_ON_CANVAS
        elif msg['button'] == 4:
            vs.zoom(msg['pos'], ZOOM_RATE)
            self.exp_regen()
            self.wait_for_update = True
        elif msg['button'] == 5:
            vs.zoom(msg['pos'], 1 / ZOOM_RATE)
            self.exp_regen()
            self.wait_for_update = True

    def mouse_up(self, msg):
        if self.status == STATUS_ON_BULLET:  # add the new arrow
            addr1 = self.addr_mouse_down
            addr2 = self.sp2addr(msg['pos'])
            if addr2 is not None and addr2[1] is not None:
                if addr1[0] != addr2[0]:
                    poly1, poly2 = self.addr2alg(addr1), self.addr2alg(addr2) 
                    self.ss.add_diff(poly1, poly2)
                    self.bullets, self.arrows = self.ss.get_bullets_and_arrows()
                    self.alg_hover_on = None
                self.status = STATUS_START
            elif addr2 is None or addr2[0] == addr1[0]:
                self.status = STATUS_START
            else:
                self.expand(addr2[0])
            self.wait_for_update = True
        elif self.status == STATUS_ON_CANVAS_MOVING:
            self.status = STATUS_START
        elif self.status == STATUS_ON_CANVAS:
            self.exp_close_all()
            self.status = STATUS_START
            self.wait_for_update = True

    def mouse_motion(self, msg):
        self.pos_current = msg['pos']
        if self.status == STATUS_ON_BULLET:
            deg = self.addr_mouse_down[0]
            if deg in self.expansion:
                rect = self.expansion[deg]['rect']
                if not rect.collidepoint(msg['pos']):
                    self.exp_close(deg)
            self.wait_for_update = True
        if self.status in {STATUS_ON_CANVAS, STATUS_ON_CANVAS_MOVING}:
            vs.translation(msg['rel'])
            self.exp_regen()
            self.status = STATUS_ON_CANVAS_MOVING
            self.wait_for_update = True
        addr = self.sp2addr(msg['pos'])
        if addr is not None and addr[1] is not None:
            if self.alg_hover_on is None or self.alg_hover_on != self.bullets[addr[0]][addr[1]][0]:
                self.alg_hover_on = self.bullets[addr[0]][addr[1]][0]
                self.wait_for_update = True
        else:
            if self.alg_hover_on is not None:
                self.alg_hover_on = None
                self.wait_for_update = True

    def key_down(self, msg):
        if msg['unicode'] == '\x03':  # Ctrl+C => output for tex
            # todo: save and open diagrams
            x_min = NUM_GRID_X
            y_min = NUM_GRID_Y
            x_max = 0
            y_max = 0
            entries = self.entries
            for pos in entries:
                if x_min > pos[0]:
                    x_min = pos[0]
                if y_min > pos[1]:
                    y_min = pos[1]
                if x_max < pos[0]:
                    x_max = pos[0]
                if y_max < pos[1]:
                    y_max = pos[1]
            w = x_max - x_min
            h = y_max - y_min
            tex_arrows = dict()
            for arrow in self.arrows:
                i = arrow[0][0] - x_min
                j = arrow[0][1] - y_min
                di = arrow[1][0] - arrow[0][0]
                dj = arrow[1][1] - arrow[0][1]
                if (i, j) in tex_arrows.keys():
                    tex_arrows[(i, j)].append((di, dj))
                else:
                    tex_arrows[(i, j)] = [(di, dj)]
            tex_string = "\\xymatrix{\n"
            for j in range(h + 1):
                for i in range(w + 1):
                    if (x_min + i, y_min + j) in entries:
                        tex_string += str(entries[(x_min + i, y_min + j)])
                    if (i, j) in tex_arrows.keys():
                        for arrow in tex_arrows[(i, j)]:
                            tex_string += " \\ar"
                            tex_string += "["
                            tex_string += "r" * arrow[0] if arrow[0] > 0 else "l" * (-arrow[0])
                            tex_string += "d" * arrow[1] if arrow[1] > 0 else "u" * (-arrow[1])
                            tex_string += "]"
                    if i < w:
                        tex_string += " & "
                tex_string += "\\\\\n"
            tex_string += "}"
            print(tex_string)
            copy_to_clipboard(tex_string)
        elif msg['key'] == 127:  # Del -> initialize every thing
            pass
        elif msg['key'] == 8:  # Backspace -> delete an arrow
            pass
        elif msg['key'] == 306:  # Ctrl
            self.ctrl_down = True
        elif msg['unicode'] == 'r':  # R -> initialize the view setting
            vs.__init__()
            self.exp_regen()
        elif msg['unicode'] == '\x0e':  # Ctrl + N -> Compute and move to the next page
            self.ss.new_page()
            self.exp_close_all()
            self.alg_hover_on = None
            self.addr_mouse_down = None
            self.wait_for_update = True

    def key_up(self, msg):
        # vir_pos = vs.surf2vir(msg['pos'])
        if msg['key'] == 306:  # Ctrl
            self.ctrl_down = False

    @staticmethod
    def bg(surface):
        """ Draw the grid """
        surface.fill(BG_COLOR)
        for i in range(0, NUM_GRID_Y + 1):
            left = vs.vir2surf((MARGIN_LEFT_RIGHT, MARGIN_TOP + GRID_HEIGHT * i))
            right = vs.vir2surf((WIN_WIDTH - MARGIN_LEFT_RIGHT, MARGIN_TOP + GRID_HEIGHT * i))
            draw_line(surface, GRID_COLOR, left, right)
            text_pos = left - array([30, 0])
            if text_pos[0] < 16:
                text_pos[0] = 16
            draw_text(surface, str(NUM_GRID_Y - i), text_pos)
        for i in range(0, NUM_GRID_X + 1):
            top = vs.vir2surf((MARGIN_LEFT_RIGHT + GRID_WIDTH * i, MARGIN_TOP))
            bottom = vs.vir2surf((MARGIN_LEFT_RIGHT + GRID_WIDTH * i, WIN_HEIGHT - MARGIN_BOTTOM))
            draw_line(surface, GRID_COLOR, top, bottom)
            text_pos = bottom + array([0, 30])
            if text_pos[1] > WIN_HEIGHT - 16:
                text_pos[1] = WIN_HEIGHT - 16
            draw_text(surface, str(i), text_pos)

    def render(self, surface):
        """ Draw on screen """
        # TODO dashed line for the region of next page
        if not self.wait_for_update:
            return
        else:
            self.wait_for_update = False

        self.bg(surface)
        # draw bullets
        for deg in self.bullets:
            n = len(self.bullets[deg])
            if n <= 9:
                for i in range(n):
                    draw_circle(surface, self.addr2color((deg, i)),
                                array(deg2sp(deg)) + array(BULLETS[n-1][i]), BULLETS_RADIUS)
            else:
                draw_circle(surface, PEN_COLOR, deg2sp(deg), BULLETS_RADIUS)
                draw_circle(surface, PEN_COLOR, deg2sp(deg), BULLETS_RADIUS * 2, False)
        # draw expansions
        for deg, v in self.expansion.items():
            n = v['n']
            rect = v['rect']
            num_edge = v['num_edge']
            draw_rect(surface, WHITE, rect)
            draw_rect(surface, PEN_COLOR, rect, 1)
            for i in range(n):
                y, x = divmod(i, num_edge)
                pos_bullet = array(rect.topleft) + array((x + 1, y + 1)) * BULLETS_SEP
                draw_circle(surface, self.addr2color((deg, i)), pos_bullet, BULLETS_RADIUS)
        # enlarge hovered-on bullet
        if self.alg_hover_on is not None:
            deg, i = self.alg2addr(self.alg_hover_on)
            n = len(self.bullets[deg])
            if n <= 9:
                draw_circle(surface, self.addr2color((deg, i)),
                            array(deg2sp(deg)) + array(BULLETS[n - 1][i]), BULLETS_RADIUS * 1.5)
            else:
                v = self.expansion[deg]
                rect = v['rect']
                num_edge = v['num_edge']
                y, x = divmod(i, num_edge)
                pos_bullet = array(rect.topleft) + array((x + 1, y + 1)) * BULLETS_SEP
                draw_circle(surface, self.addr2color((deg, i)), pos_bullet, BULLETS_RADIUS * 1.5)
        # draw arrows
        # for arrow in self.ss.get_arrows(self.expansion):
        for arrow in self.arrows:
            draw_arrow(surface, self.addr2sp(arrow[0]), self.addr2sp(arrow[1]))
        if self.status == STATUS_ON_BULLET:
            if self.addr_mouse_down is not None and self.addr_mouse_down[1] is not None:
                draw_arrow(surface, self.addr2sp(self.addr_mouse_down), self.pos_current)
        # draw label
        if self.alg_hover_on is not None:
            s = str(self.alg_hover_on)
            draw_text(surface, s, (WIN_WIDTH // 2, WIN_HEIGHT-MARGIN_BOTTOM // 2))

# 477
