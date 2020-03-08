import math
from itertools import groupby
from tkinter import Tk
from algebras.mymath import Vector
from GUI.draw import pygame, config, draw_line, draw_text, draw_arrow, draw_rect, draw_circle, myround
from GUI.specseq import SpecSeq


# TODO: change status to string
# TODO: auto expand
# TODO: improve grid lines and labels

STATUS_START = 0
STATUS_ON_BULLET = 1
STATUS_ON_CANVAS = 2
STATUS_ON_CANVAS_MOVING = 3
ZOOM_RATE = 1.3


def modify_config(grid_width):
    config["bullet.sep"] = max(min(grid_width // 5, 25), 10)
    config["bullet.radius"] = max(min(grid_width // 20, 5), 3)


class Camera:
    def __init__(self):
        """`origin_sp` is the screen position of world origin.
        `unit_screen` is the screen length of world unit length."""
        self.origin_sp = Vector((config["margin_left"], config["win_height"] - config["margin_bottom"]))
        self.unit_screen = config["grid_width"]
        
    @staticmethod
    def clip(x, interval):
        """Clip value `x` by `interval`."""
        if x < interval[0]:
            return interval[0]
        elif x > interval[1]:
            return interval[1]
        else:
            return x

    @staticmethod
    def flip(world_pos):
        return Vector((world_pos[0], -world_pos[1]))

    # setters -------------------
    def zoom(self, pos, rate):
        pos = Vector(pos)
        unit_screen1 = self.clip(self.unit_screen * rate, config["camera.unit_screen.interval"])
        rate1 = unit_screen1 / self.unit_screen
        self.unit_screen = unit_screen1
        origin_sp1 = pos + (self.origin_sp - pos) * rate1
        x1 = self.clip(origin_sp1[0], config["camera.origin_sp.x.interval"])
        y1 = self.clip(origin_sp1[1], config["camera.origin_sp.y.interval"])
        self.origin_sp = Vector((x1, y1))
        modify_config(self.unit_screen)

    def translation(self, rel):
        origin_sp1 = self.origin_sp + Vector(rel)
        x1 = self.clip(origin_sp1[0], config["camera.origin_sp.x.interval"])
        y1 = self.clip(origin_sp1[1], config["camera.origin_sp.y.interval"])
        self.origin_sp = Vector((x1, y1))

    # methods -----------------
    def wp2sp(self, world_pos):
        """Convert world position to screen position."""
        return self.origin_sp + self.flip(world_pos) * self.unit_screen

    def sp2wp(self, pos):
        """Convert screen position to world position."""
        return self.flip((Vector(pos) - self.origin_sp) / self.unit_screen)


camera = Camera()


def sp2deg(screen_pos):
    """Return the closest grid point to screen_pos."""
    wp = camera.sp2wp(screen_pos)
    return myround(wp)


def deg2sp(deg):
    """Return the screen position of the deg."""
    return camera.wp2sp(deg)


def offset_to_bullet(offset: Vector, n):
    """Need 1 <= n <= 9."""
    v = myround(offset / config["bullet.sep"] * 2)
    if v in config["bullet.patterns.offset"][n - 1]:
        return config["bullet.patterns.offset"][n - 1][v]
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
    def __init__(self, ss: SpecSeq):
        # for drawing
        self.status = STATUS_START
        self.addr_mouse_down = None
        self.id_hover_on = None
        self.expansion = {}
        self.pos_current = None
        self.ctrl_down = False
        self.wait_for_update = True

        self.ss = ss
        self.id2addr = {}
        self.addr2id = {}
        self.deg2num_bullets = {}
        self.init_from_ss(ss)

        modify_config(config["grid_width"])
        self.font = pygame.font.SysFont("Arial", 20)

    def init_from_ss(self, ss):
        ss.bullets.sort(key=lambda p: p.deg)
        for deg, group in groupby(ss.bullets, key=lambda p: p.deg):
            i = -1
            for i, b in enumerate(group):
                self.id2addr[b.id] = (deg, i)
                self.addr2id[(deg, i)] = b.id
            self.deg2num_bullets[deg] = i + 1

    def exp_collide_point(self, pos):
        for deg, v in self.expansion.items():
            rect = v['rect']
            if rect.collidepoint(pos):
                return deg
        return None

    def expand(self, deg):
        """ rect is a surface rect """
        n = self.deg2num_bullets[deg]
        num_edge = math.ceil(math.sqrt(n))
        length_edge = config["bullet.sep"] * (num_edge + 1)
        size = Vector([length_edge, length_edge])
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

    def sp2addr(self, screen_pos):
        """Return (deg, index) where deg is for the expansion."""
        deg = self.exp_collide_point(screen_pos)
        if deg is not None:
            rect = self.expansion[deg]['rect']
            n = self.expansion[deg]['n']
            num_edge = self.expansion[deg]['num_edge']
            diff = Vector(screen_pos) - rect.topleft
            x, y = myround(diff / config["bullet.sep"])
            if 0 < x <= num_edge and 0 < y <= num_edge:
                index = x - 1 + num_edge * (y - 1)
                if index < n:
                    return deg, index
            return None
        deg = sp2deg(screen_pos)
        if deg in self.deg2num_bullets:
            n = self.deg2num_bullets[deg]
            if 0 < n <= 9:
                offset = Vector(screen_pos) - Vector(deg2sp(deg))
                index = offset_to_bullet(offset, n)
                if index is not None:
                    return deg, index
                else:
                    return None
            elif n > 9:
                offset = Vector(screen_pos) - Vector(deg2sp(deg))
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
        n = self.deg2num_bullets[deg]
        if n <= 9:
            offset = config["bullet.patterns"][n - 1][index] * config["bullet.sep"] // 2
            return deg2sp(deg) + offset
        else:
            if deg in self.expansion:
                rect = self.expansion[deg]['rect']
                num_edge = self.expansion[deg]['num_edge']
                y, x = divmod(index, num_edge)
                pos_bullet = rect.topleft + Vector([x + 1, y + 1]) * config["bullet.sep"]
                return pos_bullet
            else:
                return deg2sp(deg)

    def addr2color(self, addr):  #
        return pygame.Color(self.ss.get_bullet_by_id(self.addr2id[addr]).color)

    def mouse_down(self, msg):
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
                    self.status = STATUS_ON_BULLET
        elif msg['button'] == 4:
            camera.zoom(msg['pos'], ZOOM_RATE)
            self.exp_regen()
            self.wait_for_update = True
        elif msg['button'] == 5:
            camera.zoom(msg['pos'], 1 / ZOOM_RATE)
            self.exp_regen()
            self.wait_for_update = True

    def mouse_up(self, msg):
        if self.status == STATUS_ON_BULLET:  # add the new arrow
            addr1 = self.addr_mouse_down
            addr2 = self.sp2addr(msg['pos'])
            if addr2 is not None and addr2[1] is not None:
                if addr1[0] != addr2[0]:
                    id1, id2 = self.addr2id[addr1], self.addr2id[addr2]
                    self.ss.add_arrow(id1, id2)
                    # print(len(self.ss.arrows))
                    self.id_hover_on = None
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
            camera.translation(msg['rel'])
            self.exp_regen()
            self.status = STATUS_ON_CANVAS_MOVING
            self.wait_for_update = True
        addr = self.sp2addr(msg['pos'])
        if addr is not None and addr[1] is not None:
            if self.id_hover_on is None or self.id_hover_on != self.addr2id[addr]:
                self.id_hover_on = self.addr2id[addr]
                self.wait_for_update = True
        else:
            if self.id_hover_on is not None:
                self.id_hover_on = None
                self.wait_for_update = True

    def key_down(self, msg):
        if msg['unicode'] == '\x03':  # Ctrl+C => output for tex
            # todo: save and open diagrams
            x_min = config["num_grid_x"]
            y_min = config["num_grid_y"]
            x_max = 0
            y_max = 0
            entries = self.deg2num_bullets
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
            for arrow in self.ss.arrows:
                deg1 = self.ss.get_bullet_by_id(arrow.src).deg
                deg2 = self.ss.get_bullet_by_id(arrow.tgt).deg
                i = deg1[0] - x_min
                j = deg1[1] - y_min
                di = deg2[0] - deg1[0]
                dj = deg2[1] - deg1[1]
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
            camera.__init__()
            self.exp_regen()

    def key_up(self, msg):
        # vir_pos = vs.sp2wp(msg['pos'])
        if msg['key'] == 306:  # Ctrl
            self.ctrl_down = False

    def bg(self, surface):
        """ Draw the grid """
        surface.fill(config["bg_color"])
        for i in range(0, config["num_grid_y"] + 1):
            left = camera.wp2sp((0, i))
            right = camera.wp2sp((config["num_grid_x"], i))
            draw_line(surface, config["grid_color"], left, right)
            text_pos = left - Vector([30, 0])
            if text_pos[0] < 16:
                text_pos = Vector((16, text_pos[1]))
            draw_text(surface, str(i), text_pos, self.font)
        for i in range(0, config["num_grid_x"] + 1):
            bottom = camera.wp2sp((i, 0))
            top = camera.wp2sp((i, config["num_grid_y"]))
            draw_line(surface, config["grid_color"], top, bottom)
            text_pos = bottom + Vector([0, 30])
            if text_pos[1] > config["win_height"] - 16:
                text_pos = Vector((text_pos[0], config["win_height"] - 16))
            draw_text(surface, str(i), text_pos, self.font)

    def render(self, surface):
        """Draw on the screen."""
        if not self.wait_for_update:
            return
        else:
            self.wait_for_update = False

        self.bg(surface)
        # draw bullets
        for deg in self.deg2num_bullets:
            n = self.deg2num_bullets[deg]
            if n <= 9:
                for i in range(n):
                    offset = config["bullet.patterns"][n - 1][i] * config["bullet.sep"] / 2
                    draw_circle(surface, self.addr2color((deg, i)),
                                Vector(deg2sp(deg)) + offset, config["bullet.radius"])
            else:
                draw_circle(surface, config["pen_color"], deg2sp(deg), config["bullet.radius"])
                draw_circle(surface, config["pen_color"], deg2sp(deg), config["bullet.radius"] * 2, False)
        # draw expansions
        for deg, v in self.expansion.items():
            n = v['n']
            rect = v['rect']
            num_edge = v['num_edge']
            draw_rect(surface, pygame.Color("#dddddd"), rect)
            draw_rect(surface, config["pen_color"], rect, 1)
            for i in range(n):
                y, x = divmod(i, num_edge)
                pos_bullet = Vector(rect.topleft) + Vector((x + 1, y + 1)) * config["bullet.sep"]
                draw_circle(surface, self.addr2color((deg, i)), pos_bullet, config["bullet.radius"])
        # enlarge hovered-on bullet
        if self.id_hover_on is not None:
            deg, i = self.id2addr[self.id_hover_on]
            n = self.deg2num_bullets[deg]
            if n <= 9:
                offset = config["bullet.patterns"][n - 1][i] * config["bullet.sep"] // 2
                draw_circle(surface, self.addr2color((deg, i)),
                            Vector(deg2sp(deg)) + offset, config["bullet.radius"] * 1.5)
            else:
                v = self.expansion[deg]
                rect = v['rect']
                num_edge = v['num_edge']
                y, x = divmod(i, num_edge)
                pos_bullet = Vector(rect.topleft) + Vector((x + 1, y + 1)) * config["bullet.sep"]
                draw_circle(surface, self.addr2color((deg, i)), pos_bullet, config["bullet.radius"] * 1.5)
        # draw arrows
        # for arrow in self.ss.get_arrows(self.expansion):
        for arrow in self.ss.arrows:
            addr1, addr2 = self.id2addr[arrow.src], self.id2addr[arrow.tgt]
            draw_arrow(surface, self.addr2sp(addr1), self.addr2sp(addr2))
        if self.status == STATUS_ON_BULLET:
            if self.addr_mouse_down is not None and self.addr_mouse_down[1] is not None:
                draw_arrow(surface, self.addr2sp(self.addr_mouse_down), self.pos_current)
        # draw label
        if self.id_hover_on is not None:
            s = self.ss.get_bullet_by_id(self.id_hover_on).label
            draw_text(surface, s, (config["win_width"] // 2,
                                   config["win_height"] - config["margin_bottom"] // 2), self.font)

# 477
