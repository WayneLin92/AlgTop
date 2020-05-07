import math
from itertools import groupby
from algebras.mymath import Vector, clip, interpolation
from pygame import Rect
from GUI.camera import Camera
from GUI.pygame_wrapper import config, myroundv, Paint
from GUI.specseq import SpecSeq


class App:
    """The class for drawing

    self.expansion is a dictionary where the values are dictionaries
    """
    def __init__(self, ss: SpecSeq, surface):
        # for drawing
        self.status = "start"
        self.addr_mouse_down = None
        self.id_hover_on = None
        self.spread = None
        self.expansion = set()
        self.pos_current = None
        self.ctrl_down = False
        self.wait_for_update = True  # signals changes on screen
        self.paint = Paint(surface)

        self.ss = ss
        self.camera = Camera()
        self.id2addr = {}
        self.addr2id = {}
        self.deg2num_bullets = {}
        self.init_from_ss(ss)

        self.init_spread()

    def init_from_ss(self, ss):
        ss.bullets.sort(key=lambda p: p.deg)
        for deg, group in groupby(ss.bullets, key=lambda p: p.deg):
            i = -1
            for i, b in enumerate(group):
                self.id2addr[b.id] = (deg, i)
                self.addr2id[(deg, i)] = b.id
            self.deg2num_bullets[deg] = i + 1

    def init_spread(self):
        self.spread = set()
        for deg in self.deg2num_bullets:
            n = self.deg2num_bullets[deg]
            if n > 9:
                num_edge = math.ceil(math.sqrt(n))
                length_edge = self.get_bullet_sep() * (num_edge + 1)
                if length_edge / self.camera.unit_screen < config["auto_expand_rate"]:
                    self.spread.add(deg)

    def get_bullet_radius(self):
        radius = config["bullet.radius_world"] * self.camera.unit_screen
        return clip(radius, *config["bullet.radius_screen.interval"])

    def get_bullet_sep(self):
        sep = config["bullet.sep_world"] * self.camera.unit_screen
        return clip(sep, *config["bullet.sep_screen.interval"])

    def sp2deg(self, screen_pos):
        """Return the closest grid point to screen_pos."""
        wp = self.camera.sp2wp(screen_pos)
        return myroundv(wp)
    
    def deg2sp(self, deg):
        """Return the screen position of the deg."""
        return self.camera.wp2sp(deg)
    
    def offset_to_bullet(self, offset: Vector, n):
        """Need 1 <= n <= 9."""
        v = myroundv(offset / self.get_bullet_sep() * 2)
        if v in config["bullet.patterns.offset"][n - 1]:
            return config["bullet.patterns.offset"][n - 1][v]
        else:
            return None

    def exp_collide_point(self, screen_pos, deg=None):
        """Return if `screen_pos` is in an expansion square."""
        for d in (deg,) if deg else self.expansion:
            center = self.deg2sp(d)
            n = self.deg2num_bullets[d]
            num_edge = math.ceil(math.sqrt(n))
            length_edge = self.get_bullet_sep() * (num_edge + 1)
            size = Vector([length_edge, length_edge])
            rect = Rect(myroundv(center - size / 2), myroundv(size))
            if rect.collidepoint(*screen_pos):
                return d
        return None

    def expand(self, deg):
        self.expansion.add(deg)

    def exp_close(self, deg):
        self.expansion.remove(deg)

    def exp_close_all(self):
        self.expansion = set()

    def sp2addr(self, screen_pos):
        """* Return None if not on any bullet
        * Return (deg, None) if on folded bullet
        * Return (deg, index) if on some bullet"""
        deg = self.sp2deg(screen_pos)
        deg_expand = self.exp_collide_point(screen_pos) or (deg if deg in self.spread else None)
        if deg_expand:
            center = self.deg2sp(deg_expand)
            n = self.deg2num_bullets[deg_expand]
            num_edge = math.ceil(math.sqrt(n))
            length_edge = self.get_bullet_sep() * (num_edge + 1)
            size = Vector([length_edge, length_edge])
            offset = Vector(screen_pos) - (center - size / 2)
            x, y = myroundv(offset / self.get_bullet_sep())
            if 0 < x <= num_edge and 0 < y <= num_edge:
                index = x - 1 + num_edge * (y - 1)
                if index < n:
                    return deg_expand, index
            return None
        if deg in self.deg2num_bullets:
            n = self.deg2num_bullets[deg]
            if 0 < n <= 9:
                offset = Vector(screen_pos) - Vector(self.deg2sp(deg))
                index = self.offset_to_bullet(offset, n)
                if index is not None:
                    return deg, index
                else:
                    return None
            elif n > 9:
                if math.dist(screen_pos, self.deg2sp(deg)) < 2 * self.get_bullet_radius():
                    return deg, None
                else:
                    return None
        else:
            return None

    def addr2sp(self, addr):
        deg, index = addr
        if index is None:
            return self.deg2sp(deg)
        n = self.deg2num_bullets[deg]
        if n <= 9:
            offset = config["bullet.patterns"][n - 1][index] * self.get_bullet_sep() // 2
            return self.deg2sp(deg) + offset
        else:
            if deg in self.expansion or deg in self.spread:
                center = self.deg2sp(deg)
                n = self.deg2num_bullets[deg]
                num_edge = math.ceil(math.sqrt(n))
                length_edge = self.get_bullet_sep() * (num_edge + 1)
                size = Vector([length_edge, length_edge])
                y, x = divmod(index, num_edge)
                pos_bullet = (center - size / 2) + Vector([x + 1, y + 1]) * self.get_bullet_sep()
                return pos_bullet
            else:
                return self.deg2sp(deg)

    def addr2color(self, addr):  #
        return self.ss.get_bullet_by_id(self.addr2id[addr]).color

    def mouse_down(self, msg):
        if msg['button'] == 1:
            if self.status == "start":
                self.addr_mouse_down = self.sp2addr(msg['pos'])
                if self.addr_mouse_down is None:
                    self.status = "on_canvas"
                elif self.addr_mouse_down[1] is None:
                    if not self.ctrl_down:
                        self.exp_close_all()
                    self.expand(self.addr_mouse_down[0])
                    self.wait_for_update = True
                else:
                    self.status = "on_bullet"
        elif msg['button'] == 4:
            self.camera.zoom(msg['pos'], config["camera.zoom_rate"])
            self.init_spread()
            self.wait_for_update = True
        elif msg['button'] == 5:
            self.camera.zoom(msg['pos'], 1 / config["camera.zoom_rate"])
            self.init_spread()
            self.wait_for_update = True

    def mouse_up(self, msg):
        if self.status == "on_bullet":  # add the new arrow
            addr1 = self.addr_mouse_down
            addr2 = self.sp2addr(msg['pos'])
            if addr2 is not None and addr2[1] is not None:
                if addr1[0] != addr2[0]:
                    id1, id2 = self.addr2id[addr1], self.addr2id[addr2]
                    self.ss.add_arrow(id1, id2)
                    # print(len(self.ss.arrows))
                    self.id_hover_on = None
                self.status = "start"
            elif addr2 is None or addr2[0] == addr1[0]:
                self.status = "start"
            else:
                # noinspection PyUnresolvedReferences
                self.expand(addr2[0])
            self.wait_for_update = True
        elif self.status == "on_canvas_moving":
            self.status = "start"
        elif self.status == "on_canvas":
            self.exp_close_all()
            self.status = "start"
            self.wait_for_update = True

    def mouse_motion(self, msg):
        self.pos_current = msg['pos']
        if self.status == "on_bullet":
            deg = self.addr_mouse_down[0]
            if deg in self.expansion and not self.exp_collide_point(msg['pos'], deg):
                self.exp_close(deg)
            self.wait_for_update = True
        elif self.status in {"on_canvas", "on_canvas_moving"}:
            self.camera.translation(msg['rel'])
            self.status = "on_canvas_moving"
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
            pass
        elif msg['key'] == 127:  # Del -> initialize every thing
            pass
        elif msg['key'] == 8:  # Backspace -> delete an arrow
            pass
        elif msg['key'] == 306:  # Ctrl
            self.ctrl_down = True
        elif msg['unicode'] == 'r':  # R -> initialize the view setting
            self.camera.__init__()
            self.wait_for_update = True

    def key_up(self, msg):
        # vir_pos = vs.sp2wp(msg['pos'])
        if msg['key'] == 306:  # Ctrl
            self.ctrl_down = False

    def bg(self):
        """Draw the grid lines."""
        self.paint.clear_screen()
        n_label_step = math.ceil(config["axis_text_sep_screen"] / self.camera.unit_screen)
        for i in range(0, config["y_max"] + 1):
            left = self.camera.wp2sp((0, i))
            right = self.camera.wp2sp((config["x_max"], i))
            self.paint.draw_line(config["grid_color"], left, right, 1)
            if i % n_label_step == 0:
                text_pos = left - Vector([30, 0])
                if text_pos[0] < 16:
                    text_pos = Vector((16, text_pos[1]))
                self.paint.draw_text(str(i), text_pos, config["axis_numbers_color"])
        for i in range(0, config["x_max"] + 1):
            bottom = self.camera.wp2sp((i, 0))
            top = self.camera.wp2sp((i, config["y_max"]))
            self.paint.draw_line(config["grid_color"], top, bottom, 1)
            if i % n_label_step == 0:
                text_pos = bottom + Vector([0, 30])
                if text_pos[1] > config["win_height"] - 16:
                    text_pos = Vector((text_pos[0], config["win_height"] - 16))
                self.paint.draw_text(str(i), text_pos, config["axis_numbers_color"])

    def render(self):
        """Draw on the screen."""
        if not self.wait_for_update:
            return
        else:
            self.wait_for_update = False
        
        bullet_radius = self.get_bullet_radius()

        # draw background
        self.bg()

        # draw expansions
        for deg in self.expansion:
            if deg not in self.spread:
                n = self.deg2num_bullets[deg]
                num_edge = math.ceil(math.sqrt(n))
                length_edge = self.get_bullet_sep() * (num_edge + 1)
                center = self.deg2sp(deg)
                size = Vector([length_edge, length_edge])
                rect = Rect(myroundv(center - size / 2), myroundv(size))
                self.paint.draw_rect(config["pen_color"], rect, 1)

        # draw lines
        lines = sorted(((addr1 := self.id2addr[line.src_id], addr2 := self.id2addr[line.tgt_id],
                         math.dist(addr1[0], addr2[0]), line.color)
                        for line in self.ss.lines), key=lambda _l: _l[2])
        for addr1, addr2, dist, color in reversed(lines):
            if dist < 3:
                line_color = color
            else:
                line_color = interpolation(3 / dist, color, config['bg_color'])
            self.paint.draw_line(line_color, self.addr2sp(addr1), self.addr2sp(addr2))

        # draw arrows
        for arrow in self.ss.arrows:
            addr1, addr2 = self.id2addr[arrow.src_id], self.id2addr[arrow.tgt_id]
            self.paint.draw_arrow(self.addr2sp(addr1), self.addr2sp(addr2))
        if self.status == "on_bullet":
            if self.addr_mouse_down is not None and self.addr_mouse_down[1] is not None:
                self.paint.draw_arrow(self.addr2sp(self.addr_mouse_down), self.pos_current)

        # draw bullets
        for deg in self.camera.degs_in_screen():
            if deg in self.deg2num_bullets:
                n = self.deg2num_bullets[deg]
                if n <= 9:
                    for i in range(n):
                        offset = config["bullet.patterns"][n - 1][i] * self.get_bullet_sep() / 2
                        self.paint.draw_circle(self.addr2color((deg, i)),
                                               Vector(self.deg2sp(deg)) + offset, bullet_radius)
                else:
                    num_edge = math.ceil(math.sqrt(n))
                    length_edge = self.get_bullet_sep() * (num_edge + 1)
                    center = self.deg2sp(deg)
                    size = Vector([length_edge, length_edge])
                    rect = Rect(myroundv(center - size / 2), myroundv(size))
                    if deg in self.spread:
                        for i in range(n):
                            y, x = divmod(i, num_edge)
                            pos_bullet = Vector(rect.topleft) + Vector((x + 1, y + 1)) * self.get_bullet_sep()
                            self.paint.draw_circle(self.addr2color((deg, i)), pos_bullet, bullet_radius)
                    elif deg not in self.expansion:
                        self.paint.draw_circle(config["pen_color"], self.deg2sp(deg), bullet_radius)
                        self.paint.draw_circle(config["pen_color"], self.deg2sp(deg), bullet_radius * 2, False)

        # draw bullets in expansions
        for deg in self.expansion:
            if deg not in self.spread:
                n = self.deg2num_bullets[deg]
                num_edge = math.ceil(math.sqrt(n))
                length_edge = self.get_bullet_sep() * (num_edge + 1)
                center = self.deg2sp(deg)
                size = Vector([length_edge, length_edge])
                top_left = myroundv(center - size / 2)
                for i in range(n):
                    y, x = divmod(i, num_edge)
                    pos_bullet = Vector(top_left) + Vector((x + 1, y + 1)) * self.get_bullet_sep()
                    self.paint.draw_circle(self.addr2color((deg, i)), pos_bullet, bullet_radius)

        # enlarge hovered-on bullet
        if self.id_hover_on is not None:
            addr = self.id2addr[self.id_hover_on]
            sp = self.addr2sp(addr)
            color = self.ss.get_bullet_by_id(self.id_hover_on).color
            self.paint.draw_circle(color, sp, bullet_radius * 1.5)

        # draw texts
        for text_p in self.ss.texts:
            if self.camera.is_in_screen(text_p.deg):
                s = text_p.text
                pos = self.camera.wp2sp(Vector(text_p.deg) + text_p.offset)
                self.paint.draw_text(s, pos, text_p.color)

        # draw label
        if self.id_hover_on is not None:
            s = self.ss.get_bullet_by_id(self.id_hover_on).label
            self.paint.draw_text(s, (config["win_width"] // 2, config["win_height"] - config["margin_bottom"] // 2))
