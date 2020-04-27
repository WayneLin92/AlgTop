"""The Camera class manages the mapping between world coordinates and screen coordinates."""
import math
from algebras.mymath import Vector, clip
from GUI.pygame_wrapper import config


class Camera:
    def __init__(self):
        """`origin_sp` is the screen position of world origin.
        `unit_screen` is the screen length of world unit length."""
        self.origin_sp = Vector((config["margin_left"], config["win_height"] - config["margin_bottom"]))
        self.unit_screen = config["camera.unit_screen"]

    @staticmethod
    def flip(world_pos):
        return Vector((world_pos[0], -world_pos[1]))

    # setters -------------------
    def zoom(self, pos, rate):
        pos = Vector(pos)
        unit_screen1 = clip(self.unit_screen * rate, *config["camera.unit_screen.interval"])
        rate1 = unit_screen1 / self.unit_screen
        self.unit_screen = unit_screen1
        origin_sp1 = pos + (self.origin_sp - pos) * rate1
        x_interval = (config["win_width"] / 2 - config["x_max"] * self.unit_screen, config["win_width"] / 2)
        y_interval = (config["win_height"] / 2, config["win_height"] / 2 + config["y_max"] * self.unit_screen)
        x1 = clip(origin_sp1[0], *x_interval)
        y1 = clip(origin_sp1[1], *y_interval)
        self.origin_sp = Vector((x1, y1))

    def translation(self, rel):
        origin_sp1 = self.origin_sp + Vector(rel)
        x_interval = (config["win_width"] / 2 - config["x_max"] * self.unit_screen, config["win_width"] / 2)
        y_interval = (config["win_height"] / 2, config["win_height"] / 2 + config["y_max"] * self.unit_screen)
        x1 = clip(origin_sp1[0], *x_interval)
        y1 = clip(origin_sp1[1], *y_interval)
        self.origin_sp = Vector((x1, y1))

    # methods -----------------
    def wp2sp(self, world_pos):
        """Convert world position to screen position."""
        return self.origin_sp + self.flip(world_pos) * self.unit_screen

    def sp2wp(self, pos):
        """Convert screen position to world position."""
        return self.flip((Vector(pos) - self.origin_sp) / self.unit_screen)

    def degs_in_screen(self):
        """Return degrees which regions overlap with the visible screen."""
        bottom_left_screen = (0, config["win_height"])
        top_right_screen = (config["win_width"], 0)
        deg1 = self.sp2wp(bottom_left_screen)
        deg2 = self.sp2wp(top_right_screen)
        deg1 = (math.floor(deg1[0]), math.floor(deg1[1]))
        deg2 = (math.ceil(deg2[0]), math.ceil(deg2[1]))
        return ((x, y) for x in range(deg1[0], deg2[0] + 1) for y in range(deg1[1], deg2[1] + 1))
