# -*- coding: utf-8 -*-
"""
Created on Tue Mar  2 14:09:21 2021

@author: dmg530
"""


class Data_Point:
    def __init__(self, intensity, cell_x_position, cell_y_position, wavelength, radius, angle):
        self.intensity = intensity
        self.cell_x_position = cell_x_position
        self.cell_y_position = cell_y_position
        self.wavelength = wavelength
        self.wavelength_angstrom = wavelength*0.1
        self.radius = radius
        self.angle = angle
        