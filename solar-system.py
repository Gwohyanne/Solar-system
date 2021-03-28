# -*- coding: utf-8 -*-
"""
Created on Sun Mar 28 12:41:46 2021

@author: djoan
"""

import numpy as np
import matplotlib.pyplot as plt
plt.rcParams['animation.ffmpeg_path'] = r'C:\Users\djoan\Downloads\ffmpeg-2021-03-28-git-8b2bde0494-full_build\bin\ffmpeg.exe'
import matplotlib.animation as animation
from astropy.time import Time
from astroquery.jplhorizons import Horizons

sim_start_date = "2000-01-01"     # simulating a solar system starting from this date
sim_duration = 5 * 365                # (int) simulation duration in days

class Object:                   # define the objects: the Sun, Earth, Mercury, etc
    def __init__(self, name, rad, color, r, v):
        self.name = name
        self.r    = np.array(r, dtype=np.float)
        self.v    = np.array(v, dtype=np.float)
        self.xs = []
        self.ys = []
        self.plot = ax.scatter(r[0], r[1], color=color, s=rad**2, edgecolors=None, zorder=10)
        self.line, = ax.plot([], [], color=color, linewidth=0.6)

class SolarSystem:
    def __init__(self, thesun):
        self.thesun = thesun
        self.planets = []
        self.time = None
        self.timestamp = ax.text(.03, .94, 'Date: ', color='w', transform=ax.transAxes, fontsize='x-large')
    def add_planet(self, planet):
        self.planets.append(planet)
    def evolve(self):           # evolve the trajectories
        dt = 1.0
        self.time += dt
        plots = []
        lines = []
        for p in self.planets:
            p.r += p.v * dt
            acc = -2.959e-4 * p.r / np.sum(p.r**2)**(3./2)  # in units of AU/day^2
            p.v += acc * dt
            p.xs.append(p.r[0])
            p.ys.append(p.r[1])
            p.plot.set_offsets(p.r[:2])
            p.line.set_xdata(p.xs)
            p.line.set_ydata(p.ys)
            plots.append(p.plot)
            lines.append(p.line)
        self.timestamp.set_text('Date: ' + Time(self.time, format='jd', out_subfmt='date').iso)
        return plots + lines + [self.timestamp]
    
class TerrestrialSystem:
    def __init__(self, earth, sun):
        self.earth = earth
        self.sun = sun
        self.planets = []
        self.time = None
        self.timestamp = ax.text(.03, .94, 'Date: ', color='w', transform=ax.transAxes, fontsize='x-large')
    def add_planet(self, planet):
        self.planets.append(planet)
    def evolve(self):           # evolve the trajectories
        dt = 1.0
        self.time += dt
        plots = []
        lines = []
        #trajectory of the sun
        self.sun.r += self.sun.v*dt
        acc = -2.959e-4 * self.sun.r / np.sum(self.sun.r**2)**(3./2)
        self.sun.v += acc * dt
        self.sun.xs.append(self.sun.r[0])
        self.sun.ys.append(self.sun.r[1])
        self.sun.plot.set_offsets(self.sun.r[:2])
        self.sun.line.set_xdata(self.sun.xs)
        self.sun.line.set_ydata(self.sun.ys)
        plots.append(self.sun.plot)
        lines.append(self.sun.line)
        #trajectory of the other planets
        for p in self.planets:
            p.r += p.v * dt
            acc = -2.959e-4 * p.r / np.sum(p.r**2)**(3./2)  # in units of AU/day^2
            p.v += acc * dt
            p.xs.append(p.r[0] + self.sun.r[0])
            p.ys.append(p.r[1] + self.sun.r[1])
            p.plot.set_offsets([p.r[0] + self.sun.r[0], p.r[1] + self.sun.r[1]])
            p.line.set_xdata(p.xs)
            p.line.set_ydata(p.ys)
            plots.append(p.plot)
            lines.append(p.line)
        self.timestamp.set_text('Date: ' + Time(self.time, format='jd', out_subfmt='date').iso)
        return plots + lines + [self.timestamp]

#set drawing settings
plt.style.use('dark_background')
fig = plt.figure(figsize=[6, 6])
ax = plt.axes([0., 0., 1., 1.], xlim=(-5., 5.), ylim=(-5., 5.))
ax.set_aspect('equal')
ax.axis('off')

#init earth and sun
earth = Object("Earth", 1., 'blue', [0, 0, 0], [0, 0, 0])
sunH = Horizons(id=3, location="@sun", epochs=Time(sim_start_date).jd, id_type='id').vectors()
sun = Object("Sun", 25., 'yellow', [np.double(-sunH[xi]) for xi in ['x', 'y', 'z']], 
                         [np.double(-sunH[vxi]) for vxi in ['vx', 'vy', 'vz']])
#init system
ts = TerrestrialSystem(earth, sun)
ts.time = Time(sim_start_date).jd
#add planets to system
mercury = Horizons(id=1, location="@sun", epochs=ts.time, id_type='id').vectors()
ts.add_planet(Object(1, 0.38, 'green', 
            [np.double(mercury[xi]) for xi in ['x', 'y', 'z']], 
            [np.double(mercury[vxi]) for vxi in ['vx', 'vy', 'vz']]))
venus = Horizons(id=2, location="@sun", epochs=ts.time, id_type='id').vectors()
ts.add_planet(Object(2, 0.95, 'orange', 
            [np.double(venus[xi]) for xi in ['x', 'y', 'z']], 
            [np.double(venus[vxi]) for vxi in ['vx', 'vy', 'vz']]))
mars = Horizons(id=4, location="@sun", epochs=ts.time, id_type='id').vectors()
ts.add_planet(Object(1, 0.53, 'red', 
            [np.double(mars[xi]) for xi in ['x', 'y', 'z']], 
            [np.double(mars[vxi]) for vxi in ['vx', 'vy', 'vz']]))
#legend
ax.text(-4.8, -4., "Sun", color='yellow', zorder=1000, ha='center', fontsize='small')
ax.text(-4.8, -4.2, "Mercury", color='green', zorder=1000, ha='center', fontsize='small')
ax.text(-4.8, -4.4, "Venus", color='orange', zorder=1000, ha='center', fontsize='small')
ax.text(-4.8, -4.6, "Earth", color='blue', zorder=1000, ha='center', fontsize='small')
ax.text(-4.8, -4.8, "Mars", color='red', zorder=1000, ha='center', fontsize='small')

def animate(i):
    return ts.evolve()
ani = animation.FuncAnimation(fig, animate, repeat=False, frames=sim_duration, blit=True, interval=20,)
plt.show()
 
writervideo = animation.FFMpegWriter(fps=60) 
ani.save('solar-system.mp4', writer=writervideo, dpi=150)