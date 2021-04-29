import numpy as np
from numpy import *
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from tkinter import *
import tkinter as tk
import matplotlib
matplotlib.use("TkAgg")
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
from matplotlib.figure import Figure
import glob, os
import math

files = []

f = open("MVE.txt", "r")
lines = f.readlines()
N = lines.__len__()-4
title = lines[0]
plotStyle = lines[1]
OM = lines[3].split(',')
for x in range(OM.__len__()):
    OM[x] = eval(OM[x])

years = eval(lines[2]); days = 1; h = 86400*days; endt = 31557600*years;time = 7;

bodies = np.empty((N,10),dtype=np.object)
for x in range(lines.__len__()-4):
    info = list(map(int,lines[x+4].split(',')))
    bodies[x, 0] = info[0]
    bodies[x, 1] = info[1]
    bodies[x, 2] = info[2]
    bodies[x, 3] = info[3]
    bodies[x, 4] = info[4]
    bodies[x, 5] = info[5]
    bodies[x, 6] = endt
    bodies[x, 7] = h
    bodies[x, 8] = OM
    bodies[x, 9] = info[6]

G = 6.672*(10**-11);
BX = np.zeros((365*years,N)); BY = np.zeros((365*years,N)); BZ = np.zeros((365*years,N)); BVX = np.zeros((365*years,N)); BVY = np.zeros((365*years,N)); BVZ = np.zeros((365*years,N));
xvalues = [];yvalues = [];zvalues = [];vxvalues = [];vyvalues = [];vzvalues = [];
for x in range(N):
    xvalues.append(bodies[x,0])
    yvalues.append(bodies[x,1])
    zvalues.append(bodies[x,2])
    vxvalues.append(bodies[x,3])
    vyvalues.append(bodies[x,4])
    vzvalues.append(bodies[x,5])
BX[0,:] = xvalues; BVX[0,:] = vxvalues; BGX = np.zeros((365*years,N)); BGX1 = np.zeros((365*years,N));
BY[0,:] = yvalues; BVY[0,:] = vyvalues; BGY = np.zeros((365*years,N)); BGY1 = np.zeros((365*years,N));
BZ[0,:] = zvalues; BVZ[0,:] = vzvalues; BGZ = np.zeros((365*years,N)); BGZ1 = np.zeros((365*years,N));

potentialEnergy = np.zeros((365*years,N))
kineticEnergy = np.zeros((365*years,N))


class Application(tk.Frame):
    def __init__(self, master=None):
        super().__init__(master)
        self.pack(anchor="n",expand=1,fill="both")
        self.create_widgets()

    def create_widgets(self):
        self.bottomBar = tk.Frame(self).pack(side = "bottom",expand=1,fill="both")

        self.quit = tk.Button(self.bottomBar, text="QUIT", fg="red",
                              command=root.destroy)
        self.quit.pack(side="right",anchor = "e")

        self.dynamicPlot = tk.Button(self.bottomBar, text="DYNAMIC PLOT", fg="blue", command=self.runDynamicPlot)
        self.dynamicPlot.pack(side="bottom")

        self.fig = plt.figure(figsize = (10,10))

        self.canvas = FigureCanvasTkAgg(self.fig, self)
        self.canvas.draw()
        #self.canvas.get_tk_widget().pack(side=tk.LEFT, fill=tk.BOTH, expand=False)
        self.canvas._tkcanvas.pack(side='left',fill='x', expand=0,anchor='nw')
        self.ax = Axes3D(self.fig)
        self.runStaticPlot()

        self.plotChosen = tk.Button(self,text="Run Simulation", fg="blue", command = self.runPlot)
        self.plotChosen.pack(side="top",anchor="e")


        for file in glob.glob("*.txt"):
            files.append(file)
        self.v = tk.StringVar()
        self.v.set("Select Simulation")
        self.simPicker = tk.OptionMenu(self,self.v,*files)
        self.simPicker.pack(side="top",anchor="e")

        self.yearsEntry = Entry(self)
        self.yearsEntry.pack(side="top",anchor = "w")
        self.yearsEntry.insert(0,"Enter Number of Years")

        self.timeEntry = Entry(self)
        self.timeEntry.pack(side="top", anchor="w")
        self.timeEntry.insert(0, "Enter Time Stamp")

        self.dataBox = Text(self,height=40,width=50)
        self.dataBox.pack(side="top")
        self.dataBox.config(state=DISABLED)
        self.addData()

    def getYears(self):
        global years
        years = eval(self.yearsEntry.get())

    def getTime(self):
        global time
        time = eval(self.timeEntry.get())

    def runStaticPlot(self):
            i = 365*years-1
            for plots in range(N):
                self.ax.plot(BX[0:i, plots], BY[0:i, plots], BZ[0:i, plots])
                self.ax.scatter(BX[i, plots], BY[i, plots], BZ[i, plots], s=100)
            self.canvas.draw()

    def runDynamicPlot(self):
        fig1 = plt.figure(figsize=(10, 10))
        ax = Axes3D(fig1)
        fig1.canvas.set_window_title(title)
        plt.close(self.fig)

        #max_range = np.array([BX.max() - BX.min(), BY.max() - BY.min()]).max() / 2.0

        #mid_x = (BX.max() + BX.min()) * 0.5
        #mid_y = (BY.max() + BY.min()) * 0.5

        for i in range(0,365*years-1,time):
            for plots in range(N):
                ax.plot(BX[0:i, plots], BY[0:i, plots], BZ[0:i, plots])
                ax.scatter(BX[i, plots], BY[i, plots], BZ[i, plots], s=100)
            plt.pause(.001)
            ax.clear()
            #ax.set_xlim(mid_x - max_range, mid_x + max_range)
            #ax.set_ylim(mid_y - max_range, mid_y + max_range)
            if plt.fignum_exists(fig1.number)!= True:
                break
        i = 365*years-1
        for plots in range(N):
            ax.plot(BX[0:i, plots], BY[0:i, plots], BZ[0:i, plots])
            ax.scatter(BX[i, plots], BY[i, plots], BZ[i, plots], s=100)
            #ax.set_xlim(mid_x - max_range, mid_x + max_range)
            #ax.set_ylim(mid_y - max_range, mid_y + max_range)

    def runPlanetPlot(self):
        max_range = np.array([BX.max() - BX.min(), BY.max() - BY.min(), BZ.max() - BZ.min()]).max() / 2.0

        mid_x = (BX.max() + BX.min()) * 0.5
        mid_y = (BY.max() + BY.min()) * 0.5
        mid_z = (BZ.max() + BZ.min()) * 0.5
        self.ax.set_xlim(mid_x - max_range, mid_x + max_range)
        self.ax.set_ylim(mid_y - max_range, mid_y + max_range)
        self.ax.set_zlim(mid_z - max_range, mid_z + max_range)

        u = np.linspace(0, 2 * np.pi, 100)
        v = np.linspace(0, np.pi, 100)
        xs = []; ys = []; zs = []
        for i in range(N):
            xs.append(bodies[i,9] * np.outer(np.cos(u), np.sin(v)))
            ys.append(bodies[i,9] * np.outer(np.sin(u), np.sin(v)))
            zs.append(bodies[i,9] * np.outer(np.ones(np.size(u)), np.cos(v)))

        i = 365 * years - 1
        for plots in range(N):
            self.ax.plot(BX[0:i, plots], BY[0:i, plots], BZ[0:i, plots])
            self.ax.plot_surface(xs[plots] + BX[i, plots],ys[plots] + BY[i, plots],zs[plots] + BZ[i, plots], rstride=4, cstride=4)

        self.canvas.draw()

    def runDynamicPlanetPlot(self):
        print("Dynamic Planet Plot")

    def runPlot(self):
        sim = self.v.get()
        print(sim)
        self.getYears()
        self.getTime()
        initializeSim(sim)
        print(eval(plotStyle))
        integrate()
        if eval(plotStyle) == 1:
            self.runDynamicPlot()
        elif eval(plotStyle) == 2:
            self.ax.clear()
            self.runPlanetPlot()
        else:
            self.ax.clear()
            self.runStaticPlot()
        self.addData()

    def addData(self):
        calcPotential()
        calcKinetic()
        self.dataBox.config(state=NORMAL)
        self.dataBox.delete('1.0', END)
        for i in range(N):
            self.dataBox.insert(END,"Body " + str(i+1) + " max potential energy :" + str(potentialEnergy[:, i].max()) + "\n")
            self.dataBox.insert(END,"Body " + str(i + 1) + " min potential energy :" + str(potentialEnergy[:, i].min()) + "\n")
            self.dataBox.insert(END, "Body " + str(i + 1) + " max kinetic energy :" + str(kineticEnergy[:, i].max()) + "\n")
            self.dataBox.insert(END, "Body " + str(i + 1) + " min kinetic energy :" + str(kineticEnergy[:, i].min()) + "\n")
            self.dataBox.insert(END,"\n")
        self.dataBox.config(state=DISABLED)



def initializeSim(sim):
    global f, lines, N, title, plotStyle, OM, years, days, h, endt, time, bodies, G, BX, BY, BZ, BVX, BVY, BVZ
    global xvalues, yvalues, zvalues, vxvalues, vyvalues, vzvalues, BGX, BGY, BGZ, BGX1, BGY1, BGZ1

    f = open(sim, "r")
    lines = f.readlines()
    N = lines.__len__() - 4
    title = lines[0]
    plotStyle = lines[1]
    OM = lines[3].split(',')
    for x in range(OM.__len__()):
        OM[x] = eval(OM[x])

    days = 1;
    h = 86400 * days;
    endt = 31557600 * years;

    bodies = np.empty((N, 10), dtype=np.object)
    for x in range(lines.__len__() - 4):
        info = list(map(int, lines[x + 4].split(',')))
        bodies[x, 0] = info[0]
        bodies[x, 1] = info[1]
        bodies[x, 2] = info[2]
        bodies[x, 3] = info[3]
        bodies[x, 4] = info[4]
        bodies[x, 5] = info[5]
        bodies[x, 6] = endt
        bodies[x, 7] = h
        bodies[x, 8] = OM
        bodies[x, 9] = info[6]

    G = 6.672 * (10 ** -11);
    BX = np.zeros((365 * years, N));BY = np.zeros((365 * years, N));BZ = np.zeros((365 * years, N));BVX = np.zeros((365 * years, N));BVY = np.zeros((365 * years, N));BVZ = np.zeros((365 * years, N));
    xvalues = [];yvalues = [];zvalues = [];vxvalues = [];vyvalues = [];vzvalues = [];
    for x in range(N):
        xvalues.append(bodies[x, 0])
        yvalues.append(bodies[x, 1])
        zvalues.append(bodies[x, 2])
        vxvalues.append(bodies[x, 3])
        vyvalues.append(bodies[x, 4])
        vzvalues.append(bodies[x, 5])
    BX[0, :] = xvalues;BVX[0, :] = vxvalues;BGX = np.zeros((365 * years, N));BGX1 = np.zeros((365 * years, N));
    BY[0, :] = yvalues;BVY[0, :] = vyvalues;BGY = np.zeros((365 * years, N));BGY1 = np.zeros((365 * years, N));
    BZ[0, :] = zvalues;BVZ[0, :] = vzvalues;BGZ = np.zeros((365 * years, N));BGZ1 = np.zeros((365 * years, N));

initializeSim("MVE.txt")

def integrate():
    global f, lines, N, title, plotStyle, years, days, h, endt, time, bodies, G, BX, BY, BZ, BVX, BVY, BVZ
    global xvalues, yvalues, zvalues, vxvalues, vyvalues, vzvalues, BGX, BGY, BGZ, BGX1, BGY1, BGZ1
    for t in range(365*years-1):
        for g in range(N):
            body_vx = bodies[g, 3];body_vy = bodies[g, 4];body_vz = bodies[g, 5];body_OM = bodies[g, 8];
            for OM in range(size(body_OM)-1):
                if OM != g:
                    body_x = (BX[t, g] - BX[t, OM]);
                    body_y = (BY[t, g] - BY[t, OM]);
                    body_z = (BZ[t, g] - BZ[t, OM]);
                    BGX[t, g] = BGX[t, g] + ((-G * body_OM[OM] * body_x) / (math.sqrt(body_x ** 2 + body_y ** 2 + body_z ** 2)) ** 3);
                    BGY[t, g] = BGY[t, g] + ((-G * body_OM[OM] * body_y) / (math.sqrt(body_x ** 2 + body_y ** 2 + body_z ** 2)) ** 3);
                    BGZ[t, g] = BGZ[t, g] + ((-G * body_OM[OM] * body_z) / (math.sqrt(body_x ** 2 + body_y ** 2 + body_z ** 2)) ** 3);

        for p in range(N):
            body_OM = bodies[p, 8];
            BX[t + 1, p] = BX[t, p] + ((1 / 2) * ((BVX[t, p]) + (BVX[t, p] + (h * (BGX[t, p]))))) * h;
            BY[t + 1, p] = BY[t, p] + ((1 / 2) * ((BVY[t, p]) + (BVY[t, p] + (h * (BGY[t, p]))))) * h;
            BZ[t + 1, p] = BZ[t, p] + ((1 / 2) * ((BVZ[t, p]) + (BVZ[t, p] + (h * (BGZ[t, p]))))) * h;

        for g1 in range(N):
            body_vx = bodies[g1, 3];body_vy = bodies[g1, 4];body_vz = bodies[g1, 5];body_OM = bodies[g1, 8];
            for OM in range(size(body_OM)):
                if OM != g1:
                    body_x = (BX[t + 1, g1] - BX[t + 1, OM]);
                    body_y = (BY[t + 1, g1] - BY[t + 1, OM]);
                    body_z = (BZ[t + 1, g1] - BZ[t + 1, OM]);
                    BGX1[t, g1] = BGX1[t, g1] + ((-G * body_OM[OM] * body_x) / (math.sqrt(body_x ** 2 + body_y ** 2 + body_z ** 2)) ** 3);
                    BGY1[t, g1] = BGY1[t, g1] + ((-G * body_OM[OM] * body_y) / (math.sqrt(body_x ** 2 + body_y ** 2 + body_z ** 2)) ** 3);
                    BGZ1[t, g1] = BGZ1[t, g1] + ((-G * body_OM[OM] * body_z) / (math.sqrt(body_x ** 2 + body_y ** 2 + body_z ** 2)) ** 3);

        for v in range(N):
            body_vx = bodies[v, 3];body_vy = bodies[v, 4];body_vz = bodies[v, 5];body_OM = bodies[v, 8];
            BVX[t + 1, v] = BVX[t, v] + ((1 / 2) * (BGX[t, v] + BGX1[t, v])) * h;
            BVY[t + 1, v] = BVY[t, v] + ((1 / 2) * (BGY[t, v] + BGY1[t, v])) * h;
            BVZ[t + 1, v] = BVZ[t, v] + ((1 / 2) * (BGZ[t, v] + BGZ1[t, v])) * h;

integrate()

def calcPotential():
    global potentialEnergy
    potentialEnergy = np.zeros((365 * years, N))
    for t in range(365*years):
        for g in range(N):
            body_OM = bodies[g, 8];
            for OM in range(size(body_OM)-1):
                if OM != g:
                    body_x = (BX[t, g] - BX[t, OM]);
                    body_y = (BY[t, g] - BY[t, OM]);
                    body_z = (BZ[t, g] - BZ[t, OM]);
                    potentialEnergy[t, g] = potentialEnergy[t, g] + (((-G)*body_OM[g]*body_OM[OM])/math.sqrt(body_x ** 2 + body_y ** 2 + body_z ** 2))

def calcKinetic():
    global kineticEnergy
    kineticEnergy = np.zeros((365 * years, N))
    for t in range(365 * years):
        for g in range(N):
            kineticEnergy[t, g] = .5*(OM[g])*((math.sqrt(BVX[t, g] ** 2 + BVY[t, g] ** 2 + BVZ[t, g] ** 2))**2)
calcPotential()

root = tk.Tk()
root.title("N-Body Simulation")
root.geometry('1980x1020')
root.option_add("*Font", "TkDefaultFont 24")

app = Application(master=root)
app.mainloop()

