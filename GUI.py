#!/usr/bin/env python3
import tkinter as tk

import matplotlib
matplotlib.use('TkAgg')
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg as Tkcanvas
from matplotlib.backend_bases import MouseEvent
from matplotlib.figure import Figure as Fig

import numpy as np
import json
from functools import partial

import os
os.makedirs('adj/', exist_ok = True)
if not os.path.exists('adj/nets'):
    with open('adj/nets','w') as f: json.dump([],f)
with open('adj/nets','r') as f: saved = json.load(f)

# - - - - -  master frame
root = tk.Tk()
root.focus_set()


tk.Label(root, text="Generator for adjacency matrices", font = ('Courier', '16', 'bold')).pack(side='top')
tk.Label(root, text="Bachelor's Thesis by Thomas Axmann ", font = ('Courier', '12', 'italic')).pack(side='top')

def quit(event=None):
    root.destroy()
root.bind('<Escape>', quit)




# - - - - - var and info frames 
left = tk.Frame(root, width = 500)
left.pack(side='left')


## --- info frame
infotext = """
Hello!
This is my generator for problem related data.
You can choose the state of net editing by clicking on the button beyond or by hitting keys:

:m: -- dragg&drop points
:c: -- connect points
:p: -- pick the points that will build up the subsystem.

Hitting of one of these will always lead you to the "edit mode" where right clicking on connection lines or points will induce their vanishing. Once a problem is suitably proposed you can save it via the menu. There will be a dialogue opened in the shell which can always be aborted by typing q or quit.

On the right side projects that were there befor launching this program are listed. Click to load, hit "override" to save changes (destroying the previous state).
"""
info = tk.Frame(left)
info.pack(side='top')
text = tk.Label(info, justify = 'left', wraplength=280, text = infotext)
text.pack(side='top')


## --- var frame, simulation variables int(><.get()) to use them 
var_frame = tk.Frame(left, borderwidth=2, relief='sunken')
var_frame.pack(side='top')
tk.Label(var_frame, text="Simulation Variables").pack(side='top')

#... sites, updated from plot
sites = tk.StringVar()
sites.set('0')
sites_frame=tk.Frame(var_frame)
sites_frame.pack(side='top')
tk.Label(sites_frame, text="Sites", width=6).pack(side='left')
tk.Label(sites_frame, width=4, textvariable=sites).pack(side='left')

#... particles, not visualized in plot but fiercely important
particles = tk.StringVar()
particles.set('1')
particles_frame=tk.Frame(var_frame)
particles_frame.pack(side='top')
tk.Label(particles_frame, text="Particles", width=6).pack(side='left')
tk.Entry(particles_frame, width=4, textvariable=particles).pack(side='left')

#... button for whether the particle are fermions or hard core bosons
FERM = False
fermions = tk.StringVar()
fermions.set('Hard-core Bosons')
def ferm_cb(event=None):
    global FERM
    FERM = not FERM
    if FERM:
        fermions.set("Fermions")
    else: fermions.set("Hard-core Bosons")

tk.Button(var_frame, textvariable=fermions, command = ferm_cb, width=16).pack(side='top')


## --- state of net editing
state = tk.Frame(left)
state.pack(side='top')
EDIT_MOD, MOVE_MOD, CONNECT_MOD, CHOOSE_MOD = True, False, False, False
state_var = tk.StringVar()
state_var.set('Editing Site Count')

def state_cb(event=None):
    global MOVE_MOD, CONNECT_MOD, EDIT_MOD, CHOOSE_MOD
    CHOOSE_MOD, EDIT_MOD, MOVE_MOD, CONNECT_MOD = EDIT_MOD, MOVE_MOD, CONNECT_MOD, CHOOSE_MOD 
    if EDIT_MOD:
        state_var.set("Editing Site Count")
    elif MOVE_MOD:
        state_var.set("Moving Sites")
    elif CONNECT_MOD:
        state_var.set("Connecting Sites")
    else:
        state_var.set("Choosing Subsystem")

def state_to_m(event=None):
    global EDIT_MOD, MOVE_MOD, CONNECT_MOD, CHOOSE_MOD
    CONNECT_MOD = False
    CHOOSE_MOD = False
    MOVE_MOD = not MOVE_MOD
    EDIT_MOD = not MOVE_MOD
    if MOVE_MOD:
        state_var.set("Moving Sites")
    else:
        state_var.set("Editing Site Count")

def state_to_c(event=None):
    global EDIT_MOD, MOVE_MOD, CONNECT_MOD, CHOOSE_MOD
    MOVE_MOD = False
    CHOOSE_MOD = False
    CONNECT_MOD = not CONNECT_MOD
    EDIT_MOD = not CONNECT_MOD
    if CONNECT_MOD:
        state_var.set("Connecting Sites")
    else:
        state_var.set("Editing Site Count")

def state_to_ch(event=None):
    global EDIT_MOD, MOVE_MOD, CONNECT_MOD, CHOOSE_MOD
    MOVE_MOD = False
    CONNECT_MOD = False
    CHOOSE_MOD = not CHOOSE_MOD
    EDIT_MOD = not CHOOSE_MOD
    if CHOOSE_MOD:
        state_var.set("Choosing Subsystem")
    else:
        state_var.set("Editing Site Count")


root.bind('m', state_to_m)
root.bind('c', state_to_c)
root.bind('p', state_to_ch)
tk.Label(state, text='State of Net-editing:').pack(side='left')
tk.Button(state, textvariable = state_var, command=state_cb, width = 18).pack(side = 'left')




# - - - - -  plot frame
net_frame = tk.Frame(root)
net_frame.pack(side='left')
tk.Label(net_frame, text="Net").pack(side='top')
fig = Fig(figsize=(6,6)) # ordinary matplotlib figure
Tkcanvas(fig, master=net_frame).get_tk_widget().pack(side='top') # packed as tk widget
render = fig.canvas.draw # call this to redraw


## --- set up net
ax = fig.add_axes([0,0,1,1])
ax.set_xlim([0,300])
ax.set_ylim([0,300])

#############################################
# import matplotlib.image as im             # 
# dv = im.imread('darth.jpg')               #    import darth vader 
# ax.imshow(dv,extent=(0,300,0,300))        #
#############################################


## --- net bookkeeping
#... points
class point:
    def __init__(self, x, y=None):
        if isinstance(x, MouseEvent):
            x, y = int(x.xdata), int(x.ydata)
        self.pos = (x,y)
pts = []
moving = None
connecting = None

def nearest(x, y=None):
    """
    :rtype: int, None
    :val: index + 1 if |x,y - pts[i]| < 5
    """
    if isinstance(x, MouseEvent):
        x, y = int(x.xdata), int(x.ydata)
    mindist = 1000
    idx = None
    for i in range(len(pts)):
        xx, yy = pts[i].pos
        dist = np.sqrt((xx-x)**2 + (yy-y)**2)
        if dist < mindist:
            mindist = dist
            idx = i
    if mindist < 5 : return idx + 1
    return None

#... wires
wires = np.array([0])
con_line = ()
def w_add():
    global wires
    if not pts: return
    wires = np.vstack((wires, np.zeros((1, wires.shape[0]))))
    wires = np.hstack((wires, np.zeros((wires.shape[0], 1))))

def w_rem(idx):
    global wires
    if not (len(pts) - 1): return
    wires = (wires[[not i == idx for i in range(wires.shape[0])]])[:,[not i == idx for i in range(wires.shape[0])]]

#... subsystem
subs = []
def s_rem(idx):
    global subs
    if idx in subs:
        subs.remove(idx)
    for i in range(len(subs)):
        if subs[i] > idx: subs[i] -= 1


###
def draw():
    ax.cla()
    ax.set_xlim([0,300])
    ax.set_ylim([0,300])
    #ax.imshow(dv,extent=(0,300,0,300))
    sites.set(str(len(pts)))

    if con_line:
        ax.plot(*zip(*con_line), 'k-')
    m = wires.shape[0]
    for i in range(m - 1):
        n = i + 1 
        for j in range(n, m):
            if wires[i,j]:
                ax.plot(*zip(pts[i].pos, pts[j].pos), 'k-')
    for i, p in enumerate(pts):
        if i in subs:
            ax.plot(*(p.pos), 'ro')
        else: ax.plot(*(p.pos), 'ko')

    render()

#... mouse action (((keybinds
def m_press(event):

    if event.button == 1:
        if EDIT_MOD:
            w_add()
            pts.append(point(event))
        elif MOVE_MOD:
            global moving
            moving = nearest(event)
            return
        elif CONNECT_MOD:
            global connecting
            connecting = nearest(event)
            return
        elif CHOOSE_MOD:
            idx = nearest(event) - 1
            if idx in subs: subs.remove(idx)
            else: subs.append(idx)
            
    
    if event.button == 3:
        if EDIT_MOD:
            idx = nearest(event)
            if idx:
                w_rem(idx - 1)
                pts.pop(idx - 1)
                s_rem(idx - 1)

            else:
                m = wires.shape[0]
                for i in range(m - 1):
                    n = i + 1 
                    for j in range(n, m):
                        if wires[i,j]:
                            A, B, P = np.array(pts[i].pos), np.array(pts[j].pos), np.array((event.xdata, event.ydata))
                            l = np.dot(B-A, P-A)/np.dot(A-B,A-B)
                            R = P - (A*(1-l)+B*l)
                            if (0<l<1) and (25 > np.dot(R,R)):
                                wires[i,j] = 0.
                                wires[j,i] = 0.
                                draw()
                                return
                                                        
    draw()

def m_release(event):
    global moving, connecting, con_line 
    if moving: moving = None
    if connecting:
        idx = nearest(event)
        if idx:
            wires [idx - 1, connecting - 1] = 1. 
            wires [connecting - 1, idx - 1] = 1. 
        connecting = None
        con_line = None
        draw()

def m_motion(event):
    global moving, connecting
    if moving:
        pts[moving - 1].pos = int(event.xdata), int(event.ydata)
        draw()

    if connecting:
        global con_line
        con_line = (pts[connecting - 1].pos, (event.xdata, event.ydata))
        draw()


fig.canvas.mpl_connect('button_press_event', m_press)
fig.canvas.mpl_connect('button_release_event', m_release)
fig.canvas.mpl_connect('motion_notify_event', m_motion)




# - - - - - Menue # ><.add_separator() to have a ----- in the menu

menu = tk.Menu(root)


## --- save menu
# some constants
ex = ["quit", "q"]
bar = "_____________________________"

def save_cb(event=None, tosave=''):
    if tosave: s = tosave
    else:
        s = input('Enter name to save: ')
        while s in saved + ex:
            if s in ex:
                print('Now exiting.\n' + bar)
                return
            ok = input('Existing (override or try other): ')
            if not ok: break
            s = ok
        print('Ok. Saving...\n' + bar)

    if s not in saved: saved.append(s)
    with open('adj/nets','w') as f: json.dump(sorted(saved),f)
    save_points = [p.pos for p in pts]
    try: 
        save_net = tuple(map(tuple, wires))
    except TypeError:
        save_net = []
    save_vars = [subs, int(particles.get()), int(FERM)]
    with open('adj/' + s + '_pts','w') as f: json.dump(save_points, f)
    with open('adj/' + s + '_nt', 'w') as f: json.dump(save_net, f)
    with open('adj/' + s + '_vars', 'w') as f: json.dump(save_vars, f)

def override():
    if Current: save_cb(tosave=Current)


Current = ''
def load_cb(event=None, toload=''):
    if not saved:
        print(bar + '\nNo saved projects')
        return
    s = toload
    if not toload:
        print(bar + '\nSaved projects are:')
        print(*saved, sep='\n', end='\n')
        print(bar)
        s = input('Enter name to load: ')

    while s not in saved:
        if s in ex: 
            print("Now exiting.\n" + bar)
            return
        s = input('Not Existing. Try Other: ')
    
    with open('adj/' + s + '_pts', 'r') as f: load_points = json.load(f)
    with open('adj/' + s + '_nt', 'r') as f: load_net = json.load(f)
    with open('adj/' + s + '_vars', 'r') as f: load_vars = json.load(f)

    global pts, wires, subs, FERM, Current
    pts = [point(*p) for p in load_points]
    wires = np.array(load_net)
    if not load_net: wires = np.array([0])
    subs = load_vars[0]
    particles.set(str(load_vars[1]))
    FERM = bool(load_vars[2])
    ferm_cb(),ferm_cb()
    Current = s
    draw()

fio = tk.Menu(menu)
fio.add_command(label='Save (shell)', command=save_cb)
fio.add_command(label='Load (shell)', command=load_cb)
fio.add_command(label='Quit', command=quit)
menu.add_cascade(label='I/O', menu=fio)

root.config(menu=menu) 




# - - - - - interact frame
projects = tk.Frame(root)
tk.Label(projects,text='Projects').pack(side='top')
projects.pack(side='left')
for p in saved:
    tk.Button(projects, width=15,  text = p, command=partial(load_cb, toload=p)).pack(side='top')

tk.Frame(projects, width=100, height=1, bg='black').pack(side='top')
tk.Button(projects, text='Override', command=override).pack(side='top')





# - - - - -
root.mainloop()
# - - - - - cleanup 
