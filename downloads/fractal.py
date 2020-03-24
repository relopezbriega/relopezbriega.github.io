# -*- coding: utf-8 -*-
# Fractal functions collection

# imports
import pandas as pd
import numpy as np
import random
from scipy.stats import rv_discrete as rv
from scipy import integrate
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.colors import cnames
from matplotlib import animation
import math
import numba


# Fern chaos game 
class Fern:
  A = np.array([[[0.0,0.0],[0.00,0.16]],
            [[0.85,0.04],[-0.04,0.85]],
            [[0.2,-0.26],[0.23,0.22]],
            [[-0.15,0.28],[0.26,0.24]]])

  b = np.array([[0.00,0.00],
                [0.00,1.60],
                [0.00,1.60],
                [0.00,0.44]])

  #next_transform = rv(name='Fern',values=((0,1,2,3),(0.01,0.85,0.0,0.14)))
  next_transform = rv(name='Fern',values=((0,1,2,3),(0.01,0.85,0.07,0.07)))

  def simplify_axes(self, ax):
    ax.set_yticklabels([])
    ax.set_xticklabels([])
    ax.set_yticks([])
    ax.set_xticks([])
    ax.set_xlim(-2.3,2.8)
    ax.set_ylim(-0.1,10.1)


  def next_point(self, last_point):
    transform = self.next_transform.rvs()
    return np.dot(self.A[transform],last_point) + self.b[transform]

  def plot(self, iterations=1000):
    current = np.array([0.0, 0.0])

    fig = plt.figure()
    ax = fig.add_subplot(111)
    self.simplify_axes(ax)
      
    # Don't plot the first few iterations.
    for x in range(15):
      current = self.next_point(current)

    ax.plot(current[0], current[1], 'go', markersize=1)
    xs, ys = [], []
    line_fern, = ax.plot(xs, ys, 'go', markersize=1)

    for x in range(iterations):
      current = self.next_point(current)
      xs.append(current[0])
      ys.append(current[1])
      line_fern.set_data(xs, ys)
      
    plt.show()

# Sierpinski
class Sierpinski:
  def SierpinskiTriangle(self,a, b, c, iterations):
    '''
    Recursively generated Sierpinski Triangle. 
    '''
    if iterations == 0:
      # Fill the triangle with vertices a, b, c. 
      plt.fill([a[0], b[0], c[0]], [a[1], b[1], c[1]], 'g') 
      #plt.hold(True)
    else:
      # Recursive calls for the three subtriangles. 
      self.SierpinskiTriangle(a, (a + b) / 2., (a + c) / 2., iterations - 1) 
      self.SierpinskiTriangle(b, (b + a) / 2., (b + c) / 2., iterations - 1) 
      self.SierpinskiTriangle(c, (c + a) / 2., (c + b) / 2., iterations - 1)

  def plot(self):       
    a = np.array([0, 0])
    b = np.array([1, 0])
    c = np.array([0.5, np.sqrt(3)/2.])

    iterations = 0

    fig = plt.figure(figsize=(15,15))
    plt.subplot(2,3,1).set_title("Sierpinski Triangle (iterations = 0)")

    self.SierpinskiTriangle(a, b, c, iterations)

    plt.axis('equal')
    plt.axis('off')

    iterations = 1

    plt.subplot(2,3,2).set_title("Sierpinski Triangle (iterations = 1)")

    self.SierpinskiTriangle(a, b, c, iterations)

    plt.axis('equal')
    plt.axis('off')

    iterations = 2

    plt.subplot(2,3,3).set_title("Sierpinski Triangle (iterations = 2)")

    self.SierpinskiTriangle(a, b, c, iterations)

    plt.axis('equal')
    plt.axis('off')

    iterations = 3

    plt.subplot(2,3,4).set_title("Sierpinski Triangle (iterations = 3)")

    self.SierpinskiTriangle(a, b, c, iterations)

    plt.axis('equal')
    plt.axis('off')

    iterations = 4

    plt.subplot(2,3,5).set_title("Sierpinski Triangle (iterations = 4)")

    self.SierpinskiTriangle(a, b, c, iterations)

    plt.axis('equal')
    plt.axis('off')

    iterations = 5

    plt.subplot(2,3,6).set_title("Sierpinski Triangle (iterations = 5)")

    self.SierpinskiTriangle(a, b, c, iterations)

    plt.axis('equal')
    plt.axis('off')

    iterations = 6

    plt.figure(figsize=(25,25))

    self.SierpinskiTriangle(a, b, c, iterations)

    plt.title("Sierpinski Triangle (iterations = 6)")
    plt.axis('equal')
    plt.axis('off')
    plt.show()

# Vicsek
class Vicsek:
  def vicsek(self, a, b, c, d, iterations, offset=np.array([0,0])):
    ab = (a + b)/3.
    ba = 2*(a + b)/3.
    bc = (2*b + c)/3.
    cb = (b + 2*c)/3.
    dc = (c + 2*d)/3.
    cd = (2*c + d)/3.
    ad = (d + a)/3.
    da = 2*(d + a)/3.

    abd = 2*a/3. + (b + d)/3.
    bac = a + (2*b + d)/3.
    cbd = 4*a/3. + 2*(b + d)/3.
    dac = a + (b + 2*d)/3.


    if iterations == 0:
      plt.fill([ab[0]+offset[0],ba[0]+offset[0],bac[0]+offset[0],bc[0]+offset[0],
      cb[0]+offset[0],cbd[0]+offset[0],cd[0]+offset[0],dc[0]+offset[0],dac[0]+offset[0],
      da[0]+offset[0],ad[0]+offset[0],abd[0]+offset[0]],
      [ab[1]+offset[1],ba[1]+offset[1],bac[1]+offset[1],bc[1]+offset[1],
      cb[1]+offset[1],cbd[1]+offset[1],cd[1]+offset[1],dc[1]+offset[1],
      dac[1]+offset[1],da[1]+offset[1],ad[1]+offset[1],abd[1]+offset[1]],
        'saddlebrown')
      #plt.hold(True)
    else:
      abd_m =np.array([0,0])
      bac_m = bac - abd
      cbd_m = cbd - abd
      dac_m = dac - abd
      offset1= offset +abd
      self.vicsek(abd_m, bac_m, cbd_m, dac_m, iterations - 1,offset1)

      ab_m =np.array([0,0])
      ba_m = ba - ab
      bac_m = bac - ab
      abd_m = abd - ab
      offset2= offset +ab
      self.vicsek(ab_m, ba_m, bac_m, abd_m, iterations - 1,offset2)

      bac_m =np.array([0,0])
      bc_m = bc - bac
      cb_m = cb - bac
      cbd_m = cbd - bac
      offset4= offset +bac
      self.vicsek(bac_m, bc_m, cb_m, cbd_m, iterations - 1,offset4)

       
      dac_m = np.array([0, 0])
      cbd_m = cbd - dac
      cd_m = cd - dac
      dc_m = dc - dac
      offset6= offset +dac
      self.vicsek(dac_m, cbd_m, cd_m, dc_m, iterations - 1,offset6)

      ad_m = np.array([0, 0])
      abd_m = abd - ad
      dac_m = dac - ad
      da_m = da - ad
      offset8= offset +ad
      self.vicsek(ad_m, abd_m, dac_m, da_m, iterations - 1,offset8)

  def plot(self):
    a = np.array([0, 0])
    b = np.array([3, 0])
    c = np.array([3, 3])
    d = np.array([0, 3])

    fig = plt.figure(figsize=(10,10))

    iterations = 0

    plt.subplot(2,2,1).set_title("Vicsek Fractal (iterations = 0)")

    self.vicsek(a, b, c, d, iterations)

    plt.axis('equal')
    plt.axis('off')


    iterations = 1

    plt.subplot(2,2,2).set_title("Vicsek Fractal (iterations = 1)")

    self.vicsek(a, b, c, d, iterations)

    plt.axis('equal')
    plt.axis('off')


    iterations = 2

    plt.subplot(2,2,3).set_title("Vicsek Fractal (iterations = 2)")

    self.vicsek(a, b, c, d, iterations)

    plt.axis('equal')
    plt.axis('off')


    iterations = 3

    plt.subplot(2,2,4).set_title("Vicsek Fractal (iterations = 3)")

    self.vicsek(a, b, c, d, iterations)

    plt.axis('equal')
    plt.axis('off')

    iterations = 5

    plt.figure(figsize=(20,20))
    self.vicsek(a, b, c, d, iterations)
    #plt.hold(False)
    plt.title("Vicsek Fractal (iterations = 5)")
    plt.axis('equal')
    plt.axis('off')

    plt.show()

## Fractal Tree
class Tree:

  def drawTree(self, x1, y1, angle, depth):

    if depth:
      x2 = x1 + int(math.cos(math.radians(angle)) * depth * 10.0)
      y2 = y1 + int(math.sin(math.radians(angle)) * depth * 10.0)
      plt.plot([x1,x2],[y1,y2],'-',color='darkgreen',lw=3)
      self.drawTree(x2, y2, angle - 20, depth - 1)
      self.drawTree(x2, y2, angle + 20, depth - 1)

  def plot(self):
    plt.figure(figsize=(20,3))

    plt.subplot(1,5,1)
    depth=1
    self.drawTree(300, 550, 90, depth)
    s="Fractal Tree (iterations = %i)" %depth
    plt.title(s)
    plt.axis('off')

    plt.subplot(1,5,2)
    depth=2
    self.drawTree(300, 550, 90, depth)
    s="Fractal Tree (iterations = %i)" %depth
    plt.title(s)
    plt.axis('off')

    plt.subplot(1,5,3)
    depth=3
    self.drawTree(300, 550, 90, depth)
    s="Fractal Tree (iterations = %i)" %depth
    plt.title(s)
    plt.axis('off')

    plt.subplot(1,5,4)
    depth=4
    self.drawTree(300, 550, 90, depth)
    s="Fractal Tree (iterations = %i)" %depth
    plt.title(s)
    plt.axis('off')

    plt.subplot(1,5,5)
    depth=5
    self.drawTree(300, 550, 90, depth)
    s="Fractal Tree (iterations = %i)" %depth
    plt.title(s)
    plt.axis('off')

    plt.figure(figsize=(25,20))
    depth=12
    self.drawTree(300, 550, 90, depth)
    #self.drawTree(300, 550, 45, depth)
    #self.drawTree(300, 550, 135, depth)
    s="Fractal Tree (iterations = %i)" %depth
    plt.title(s)
    plt.axis('off')

    plt.show()


# Dragon
class Dragon:
  def L_system(self, level, initial_state, trgt, rplcmnt, trgt2, rplcmnt2):
    state = initial_state
   
    for counter in range(level):
      state2 = ''
      for character in state:
        if character == trgt:
          state2 += rplcmnt
        elif character == trgt2:
          state2 += rplcmnt2
        else:
          state2 += character
      state = state2
    return state

  def L( angle,coords,jump):
    return angle + math.radians(45)
  def R( angle,coords,jump):
    return angle - math.radians(45)
  def l( angle,coords,jump):
    return angle + math.radians(90)
  def r( angle,coords,jump):
    return angle - math.radians(90)

  def F( angle, coords, jump):
    coords.append(
        (coords[-1][0] + jump * math.cos(angle),
          coords[-1][1] + jump * math.sin(angle)))
    return angle

  def G( angle,coords,jump):
    coords.append(
      (coords[-1][0] + cosin[angle],
          coords[-1][1] +sines[angle]))
    return angle

  decode = dict(L=L, R=R, F=F, 
                G=G,l=l,r=r)

  def dragon(self, steps, length=200, startPos=(0,0)):
    starting= 'R'*steps+'FX'
    pathcodes = self.L_system(steps,  starting, 'X', 'XlYFl', 'Y', 'rFXrY')
    jump = float(length) / (2 ** steps)
    coords = [startPos]
    angle = 0
    for move in pathcodes:
      if move == 'F' or move =='r' or move== 'l' or move == 'R':
          angle= self.decode[move](angle,coords,jump)
    return coords


  def plot(self):  
    totalwidth=100
    iterations = 0

    fig = plt.figure(figsize=(17,10))
    points = self.dragon(iterations,totalwidth,(-totalwidth/2,0))

    plt.subplot(2,3,1).set_title("Dragon Curve (iterations = 0)")
        
    plt.plot([p[0] for p in points], [p[1] for p in points], '-',lw=3)#,lw=5)
    plt.axis('equal')
    plt.axis('off')

    iterations = 1

    plt.subplot(2,3,2).set_title("Dragon Curve (iterations = 1)")

    points = self.dragon(iterations,totalwidth,(-totalwidth/2,0))
    plt.plot([p[0] for p in points], [p[1] for p in points], '-',lw=3)#,lw=5)
    plt.axis('equal')
    plt.axis('off')

    iterations = 2

    plt.subplot(2,3,3).set_title("Dragon Curve (iterations = 2)")

    points = self.dragon(iterations,totalwidth,(-totalwidth/2,0))
    plt.plot([p[0] for p in points], [p[1] for p in points], '-',lw=3)#,lw=5)
    plt.axis('equal')
    plt.axis('off')

    iterations = 3

    plt.subplot(2,3,4).set_title("Dragon Curve (iterations = 3)")

    points = self.dragon(iterations,totalwidth,(-totalwidth/2,0))
    plt.plot([p[0] for p in points], [p[1] for p in points], '-',lw=3)#,lw=5)
    plt.axis('equal')
    plt.axis('off')

    iterations = 4

    plt.subplot(2,3,5).set_title("Dragon Curve (iterations = 4)")

    points = self.dragon(iterations,totalwidth,(-totalwidth/2,0))
    plt.plot([p[0] for p in points], [p[1] for p in points], '-',lw=3)#,lw=5)
    plt.axis('equal')
    plt.axis('off')

    iterations = 5

    plt.subplot(2,3,6).set_title("Dragon Curve (iterations = 5)")

    points = self.dragon(iterations,totalwidth,(-totalwidth/2,0))
    plt.plot([p[0] for p in points], [p[1] for p in points], '-',lw=3)#,lw=5)
    plt.axis('equal')
    plt.axis('off')

    iterations = 20

    plt.figure(figsize=(20,16))

    points = self.dragon(iterations,totalwidth,(-totalwidth/2,0))
    plt.plot([p[0] for p in points], [p[1] for p in points], '-',color='darkslateblue')#,lw=5)
    plt.axis('equal')
    plt.axis('off')

    plt.title("Dragon Curve (iterations = 20)")

    plt.show()

# Levy
class Levy:
  def L_system(self, level, initial_state, trgt, rplcmnt, trgt2, rplcmnt2):
    state = initial_state
     
    for counter in range(level):
      state2 = ''
      for character in state:
        if character == trgt:
          state2 += rplcmnt
        elif character == trgt2:
          state2 += rplcmnt2
        else:
          state2 += character
      state = state2
    return state

  def L(angle,coords,jump):
    return angle + math.radians(45)
  def R(angle,coords,jump):
    return angle - math.radians(45)
  def l(angle,coords,jump):
    return angle + math.radians(45)
  def r(angle,coords,jump):
    return angle - math.radians(45)

  def F(angle, coords, jump):
    coords.append(
      (coords[-1][0] + jump * math.cos(angle),
        coords[-1][1] + jump * math.sin(angle)))
    return angle

  def G(angle,coords,jump):
    coords.append(
      (coords[-1][0] + cosin[angle],
        coords[-1][1] +sines[angle]))
    return angle

  decode = dict(L=L, R=R, F=F, G=G,l=l,r=r)

  def levyc(self, steps, length=200, startPos=(0,0)):
    starting= 'R'*steps+'FX'
    pathcodes = self.L_system(steps,  'F', 'F', 'rFllFr', '', '')
    jump = float(length) / (math.sqrt(2) ** steps)
    coords = [startPos]
    angle = 0
    for move in pathcodes:
      if move == 'F' or move =='r' or move== 'l' or move == 'R':
        angle= self.decode[move](angle,coords,jump)
    return coords

  def plot(self):
    totalwidth=100
    iterations = 0

    fig = plt.figure(figsize=(17,10))
    points = self.levyc(iterations,totalwidth,(-totalwidth/2,0))

    plt.subplot(2,3,1).set_title("Levy's C Curve (iterations = 0)")
        
    plt.plot([p[0] for p in points], [p[1] for p in points], '-',lw=3,color='darkgreen')
    plt.axis('equal')
    plt.axis('off')

    iterations = 1

    plt.subplot(2,3,2).set_title("Levy's C Curve (iterations = 1)")

    points = self.levyc(iterations,totalwidth,(-totalwidth/2,0))
    plt.plot([p[0] for p in points], [p[1] for p in points], '-',lw=3,color='darkgreen')
    plt.axis('equal')
    plt.axis('off')

    iterations = 2

    plt.subplot(2,3,3).set_title("Levy's C Curve (iterations = 2)")

    points = self.levyc(iterations,totalwidth,(-totalwidth/2,0))
    plt.plot([p[0] for p in points], [p[1] for p in points], '-',lw=3,color='darkgreen')
    plt.axis('equal')
    plt.axis('off')

    iterations = 3

    plt.subplot(2,3,4).set_title("Levy's C Curve (iterations = 3)")

    points = self.levyc(iterations,totalwidth,(-totalwidth/2,0))
    plt.plot([p[0] for p in points], [p[1] for p in points], '-',lw=3,color='darkgreen')
    plt.axis('equal')
    plt.axis('off')

    iterations = 4

    plt.subplot(2,3,5).set_title("Levy's C Curve (iterations = 4)")

    points = self.levyc(iterations,totalwidth,(-totalwidth/2,0))
    plt.plot([p[0] for p in points], [p[1] for p in points], '-',lw=3,color='darkgreen')
    plt.axis('equal')
    plt.axis('off')

    iterations = 5

    plt.subplot(2,3,6).set_title("Levy's C Curve (iterations = 5)")

    points = self.levyc(iterations,totalwidth,(-totalwidth/2,0))
    plt.plot([p[0] for p in points], [p[1] for p in points], '-',lw=3,color='darkgreen')
    plt.axis('equal')
    plt.axis('off')

    iterations = 20

    plt.figure(figsize=(20,16))

    points = self.levyc(iterations,totalwidth,(-totalwidth/2,0))
    plt.plot([p[0] for p in points], [p[1] for p in points], '-',color='olivedrab')
    plt.axis('equal')
    plt.axis('off')

    plt.title("Levy's C Curve (iterations = 20)")

    plt.show()

# Koch
class Koch:

  def koch(self, a,b,iterations):
    a1=a[0]
    a2=a[1]
    
    b1=b[0]
    b2=b[1]
    
    theta = np.arctan((b2-a2)/(b1-a1))
    length = np.sqrt((a1-b1)**2+(a2-b2)**2)
    
    c1 = (2*a1+b1)/3.
    c2 = (2*a2+b2)/3.
    c = [c1,c2]
    
    d1 = (a1+2*b1)/3.
    d2 = (a2+2*b2)/3.
    d = [d1,d2]
    
    if c1 >= a1:
      m1 = c1 + (length/3.)*math.cos(theta+math.pi/3.)
      m2 = c2 + (length/3.)*math.sin(theta+math.pi/3.)
    else:
      m1 = c1 + (length/3.)*math.cos(theta-2*math.pi/3.)
      m2 = c2 + (length/3.)*math.sin(theta-2*math.pi/3.)
    m = [m1,m2]
    
    c = np.array(c)
    d = np.array(d)
    m = np.array(m)
    
    points = []
    
    if iterations == 0:
        points.extend([a,b])
    elif iterations == 1:
        points.extend([a, c, m, d, b])
    else:
        points.extend(self.koch(a,c,iterations-1))
        points.extend(self.koch(c,m,iterations-1))
        points.extend(self.koch(m,d,iterations-1))
        points.extend(self.koch(d,b,iterations-1))  
                        
    return points


  def plot(self):
    fig = plt.figure(figsize=(15,5))

    plt.subplot(2,3,1).set_title("Koch Line (iterations = 0)")
    points = self.koch(a=np.array([0, 0]),b=np.array([1,0]),iterations=0)
    ptsx=[]
    ptsy=[]
    for i in range(len(points)):
        ptsx.append(points[i][0])
        ptsy.append(points[i][1])
    plt.plot(ptsx, ptsy, '-')
    plt.axis('equal')
    plt.axis('off')

    plt.subplot(2,3,2).set_title("Koch Line (iterations = 1)")
    points = self.koch(a=np.array([0, 0]),b=np.array([1,0]),iterations=1)
    ptsx=[]
    ptsy=[]
    for i in range(len(points)):
        ptsx.append(points[i][0])
        ptsy.append(points[i][1])
    plt.plot(ptsx, ptsy, '-')
    plt.axis('equal')
    plt.axis('off')

    plt.subplot(2,3,3).set_title("Koch Line (iterations = 2)")
    points = self.koch(a=np.array([0, 0]),b=np.array([1,0]),iterations=2)
    ptsx=[]
    ptsy=[]
    for i in range(len(points)):
        ptsx.append(points[i][0])
        ptsy.append(points[i][1])
    plt.plot(ptsx, ptsy, '-')
    plt.axis('equal')
    plt.axis('off')

    plt.subplot(2,3,4).set_title("Koch Line (iterations = 3)")
    points = self.koch(a=np.array([0, 0]),b=np.array([1,0]),iterations=3)
    ptsx=[]
    ptsy=[]
    for i in range(len(points)):
        ptsx.append(points[i][0])
        ptsy.append(points[i][1])
    plt.plot(ptsx, ptsy, '-')
    plt.axis('equal')
    plt.axis('off')


    plt.subplot(2,3,5).set_title("Koch Line (iterations = 4)")
    points = self.koch(a=np.array([0, 0]),b=np.array([1,0]),iterations=4)
    ptsx=[]
    ptsy=[]
    for i in range(len(points)):
        ptsx.append(points[i][0])
        ptsy.append(points[i][1])
    plt.plot(ptsx, ptsy, '-')
    plt.axis('equal')
    plt.axis('off')

    plt.subplot(2,3,6).set_title("Koch Line (iterations = 5)")
    points = self.koch(a=np.array([0, 0]),b=np.array([1,0]),iterations=5)
    ptsx=[]
    ptsy=[]
    for i in range(len(points)):
        ptsx.append(points[i][0])
        ptsy.append(points[i][1])
    plt.plot(ptsx, ptsy, '-')
    plt.axis('equal')
    plt.axis('off')

    h = np.sqrt(3)/2.
    a = np.array([0, 0])
    b = np.array([1, 0])
    c = np.array([0.5, h])

    iterations = 7

    points1 = self.koch(a=np.array([0, 0]),b=np.array([1,0]),iterations=iterations)
    points2 = self.koch(a=np.array([1, 0]),b=np.array([0.5,-h]),iterations=iterations)
    points3 = self.koch(a=np.array([0.5, -h]),b=np.array([0,0]),iterations=iterations)

    points = []
    for i in range(len(points1)):
        points.append(np.array(points1[i]))
    for i in range(len(points2)):
        points.append(np.array(points2[i]))
    for i in range(len(points3)):
        points.append(np.array(points3[i]))

    ptsx=[]
    ptsy=[]
    for i in range(len(points)):
        ptsx.append(points[i][0])
        ptsy.append(points[i][1])

    plt.figure(figsize=(25,25))

    plt.title("Koch Triangle (iterations = 7)")
    plt.plot(ptsx, ptsy, '-')
    plt.fill(ptsx, ptsy, color='lightcyan',alpha=0.7)
    plt.axis('equal')
    plt.axis('off')
    plt.show()

# Hilbert
class Hilbert:
  def hilbert(self, x0, y0, xi, xj, yi, yj, n,points):
    if n <= 0:
      X = x0 + (xi + yi)/2
      Y = y0 + (xj + yj)/2
      points.append((X,Y))
    else:
      self.hilbert(x0, y0, yi/2, yj/2, xi/2, xj/2, n - 1,points)
      self.hilbert(x0 + xi/2, y0 + xj/2,  xi/2, xj/2, yi/2, yj/2, n - 1,points)
      self.hilbert(x0 + xi/2 + yi/2, y0 + xj/2 + yj/2, xi/2, xj/2, yi/2, yj/2, n - 1,points)
      self.hilbert(x0 + xi/2 + yi,   y0 + xj/2 + yj,  -yi/2,-yj/2,-xi/2,-xj/2, n - 1,points)
      return points

  def plot(self):
    a = np.array([0, 0])
    b = np.array([1, 0])
    c = np.array([1, 1])
    d = np.array([0, 1])
    ab = (a + b)/2.
    bc = (b + c)/2.
    cd = (c + d)/2.
    ad = (d + a)/2.
    aab = (a + ab)/2.
    bba = (b + ab)/2.
    aad = (a + ad)/2.
    dda = (d + ad)/2.
    ccb = (c + bc)/2.
    bbc = (b + bc)/2.
    ccd = (c + cd)/2.
    ddc = (d + cd)/2.

    iterations = 1

    fig = plt.figure(figsize=(17,17))
    plt.subplot(2,3,1).set_title("Hilbert Curve (iterations = 1)")

    points = self.hilbert(0.0, 0.0, 1.0, 0.0, 0.0, 1.0, iterations,[])
    plt.plot([p[0] for p in points], [p[1] for p in points], '-',lw=3,color='darkmagenta')

    plt.plot([a[0],b[0],c[0],d[0],a[0]],[a[1],b[1],c[1],d[1],a[1]],'k-',lw=1)
    plt.plot([ab[0],cd[0]],[ab[1],cd[1]],'k--',lw=1)
    plt.plot([ad[0],bc[0]],[ad[1],bc[1]],'k--',lw=1)
    plt.plot([b[0],c[0]],[b[1],c[1]],'k-',lw=3)

    plt.axis('equal')
    plt.axis('off')

    iterations = 2

    plt.subplot(2,3,2).set_title("Hilbert Curve (iterations = 2)")

    points = self.hilbert(0.0, 0.0, 1.0, 0.0, 0.0, 1.0, iterations,[])
    plt.plot([p[0] for p in points], [p[1] for p in points], '-',lw=3,color='darkmagenta')#,lw=5)

    plt.plot([a[0],b[0],c[0],d[0],a[0]],[a[1],b[1],c[1],d[1],a[1]],'k-',lw=1)
    plt.plot([ab[0],cd[0]],[ab[1],cd[1]],'k--',lw=1)
    plt.plot([ad[0],bc[0]],[ad[1],bc[1]],'k--',lw=1)
    plt.plot([b[0],c[0]],[b[1],c[1]],'k-',lw=3)

    plt.axis('equal')
    plt.axis('off')

    iterations = 3

    plt.subplot(2,3,3).set_title("Hilbert Curve (iterations = 3)")

    points = self.hilbert(0.0, 0.0, 1.0, 0.0, 0.0, 1.0, iterations,[])
    plt.plot([p[0] for p in points], [p[1] for p in points], '-',lw=3,color='darkmagenta')#,lw=5)

    plt.plot([a[0],b[0],c[0],d[0],a[0]],[a[1],b[1],c[1],d[1],a[1]],'k-',lw=1)
    plt.plot([ab[0],cd[0]],[ab[1],cd[1]],'k--',lw=1)
    plt.plot([ad[0],bc[0]],[ad[1],bc[1]],'k--',lw=1)
    plt.plot([b[0],c[0]],[b[1],c[1]],'k-',lw=3)
    plt.plot([aab[0],ddc[0]],[aab[1],ddc[1]],'k--',lw=1)
    plt.plot([bba[0],ccd[0]],[bba[1],ccd[1]],'k--',lw=1)
    plt.plot([aad[0],bbc[0]],[aad[1],bbc[1]],'k--',lw=1)
    plt.plot([dda[0],ccb[0]],[dda[1],ccb[1]],'k--',lw=1)

    plt.axis('equal')
    plt.axis('off')

    iterations = 4

    plt.subplot(2,3,4).set_title("Hilbert Curve (iterations = 4)")

    points = self.hilbert(0.0, 0.0, 1.0, 0.0, 0.0, 1.0, iterations,[])
    plt.plot([p[0] for p in points], [p[1] for p in points], '-',lw=3,color='darkmagenta')#,lw=5)

    plt.axis('equal')
    plt.axis('off')

    iterations = 5

    plt.subplot(2,3,5).set_title("Hilbert Curve (iterations = 5)")

    points = self.hilbert(0.0, 0.0, 1.0, 0.0, 0.0, 1.0, iterations,[])
    plt.plot([p[0] for p in points], [p[1] for p in points], '-',lw=3,color='darkmagenta')#,lw=5)

    plt.axis('equal')
    plt.axis('off')
    iterations = 6

    plt.subplot(2,3,6).set_title("Hilbert Curve (iterations = 6)")

    points = self.hilbert(0.0, 0.0, 1.0, 0.0, 0.0, 1.0, iterations,[])
    plt.plot([p[0] for p in points], [p[1] for p in points], '-',lw=3,color='darkmagenta')#,lw=5)

    plt.axis('equal')
    plt.axis('off')

    fig = plt.figure(figsize=(25,25))

    points = self.hilbert(0.0, 0.0, 1.0, 0.0, 0.0, 1.0, iterations,[])
    plt.plot([p[0] for p in points], [p[1] for p in points], '-',lw=3,color='darkmagenta')#,lw=5)
    plt.title("Hilbert Curve (iterations = 7)")
    plt.axis('equal')
    plt.axis('off')

    plt.show()

# Lorentz Attactor
class Lorentz:
  def __init__(self, ntrajectories):
    self.ntrajectories = ntrajectories
  
  def lorentz_deriv(self, x_y_z, t0, sigma=10., beta=8./3, rho=28.0):
    """Compute the time-derivative of a Lorentz system."""
    (x, y, z) = x_y_z
    return [sigma * (y - x), x * (rho - z) - y, x * y - beta * z]

  def plot(self):
    # Choose random starting points, uniformly distributed from -15 to 15
    np.random.seed(1)
    x0 = -15 + 30 * np.random.random((self.ntrajectories, 3))

    # Solve for the trajectories
    t = np.linspace(0, 4, 1000)
    x_t = np.asarray([integrate.odeint(self.lorentz_deriv, x0i, t)
                      for x0i in x0])

    # Set up figure & 3D axis for animation
    fig = plt.figure()
    ax = fig.add_axes([0, 0, 1, 1], projection='3d')
    ax.axis('off')

    # choose a different color for each trajectory
    colors = plt.cm.jet(np.linspace(0, 1, self.ntrajectories))

    # set up lines and points
    lines = sum([ax.plot([], [], [], '-', c=c)
                for c in colors], [])
    pts = sum([ax.plot([], [], [], 'o', c=c)
              for c in colors], [])

    # prepare the axes limits
    ax.set_xlim((-25, 25))
    ax.set_ylim((-35, 35))
    ax.set_zlim((5, 55))

    # set point-of-view: specified by (altitude degrees, azimuth degrees)
    ax.view_init(30, 0)

    # initialization function: plot the background of each frame
    def init():
        for line, pt in zip(lines, pts):
            line.set_data([], [])
            line.set_3d_properties([])

            pt.set_data([], [])
            pt.set_3d_properties([])
        return lines + pts

    # animation function.  This will be called sequentially with the frame number
    def animate(i):
        # we'll step two time-steps per frame.  This leads to nice results.
        i = (2 * i) % x_t.shape[1]

        for line, pt, xi in zip(lines, pts, x_t):
            x, y, z = xi[:i].T
            line.set_data(x, y)
            line.set_3d_properties(z)

            pt.set_data(x[-1:], y[-1:])
            pt.set_3d_properties(z[-1:])

        ax.view_init(30, 0.3 * i)
        fig.canvas.draw()
        return lines + pts

    # instantiate the animator.
    anim = animation.FuncAnimation(fig, animate, init_func=init,
                                  frames=500, interval=30, blit=True)

    # Save as mp4. This requires mplayer or ffmpeg to be installed
    anim.save('lorentz_attractor.mp4', fps=15, extra_args=['-vcodec', 'libx264'])

    plt.show()

# Julia Set
class Julia:
  def py_julia_fractal(self, z_re, z_im, j):
    '''Crea el grafico del fractal de Julia.'''
    for m in range(len(z_re)):
        for n in range(len(z_im)):
            z = z_re[m] + 1j * z_im[n]
            for t in range(256):
                z = z ** 2 - 0.05 + 0.68j
                if np.abs(z) > 2.0:
                    j[m, n] = t
                    break
  
  def plot(self):
    jit_julia_fractal = numba.jit(nopython=True)(self.py_julia_fractal)

    N = 1024
    j = np.zeros((N, N), np.int64)
    z_real = np.linspace(-1.5, 1.5, N)
    z_imag = np.linspace(-1.5, 1.5, N)
    jit_julia_fractal(z_real, z_imag, j)

    fig, ax = plt.subplots(figsize=(8, 8))
    ax.imshow(j, cmap=plt.cm.RdBu_r, extent=[-1.5, 1.5, -1.5, 1.5])
    ax.set_xlabel("$\mathrm{Re}(z)$", fontsize=18)
    ax.set_ylabel("$\mathrm{Im}(z)$", fontsize=18)
    plt.show()

def mandelbrot( h,w, maxit=40): 
  '''
  Returns an image of the Mandelbrot fractal of size (h,w).
  '''
  y,x = np.ogrid[ -1.4:1.4:h*1j, -2:0.8:w*1j ]
  #y,x = ogrid[ -0.25:0.25:h*1j, -1.8:-1.3:w*1j ]
  #y,x = ogrid[ -0.05:0.05:h*1j, -1.5:-1.4:w*1j ]
  c = x+y*1j
  z = c
  divtime = maxit + np.zeros(z.shape, dtype=int)
  for i in range(maxit):
    z  = z**2 + c
    diverge = z*np.conj(z) > 2**2            # who is diverging
    div_now = diverge & (divtime==maxit)  # who is diverging now
    divtime[div_now] = i  +100                # note when
    z[diverge] = 2                        # avoid diverging too much

  return divtime
    
def mandelbrot2( h,w,a=-2.,b=.8,c=-1.4,d=1.4, maxit=40): 
  '''
  Returns an image of the Mandelbrot fractal of size (h,w).
  '''
  y,x = np.ogrid[ c:d:h*1j, a:b:w*1j ]
  # y,x = ogrid[ -0.25:0.25:h*1j, -1.8:-1.3:w*1j ]
  #y,x = ogrid[ -0.05:0.05:h*1j, -1.5:-1.4:w*1j ]
  c = x+y*1j
  z = c
  divtime = maxit + np.zeros(z.shape, dtype=int)
  for i in range(maxit):
    z  = z**2 + c
    diverge = z*np.conj(z) > 2**2            # who is diverging
    div_now = diverge & (divtime==maxit)  # who is diverging now
    divtime[div_now] = i  +100                # note when
    z[diverge] = 2                        # avoid diverging too much

  return divtime
