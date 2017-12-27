# -*- coding: utf-8 -*-
"""
    dxfsl
    ~~~~~

    renders geometry




"""
import numpy as np
import numpy.linalg as la

import matplotlib.pylab as pl
import matplotlib.patches as pch
import femagtools.dxfsl.geom as g
import logging
import io

logger = logging.getLogger(__name__)


#############################
#       PlotRenderer        #
#############################

class PlotRenderer(object):
    def __init__(self):
        self.fig = None
        self.background = '#eeeeee'
        pass

    def figure(self):
        if self.fig == None:
            self.fig = pl.figure(figsize=(9, 10), facecolor='lightblue')
            #self.fig.subplots_adjust(left=0.0, right=0.5, top=1.0, bottom=0.0)
        return self.fig
        
    def show_plot(self):
        if self.fig != None:
            pl.show()
            self.fig = None
        
    def arc(self, startangle, endangle, center, radius, color = 'blue'):
        self.ax.add_patch(pch.Arc(center, 2*radius, 2*radius,
                                  angle=0,
                                  theta1=startangle*180/np.pi,
                                  theta2=endangle*180/np.pi,
                                  color=color))

    def line(self, p1, p2, color='blue'):
        self.ax.add_line(pl.Line2D((p1[0], p2[0]),
                                   (p1[1], p2[1]), color=color))
    
    def circle(self, center, radius, color='blue'):
        self.ax.add_patch(pch.Circle(center, radius,
                                     fill=False, color=color))

    def point(self, p, col, color='blue'):
        pl.plot([p[0]], [p[1]], col, color=color)
        
    def render(self, geom, filename=None, **kwargs):
        draw_center = kwargs.get('draw_center', False)
        draw_corners = kwargs.get('draw_corners', False)
        draw_hull = kwargs.get('draw_hull', False)
        incl_bnd = kwargs.get('incl_bnd', False)

        fig = pl.figure()
        self.ax = fig.add_subplot(111)

        if draw_hull:
            for h in g.convex_hull(geom.g.nodes()):
                pl.plot([h[0]], [h[1]], 'ro')
                
        id = 0
        for area in geom.areas(incl_bnd):
            id += 1
            logger.info("area %d (size %d)", id, len(area))
            if len(area) > 1:
                for s in area:
                    s.render(self)
                    if draw_center:
                        c = s.center_of_connection()
                        pl.plot([c[0]], [c[1]], 'gs')

#                p = get_point_inside(area)
#                if p:
#                    pl.plot([p[0]], [p[1]], 'ro')
#                    self.ax.text(p[0], p[1], str(id))

        geom.remove_areas(incl_bnd)

        # draw all remaining edges and circles
        [circle.render(self) for circle in geom.circles()]
        [attr['object'].render(self)
         for e1, e2, attr in geom.g.edges(data=True)]
        if draw_center:
            for e1, e2, attr in geom.g.edges(data=True):
                c = attr['object'].center_of_connection()
                pl.plot([c[0]], [c[1]], 'gs')

        if draw_corners:
            for c in geom.find_corners(geom.g):
                pl.plot([c[0]], [c[1]], 'bs')

        geom.render_cut_lines(self)
        self.ax.axis('scaled', aspect='equal')
        if filename:
            pl.savefig(filename)
        else:
            pl.show()

    def print_convex_hull(self, nodes):
        points = list(nodes)
        points.sort()
        if len(points) <= 1:
            return points
    
        def cross(o, a, b):
            dx = a[0] - o[0], a[1] - o[1]
            dy = b[0] - o[0], b[1] - o[1]
            return dx[0] * dy[1] - dx[1] * dy[0]
    
        print("=== Print Convex Hull ===")
        lower = []
        for p in points:
            self.point(p, 'go')
            while len(lower) >= 2 and cross(lower[-2], lower[-1], p) <= 0:
                lower.pop()
            lower.append(p)
        print("lower={}".format(lower))

    def render_elements(self, geom, type, **kwargs):
        with_nodes = kwargs.get('with_nodes', False)
        with_hull = kwargs.get('with_hull', False)
        with_corners = kwargs.get('with_corners', False)
        single_view = kwargs.get('single_view', False)
        neighbors = kwargs.get('neighbors', False)
        draw_center = kwargs.get('draw_center', False)
        draw_inside = kwargs.get('draw_inside', False)
        title = kwargs.get('title', "")
        show  = kwargs.get('show', True)
        rows  = kwargs.get('rows', 1)
        cols  = kwargs.get('cols', 1)
        num  = kwargs.get('num', 1)
        points = kwargs.get('points', [])

        if show:
            rows = 1
            cols = 1
            num = 1
            
        mm = geom.minmax()
        
        if single_view:
            count=0            
            for e in geom.elements(type):
                print("Render Element {}".format(e))                
                if count == 0:
                    self.ax = self.figure().add_subplot(111)
                    self.ax.axis('scaled', aspect='equal')
                    
                e.render(self, 'blue', True)

                count += 1
                if count == 3:
                    self.point((mm[0]-5, mm[2]-5), 'ro', color='red')
                    self.point((mm[0]-5, mm[3]+5), 'ro', color='red')
                    self.point((mm[1]+5, mm[2]-5), 'ro', color='red')
                    self.point((mm[1]+5, mm[3]+5), 'ro', color='red')
                    pl.show()
                    count = 0
                    
            if count != 0:
                self.point((mm[0]-5, mm[2]-5), 'ro', color='red')
                self.point((mm[0]-5, mm[3]+5), 'ro', color='red')
                self.point((mm[1]+5, mm[2]-5), 'ro', color='red')
                self.point((mm[1]+5, mm[3]+5), 'ro', color='red')
                self.ax.axis('scaled', aspect='equal')
                self.show_plot()
            return
        
        self.ax = self.figure().add_subplot(rows, cols, num, axisbg=self.background)
        if len(title) > 0:
            self.ax.set_title(title, size=14)
        self.ax.grid(color='blue',linewidth=0.5)
                  
        for e in geom.elements(type):
            e.render(self, 'blue', with_nodes)
        
        if with_hull:
            for h in g.convex_hull(geom.virtual_nodes()):
                pl.plot([h[0]], [h[1]], 'ro')

        if with_corners:
            for c in geom.start_corners:
                self.point(c, 'rs')
            for c in geom.end_corners:
                self.point(c, 'rs')

        if draw_center:
            for area in geom.list_of_areas():
                for s in area.elements():
                    c = s.center_of_connection()
                    pl.plot([c[0]], [c[1]], 'gs')

        if draw_inside:
            for area in geom.list_of_areas():
                p = area.get_point_inside(geom)
                if p:
                    pl.plot([p[0]], [p[1]], 'ro', color='yellow')
            
        geom.render_cut_lines(self)
        geom.render_airgaps(self)
        if neighbors:
            geom.render_neighbors(self)
        
        if len(points) > 0:
            for p in points:
                self.point(p, 'ro', color='red')
                
        if geom.center:
            self.point(geom.center, 'ro', color='darkgreen')

        self.point((mm[0]-5, mm[2]-5), 'ro', color=self.background)
        self.point((mm[0]-5, mm[3]+5), 'ro', color=self.background)
        self.point((mm[1]+5, mm[2]-5), 'ro', color=self.background)
        self.point((mm[1]+5, mm[3]+5), 'ro', color=self.background)
        self.ax.axis('scaled', aspect='equal')        
        
        if show:
            self.show_plot()

    def render_areas(self, geom, **kwargs):
        with_nodes = kwargs.get('with_nodes', False)
        single_view = kwargs.get('single_view', False)
        
        if not single_view:
            fig = pl.figure()
            self.ax = fig.add_subplot(111)

        colors = ('red', 'green', 'blue', 'magenta', 'orange', 'grey', 'darkgreen')
        
        c = -1
        for area in geom.list_of_areas():
            if area.number_of_elements() > 1:
                c += 1
                if c >= len(colors):
                    c = 0
                if single_view:
                    print("*** AREA has {} elements ***".format(area.number_of_elements()))
                    fig = pl.figure()
                    self.ax = fig.add_subplot(111)

                area.render(self, colors[c], with_nodes)
                p = area.get_point_inside(geom)
                if p:
                    self.point(p, 'ro', color='magenta')
                
                if single_view:
                    self.ax.axis('scaled', aspect='equal')
                    pl.show()

        if not single_view:
            self.ax.axis('scaled', aspect='equal')
            pl.show()

    def render_area(self, area):
        fig = pl.figure()
        self.ax = fig.add_subplot(111)
       
        area.render(self, 'red', with_nodes=True)

        self.ax.axis('scaled', aspect='equal')
        pl.show()
            
    def draw_slot(self, id, slot, ax):
        poly = pch.Polygon(slot, fill=True, color='#ee82ee')

        center = (sum(list(zip(*slot))[0])/len(slot),
                  sum(list(zip(*slot))[1])/len(slot))
        ax.text(center[0], center[1], str(id))
        ax.add_patch(poly)
        
#############################
#       DumpRenderer        #
#############################

class DumpRenderer(object):
    def __init__(self, name):
        self.model = name
        
    def circle(self, center, radius):
        self.content.append(u'Crcle({}, {}, {})'.format(center[0],
                                                        center[1],
                                                        radius))

    def arc(self, startangle, endangle, center, radius):
        p1 = (center[0] + radius*np.cos(startangle),
              center[1] + radius*np.sin(startangle))
        p2 = (center[0] + radius*np.cos(endangle),
              center[1] + radius*np.sin(endangle))
        self.content.append(
            u"Arc({}, {}, {}, {}, {}, {})".format(
                p1[0], p1[1], p2[0], p2[1],
                center[0], center[1]))

    def line(self, p1, p2):
        self.content.append(
            u"Line({}, {}, {}, {})".format(
                p1[0], p1[1], p2[0], p2[1]))
        
    def render(self, geom):
        '''create file with elements'''
        incl_bnd = True
        self.content = []
        handled = set()
        
        # collect all areas
        for i, area in enumerate(geom.areas(incl_bnd)):
            self.content.append('-- {}'.format(i))
            for e in area:
                if not e in handled:
                    e.render(self)
                    handled.add(e)
        self.content.append('--')
        geom.remove_areas(incl_bnd)
        # and all remaining edges and circles
        for circle in geom.circles():
            circle.render(self)

        for e1, e2, attr in geom.g.edges(data=True):
            attr['object'].render(self)
    
        return self.content


#############################
#       NewFslRenderer      #
#############################

class NewFslRenderer(object):
    """a model that can created by FSL"""
    header = ['exit_on_error = false',
              'exit_on_end = false',
              'verbosity = 2',
              'pickdist = {}',
              'fm_nlin_colour=blue',
              'fm_nlin_color=fm_nlin_colour',
              'fm_nlin_mcvfile="M470P65A-AM-50Hz-96"',
              'fm_nlin_mcvfile_shft="ST52-3"',
              'fm_nlin_rlen = 100',
              'shaft_mat=1',
              'agndst = 0.5',
              'new_model_force("{}","Test")',
              'blow_up_wind(0,0, {}, {})']
    
    def __init__(self, name):
        self.model = name
        self.mirror_axis = None
        self.fm_nlin = None
        self.shaft = None

    def mirror_nodechains(self, p0, p1):
        self.mirror_axis = np.array((p0, p1)).ravel().tolist()

    def material(self, p0):
        self.fm_nlin = p0

    def circle(self, center, radius, color='blue'):
        num = int(2*np.pi*radius)
        if num < 8:
            num = 8
        circle = [u'cx, cy = {}, {}'.format(center[0],
                                            center[1]),
                  u'nc_circle_m(cx,cy,{}, {})'.format(radius, num),
                  u'create_mesh_se(cx, cy)\n']
        self.content += circle

    def arc(self, startangle, endangle, center, radius, color='blue'):
        num = 0
#        if self.nodedist > 0:
#            s = startangle
#            d = endangle - s
#            n = int(d/(2*np.pi))
#            if n < 0:
#                n -= 1
#            d -= n*2*np.pi
#            num = int(radius*d/self.nodedist + 1)
#            if num < 3 and radius*d > self.nodedist:
#                num = 3

        p1 = (center[0] + radius*np.cos(startangle),
              center[1] + radius*np.sin(startangle))
        p2 = (center[0] + radius*np.cos(endangle),
              center[1] + radius*np.sin(endangle))
        self.content.append(
            u"nc_circle_m({}, {}, {}, {}, {}, {}, {})".format(
                p1[0], p1[1], p2[0], p2[1],
                center[0], center[1], num))

    def line(self, p1, p2, color='blue'):
        num = 0
#        if self.nodedist > 0:
#            l = la.norm(np.asarray(p1)-p2)
#            num = int(l/self.nodedist + 1)
        self.content.append(
            u"nc_line({}, {}, {}, {}, {})".format(
                p1[0], p1[1], p2[0], p2[1], num))

    def sorted_elements(self, geom, inner=False):
        if inner:
            # Abstand von airgap Richtung Nullpunkt
            el_sorted = [(geom.max_radius - e.minmax_from_center((0.0, 0.0))[1], e)
                        for e in geom.elements(g.Shape)]
        else:
            # Abstand von airgap Richtung Aussen
            el_sorted = [(e.minmax_from_center((0.0, 0.0))[0] - geom.min_radius, e)
                        for e in geom.elements(g.Shape)]
            
        el_sorted.sort()
        return el_sorted
        
    def render(self, geom, filename, inner=False, outer=False):
        '''create file with nodechains'''

        self.content = []

        self.content.append(u'\n\nndt(agndst)\n')

        ndt_list = [(0.25, 1.25), (0.5, 2), (0.75, 3.0), (1.1, 3.0)]
        dist = geom.max_radius - geom.min_radius
            
        el_sorted = self.sorted_elements(geom, inner)        
        
#        self.content.append(u'-- all elements')
        x = 0
        for d, e in el_sorted:
            d_percent = d / dist
            if ndt_list[x][0] < d_percent:
                self.content.append(u'\nndt({}*agndst)\n'.format(ndt_list[x][1]))
                while ndt_list[x][0] < d_percent:
                    x += 1
#            self.content.append(u'-- d={} / dist={} == {}'.format(d, dist, d_percent))
            e.render(self)
            
        self.content.append(u'\n')

        for area in geom.list_of_areas():
            if area.number_of_elements() > 1:
                p = area.get_point_inside(geom)
                if p:
                    self.content.append(u"x0, y0 = {}, {}".format(p[0], p[1]))
                    self.content.append(u"point(x0, y0, red, 4)")
                    self.content.append(u"create_mesh_se(x0, y0)\n")
      
                if area.type == 1:
                    self.content.append(u'def_new_subreg(x0, y0, "iron", blue)')
                elif area.type == 2:
                    self.content.append(u'def_new_subreg(x0, y0, "windings", green)')
                elif area.type == 3:
                    self.content.append(u'def_new_subreg(x0, y0, "windings", red)')
#        ag = 0
#        if ag > 0:
#            self.content.append(
#                u'\n-- Airgap\nnc_circle_m({}, 0, 0, {}, 0, 0, 0)'.format(
#                    da2/2+2*ag/3, da2/2+2*ag/3))
#            self.content.append(
#                u'nc_circle_m({}, 0, 0, {}, 0, 0, 0)'.format(
#                    da1/2-2*ag/3, da1/2-2*ag/3))
#            self.content.append(
#                u'nc_line({}, 0, {}, 0, 0)'.format(
#                    da1/2, da2/2))
#            self.content.append(
#                u'nc_line(0, {}, 0, {}, 0)'.format(
#                    da1/2, da2/2))
#        else:
#            self.content.append(u'\nx0, y0 = {}, {}'. format(
#                dy2/2+0.1, 0.1))
#            self.content.append(u'create_mesh_se(x0, y0)')
#            self.content.append(u'def_new_subreg(x0, y0, "Rotor", green)')
#
#        if dy2 > 0:
#            self.content.append(u'\nx0, y0 = pr2c({}, {})'. format(
#                dy2/2, geom.alpha))
#            self.content.append(u'nc_line(0, 0, {}, 0, 0)'.format(
#                dy2/2))
#            self.content.append(u'nc_line(0, 0, x0, y0, 0)')
#            self.content.append(u'create_mesh_se(0.1, 0.1)')
#            self.content.append(
#                u'def_new_subreg(0.1, 0.1, "Shaft", lightgrey)')

        if geom.is_mirrored():
            self.content.append(u'-- mirror')
            self.content.append(u'mirror_nodechains({}, {}, {}, {})\n'.format(
                geom.mirror_corners[1][0], # max x1
                geom.mirror_corners[1][1], # max y1
                geom.mirror_corners[0][0], # min x2
                geom.mirror_corners[0][1]))# min y2
                
        # Winkel nach allfälligem Spiegeln
        self.content.append(u'alfa = {}\n'.format(geom.get_alfa()))
            
        self.content.append(u'-- rotate')
        self.content.append(u'x1, y1 = {}, {}'.format(
                geom.start_corners[0][0], geom.start_corners[0][1])) # min xy1
        self.content.append(u'x2, y2 = {}, {}'.format(
                geom.start_corners[1][0], geom.start_corners[1][1])) # max xy2
            
        if geom.is_mirrored():
            self.content.append(u'x3, y3 = pr2c(x2, alfa)')
            self.content.append(u'x4, y4 = pr2c(x1, alfa)')
        else:
            self.content.append(u'x3, y3 = {}, {}'.format(
                    geom.end_corners[1][0], geom.end_corners[1][1])) # max xy3
            self.content.append(u'x4, y4 = {}, {}'.format(
                    geom.end_corners[0][0], geom.end_corners[0][1])) # min xy4

        self.content.append(u'if m.{}_ncopies > 1 then'.format(geom.kind))
        self.content.append(u'  rotate_copy_nodechains(x1,y1,x2,y2,x3,y3,x4,y4,m.{}_ncopies-1)'
                            .format(geom.kind))
        self.content.append(u'  end')

        if self.fm_nlin:
            self.content.append(u'\nx0, y0 = {}, {}'. format(
                self.fm_nlin[0], self.fm_nlin[1]))

            mat = [u"if fm_nlin_mcvfile ~= 'dummy' then",
                   u"  if fm_nlin_mcvfile == 'air' then",
                   u"    def_mat_fm(x0,y0, 1.0, fm_nlin_rlen)",
                   u"  else",
                   u"    def_mat_fm_nlin(x0,y0, fm_nlin_colour, fm_nlin_mcvfile, fm_nlin_rlen)",
                   u"  end",
                   u"else",
                   u"  def_mat_fm(x0,y0, 1000.0, fm_nlin_rlen)",
                   u"end"]
            self.content.append(u'\n'.join(mat))

        if self.shaft:
            mat = [u'\nif shaft_mat==1 then',
                   u'  def_mat_fm_nlin(0.1,0.1,lightgrey,fm_nlin_mcvfile_shft,fm_nlin_rlen)',
                   u'end']
            self.content.append(u'\n'.join(mat))

        #f.write(u"\nadapt_window()\n")
        with io.open(filename, 'w', encoding='utf-8') as f:
            f.write('\n'.join(self.content))
                
    def render_main(self, motor, geom_inner, geom_outer,
                    filename, with_header=False):
        '''create main file'''

        self.content = []
        self.content.append(u'exit_on_error = false')
        self.content.append(u'exit_on_end = false')
        self.content.append(u'verbosity = 2')
        self.content.append(u'pickdist = 0.001\n')
        
        if geom_inner and geom_outer:
            n = [int(round(np.pi/x)) for x in [geom_outer.alfa,
                                               geom_inner.alfa]]
            num_poles = min(n)
            num_slots = max(n)
         
            self.content.append(u'm.num_poles = {}'.format(num_poles))
            self.content.append(u'm.tot_num_slot = {}'.format(num_slots))
            if num_poles == n[0]:
                npols_gen = int(geom_outer.get_symmetry_copies())+1
            else:
                npols_gen = int(geom_inner.get_symmetry_copies())+1
            self.content.append(u'm.npols_gen = {}'.format(npols_gen))
            self.content.append(
                u'm.num_slots = m.tot_num_slot*m.npols_gen/m.num_poles')
                
        self.content.append(u'da1 = {}'.format(
            2*geom_outer.min_radius))
        self.content.append(u'da2 = {}'.format(
            2*geom_inner.max_radius))
        self.content.append(u'ag = (da1 - da2)/2\n')
        self.content.append(u'agndst = 0.75')

        self.content.append(u'new_model_force("{}","Test")\n'.format(
            self.model))

        if geom_inner:
            ncopies = 'm.npols_gen'
            if num_poles == n[0]:
                ncopies = 'm.num_slots'
            self.content.append(
                u'm.{}_ncopies = {}'.format(
                    geom_inner.kind, ncopies))

        if geom_outer:
            ncopies = 'm.num_slots'
            if num_poles == n[0]:
                ncopies = 'm.npols_gen'
            self.content.append(
                u'm.{}_ncopies = {}'.format(
                    geom_outer.kind, ncopies))
            
        if geom_inner:
            self.content.append(
                u'dofile("{}_{}.fsl")'.format(
                    self.model, geom_inner.kind))
        if geom_outer:
            self.content.append(
                u'dofile("{}_{}.fsl")'.format(
                    self.model, geom_outer.kind))
            
        alfa = geom_inner.get_alfa() * (geom_inner.get_symmetry_copies()+1)
        alfa2 = geom_outer.get_alfa() * (geom_outer.get_symmetry_copies()+1)
        assert(np.isclose(alfa, alfa2))

        self.content.append(u'alfa = {}\n'.format(alfa))

        self.content.append(u'-- airgap')
        self.content.append(u'r1 = da2/2 + ag/3')
        self.content.append(u'x1, y1 = pr2c(r1, alfa)')
        self.content.append(u'n = r1*alfa/agndst + 1')
        self.content.append(u'nc_circle_m(r1, 0, x1, y1, 0.0, 0.0, n)\n')
        
        self.content.append(u'r2 = da2/2 + 2*ag/3')
        self.content.append(u'x2, y2 = pr2c(r2, alfa)')
        self.content.append(u'nc_circle_m(r2, 0, x2, y2, 0.0, 0.0, n)\n')

        self.content.append(u'nc_line(da2/2, 0, r2, 0, 0)')
        self.content.append(u'x3, y3 = pr2c(da2/2, alfa)')
        self.content.append(u'x4, y4 = pr2c(r2, alfa)')
        self.content.append(u'nc_line(x3, y3, x4, y4, 0, 0)\n')

        self.content.append(u'x0, y0 = pr2c(da2/2+ag/6, alfa/2)')
        self.content.append(u'create_mesh_se(x0, y0)')

        self.content.append(u'x0, y0 = pr2c(da2/2+3*ag/6, alfa/2)')
        self.content.append(u'create_mesh_se(x0, y0)\n')

        self.content.append(u'connect_models()\n')
           
        with io.open(filename, 'w', encoding='utf-8') as f:
            f.write('\n'.join(self.content))
