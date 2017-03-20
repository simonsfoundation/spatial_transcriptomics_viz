import numpy as np 
from numpy.linalg import norm
from jp_svg_canvas import canvas, cartesian_svg
from jp_gene_viz import proxy_html5_canvas
from jp_gene_viz import js_proxy
import contourist.lasso
from threading import Timer

canvas.load_javascript_support()

import ipywidgets as widgets
from IPython.display import display

SLIDE = 'SLIDE'
SLIDE_x = 'SLIDE_x'
SLIDE_y = 'SLIDE_y'
TSNE_x = 'TSNE_x'
TSNE_y = 'TSNE_y'
INDICATOR = 'INDICATOR'
NAME = "NAME"

#SLIDE_SIZE = 40
TSNE_RADIUS = 11
BACKGROUND_COLOR = "#FFFFDD"


def canvas_doodle(xmin, ymin, xmax, ymax, html_width=500, html_height=None, margin=50):
    canvas = proxy_html5_canvas.HTML5CanvasProxy()
    result = cartesian_svg.doodle(xmin, ymin, xmax, ymax, html_width, html_height, margin, svg=canvas)
    result.buffered = True
    return result

def color(i):
    c = (0x0d1f43 * (i+1)) % (0x1000000)
    return "#%06x" % c

class Data(object):

    def __init__(self, fn="spatial_transcriptomics_example_data_02222017.tsv"):
        self.fn = fn
        f = open(fn)
        headers = f.readline().split()
        data = []
        for line in f.readlines():
            D = {}
            for (name, value) in zip(headers, line.split()):
                try:
                    value = float(value)
                except ValueError:
                    pass
                D[name] = value
            data.append(D)
        name_to_datum = {}
        indicator_to_data = {}
        slide_to_data = {}
        for d in data:
            indicator = d[INDICATOR]
            slide = d[SLIDE]
            i_list = indicator_to_data.get(indicator, [])
            i_list.append(d)
            indicator_to_data[indicator] = i_list
            s_list = slide_to_data.get(slide, [])
            s_list.append(d)
            slide_to_data[slide] = s_list
            name = indicator + "_" + str(len(i_list))
            name_to_datum[name] = d
            d[NAME] = name
        self.data_list = data
        self.name_to_datum = name_to_datum
        self.indicator_to_data = indicator_to_data
        self.slide_to_data = slide_to_data
        self.indicator_to_color = {ind: color(i) for (i, ind) in enumerate(indicator_to_data.keys())}
        self.slides = list(sorted(slide_to_data.keys()))
        self.slide_height = max(d[SLIDE_y] for d in data)
        self.slide_width = max(d[SLIDE_x] for d in data)

    def data_names(self):
        return [d[NAME] for d in self.data_list]

    def data_tsne_positions(self):
        return [(d[TSNE_x], d[TSNE_y]) for d in self.data_list]

class Slide(object):

    #radius = SLIDE_SIZE * 0.01
    #selected_radius = radius * 2
    xy = None  # slide location assigned by container
    
    def __init__(self, name, data, drawing=None):
        self.name = name
        self.data = data
        self.radius = self.data.slide_width * 0.01
        self.drawing = drawing

    def draw_on(self, dx, dy, data, drawing=None, selected=None, indicator=None):
        if drawing is None:
            drawing = self.drawing
        assert drawing is not None
        #drawing.empty()
        drawing.rect(None, dx, dy, data.slide_width, data.slide_height, BACKGROUND_COLOR)
        slide_name = self.name
        data_list = self.data.slide_to_data[slide_name]
        indicator_to_color = self.data.indicator_to_color
        for d in data_list:
            if d[SLIDE] == slide_name:
                name = d[NAME]
                ind = d[INDICATOR]
                radius = self.radius
                if selected is not None:
                    radius = self.radius * 0.5
                    if name in selected:
                        radius = self.radius * 2
                elif indicator is not None:
                    if ind == indicator:
                        radius = self.radius * 2
                    else:
                        radius = self.radius * 0.5
                color = indicator_to_color [d["INDICATOR"]]
                x = d["SLIDE_x"] + dx
                y = d["SLIDE_y"] + dy
                drawing.circle(name, x, y, radius, color)
        drawing.text(None, dx, dy + data.slide_height, slide_name)
        drawing.flush()

class TSNE(object):
    
    radius = TSNE_RADIUS * 0.005

    def __init__(self, data, drawing=None):
        self.data = data
        self.drawing = drawing

    def draw_on(self, dx, dy, data, drawing=None, indicator=None, slide=None, selected=None):
        if drawing is None:
            drawing = self.drawing
        assert drawing is not None
        drawing.empty()
        drawing.rect("background", dx-TSNE_RADIUS, dy-TSNE_RADIUS, dx+2*TSNE_RADIUS, dy+2*TSNE_RADIUS, BACKGROUND_COLOR)
        data = self.data
        ns, xs, ys, rs, cs = [],[],[],[],[]
        indicator_to_color = data.indicator_to_color
        for d in data.data_list:
            ind = d[INDICATOR]
            color = indicator_to_color[ind]
            x = d[TSNE_x] + dx
            y = d[TSNE_y] + dy
            n = d[NAME]
            xs.append(x)
            ys.append(y)
            cs.append(color)
            radius = self.radius
            if selected is not None:
                radius = self.radius * 0.5
                if d[NAME] in selected:
                    radius = self.radius * 2
            elif slide is not None:
                radius = self.radius * 0.5
                if d[SLIDE] == slide:
                    radius = self.radius * 2
            elif indicator is not None:
                radius = self.radius * 0.5
                if ind == indicator:
                    radius = self.radius * 2
            rs.append(radius)
            ns.append(n)
        drawing.circles(ns, xs, ys, rs, cs)
        drawing.flush()


class Slides(object):

    small_width = 80
    big_width = small_width * 6
    row_size = 8
    
    def __init__(self, data=None, slide_columns=8):
        self.row_size = slide_columns
        if data is None:
            data = Data()
        self.data = data
        slides = data.slides
        slide_displays = []
        slide_name_to_controller = {}
        slide_controllers = []
        slide_rows = int(len(slides) / slide_columns)
        if slide_rows * slide_columns < len(slides):
            slide_rows += 1
        self.slide_margin = margin = 5
        row_size = 2*margin + slide_rows * (data.slide_height + margin)
        col_size = 2*margin + (slide_columns) * (data.slide_width + margin)
        slides_width = 2 * self.big_width
        #slide_drawing = cartesian_svg.doodle(-1, -1, col_size+1, row_size+1, slides_width, margin=2)
        #print "rows, cols", slide_rows, slide_columns
        #print "slides dimensions", len(slides), col_size, row_size, "html width", slides_width
        slide_drawing = canvas_doodle(-1, -1, col_size+1, row_size+1, slides_width, margin=2)
        self.slide_drawing = slide_drawing
        #slide_drawing.enable_events("click", self.slide_drawing_click)
        count = 0
        for slide_name in slides:
            count += 1
            #if count > 30:
            #    break
            controller = Slide(slide_name, data, slide_drawing)
            slide_controllers.append(controller)
            assert slide_name not in slide_name_to_controller
            slide_name_to_controller[slide_name] = controller
        self.slide_controllers = slide_controllers
        controller_rows = []
        available_controllers = list(slide_controllers)
        while available_controllers:
            row_list = available_controllers[:self.row_size]
            available_controllers = available_controllers[self.row_size:]
            controller_rows.append(row_list)
        self.controller_rows = controller_rows
        self.slide_name_to_controller = slide_name_to_controller
        #tsne_drawing = cartesian_svg.doodle(-TSNE_RADIUS-1, -TSNE_RADIUS-1, TSNE_RADIUS+1, TSNE_RADIUS+1, self.big_width, margin=4)
        tsne_drawing = canvas_doodle(-TSNE_RADIUS-1, -TSNE_RADIUS-1, TSNE_RADIUS+1, TSNE_RADIUS+1, self.big_width, margin=4)
        self.tsne_drawing = tsne_drawing
        #tsne_events = "click mouseover mouseout"
        #tsne_events = "click"
        #print "enabling tsne events", tsne_events
        #tsne_drawing.enable_events(tsne_events, self.tsne_event_callback)
        #big_slide_drawing = cartesian_svg.doodle(-1, -1, data.slide_width+1, data.slide_height+1, self.big_width, margin=4)
        big_slide_drawing = canvas_doodle(-1, -1, data.slide_width+1, data.slide_height+1, self.big_width, margin=4)
        self.big_slide_drawing = big_slide_drawing
        # make dialog
        w = js_proxy.ProxyWidget()
        e = w.element()
        w(e.html("Temporary content for dialog").dialog())
        w(e.dialog("close"))
        w.flush()
        self.dialog = w
        big_assembly = widgets.HBox(children=[tsne_drawing.target.widget, big_slide_drawing.target.widget])
        self.assembly = widgets.VBox(children=[big_assembly, slide_drawing.target.widget, self.dialog])
        self.tsne_controller = TSNE(data, tsne_drawing)
        self.slide_controller_rows = slide_controllers
        self.chosen_slide = slides[0]
        self.chosen_indicator = None
        self.events_enabled = False
        self.reset_bookkeeping()

    def reset_bookkeeping(self):
        self.lasso_points = []
        self.highlight_slide = None
        self.selected_data = None

    def enable_events(self):
        self.slide_drawing.enable_events("click", self.slide_drawing_click)
        self.events_enabled = True
        tsne_events = "click mousedown mouseup mousemove mouseout"
        #tsne_events = "click"
        #print "enabling tsne events", tsne_events
        self.tsne_drawing.enable_events(tsne_events, self.tsne_event_callback)
        self.tsne_drawing.flush()
        self.slide_drawing.flush()

    def slide_drawing_click(self, info):
        #import pprint
        #pprint.pprint(info)
        (x, y) = info["point"]
        #print "slide click at", (x,y), info.get("type")
        for controller in self.slide_controllers:
            (sx, sy) = controller.xy
            if x > sx and x < sx + self.data.slide_width and y > sy and y < sy + self.data.slide_height:
                self.chosen_slide = controller.name
                self.highlight_slide = controller.name
                self.draw_later()
                return

    def draw_later(self):
        "Draw after delay in another thread to prevent callback handshake delays."
        #print "delayed draw"
        t = Timer(0.1, self.draw)
        t.start()

    def tsne_event_callback(self, info):
        try:
            point = info.get("point")
            (px, py) = point
            ty = info.get("type")
            #print "TSNE event at", (point, ty)
            pageX = info["pageX"]
            pageY = info["pageY"]
            data = self.data
            lasso = self.lasso_points
            if ty == "mouseout":
                self.lasso_points = []
            elif ty == "mousedown":
                #print "lasso_start", point
                self.lasso_points = [point]
            elif ty == "mousemove":
                addit = True
                if len(lasso) > 0:
                    (lx, ly) = lasso[-1]
                    if lx != px and ly != py:
                        self.tsne_drawing.line(None, px, py, lx, ly, "red")
                        self.tsne_drawing.flush()
                    else:
                        addit = False
                    if addit:
                        lasso.append(point)
            elif ty == "mouseup":
                #print "lasso complete"
                if len(lasso) > 2:
                    (px, py) = lasso[0]
                    (lx, ly) = lasso[-1]
                    if lx != px and ly != py:
                        self.tsne_drawing.line(None, px, py, lx, ly, "red")
                        self.tsne_drawing.flush()
                    names = data.data_names()
                    positions = data.data_tsne_positions()
                    inside = contourist.lasso.inside_lasso(positions, lasso)
                    #print "inside", inside
                    selected = set(names[i] for i in inside)
                    #print "selected", selected
                    self.selected_data = selected
                    self.draw_later()
                #print lasso
                self.lasso_points = []
        except Exception as e:
            #print "exception in tsne event", e
            raise

    def show(self):
        display(self.assembly)
        self.draw()

    def draw(self):
        self.slide_drawing.empty()
        self.tsne_drawing.empty()
        self.slide_drawing.flush()
        self.tsne_drawing.flush()
        #return
        self.tsne_controller.draw_on(0, 0, slide=self.highlight_slide, data=self.data, selected=self.selected_data)
        margin = self.slide_margin
        #offset = margin + SLIDE_SIZE
        xoffset = margin + self.data.slide_width
        yoffset = margin + self.data.slide_height
        #dy = margin
        dy = yoffset * (len(self.controller_rows) - 1) + margin
        for row in self.controller_rows:
            dx = margin
            for controller in row:
                controller.draw_on(dx, dy, self.data, indicator=self.chosen_indicator, selected=self.selected_data)
                controller.xy = (dx, dy)
                dx += xoffset
                #break
            dy -= yoffset
        slide_name = self.chosen_slide
        controller = self.slide_name_to_controller[slide_name]
        self.big_slide_drawing.empty()
        controller.draw_on(0,0,self.data,self.big_slide_drawing, selected=self.selected_data)
        if not self.events_enabled:
            self.enable_events()
        self.reset_bookkeeping()
