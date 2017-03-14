
import numpy as np 
from numpy.linalg import norm
from jp_svg_canvas import canvas, cartesian_svg
from jp_gene_viz import proxy_html5_canvas
from jp_gene_viz import js_proxy
import contourist.lasso
import os
from PIL import Image

canvas.load_javascript_support()

import ipywidgets as widgets
from IPython.display import display

Vent_Med_White = 'Vent_Med_White'
Vent_Horn = 'Vent_Horn'
Vent_Lat_White = 'Vent_Lat_White'
Med_Grey = 'Med_Grey'
Dors_Horn = 'Dors_Horn'
Dors_Edge = 'Dors_Edge'
Med_Lat_White = 'Med_Lat_White'
Vent_Edge = 'Vent_Edge'
Med_Edge = 'Med_Edge'
Dors_Med_White = 'Dors_Med_White'
Cent_Can = 'Cent_Can'
Lat_Edge = 'Lat_Edge'
Undefined = 'Undefined'

LABELS = [
    Vent_Med_White,
    Vent_Horn,
    Vent_Lat_White,
    Med_Grey,
    Dors_Horn,
    Dors_Edge,
    Med_Lat_White,
    Vent_Edge,
    Med_Edge,
    Dors_Med_White,
    Cent_Can,
    Lat_Edge,
    Undefined,
]

COLORS = {
    Vent_Med_White: '#7722dd',
    Vent_Horn: '#ee45ba',
    Vent_Lat_White: '#656897',
    Med_Grey: '#dc8b74',
    Dors_Horn: '#53ae51',
    Dors_Edge: '#cad12e',
    Med_Lat_White: '#41f40b',
    Vent_Edge: '#b916e8',
    Med_Edge: '#3039c5',
    Dors_Med_White: '#a75ca2',
    Cent_Can: '#1e7f7f',
    Lat_Edge: '#95a25c',
    Undefined: '#0cc539',
}

def add_image(svg_canvas, image_href, x=0, y=0, width=None, height=None, name="background",
        opacity=1.0):
    tag = "image"
    if width is None:
        width = svg_canvas.svg_width
    if height is None:
        height = svg_canvas.svg_height
    atts = {
        "width": width,
        "height": height,
        "href": image_href,
        "opacity": opacity
    }
    svg_canvas.add_element(name, tag, atts)
    svg_canvas.send_commands()
    return (width, height)

class SpotAnnotator(object):

    selected_label = None

    def __init__(self, image_href, tsv_data_path, annotation_path=None, 
        mins=(2,2), maxes=(32,34), radius=0.3):
        self.radius = radius
        self.mins = np.array(mins, dtype=np.float)
        self.maxes = np.array(maxes, dtype=np.float)
        self.offset = self.maxes - self.mins
        assert os.path.exists(tsv_data_path)
        assert os.path.exists(image_href)
        self.image_href = image_href
        self.tsv_data_path = tsv_data_path
        if annotation_path is None:
            annotation_path = tsv_data_path + ".annotations.tsv"
        self.annotation_path = annotation_path
        self.get_coordinates()
        # get image dimensions
        im = Image.open(image_href)
        self.image_size = np.array(im.size)

    def show(self):
        ls = self.label_selection = self.label_selector()
        sd = self.spot_display = self.display_spot_canvas()
        #detail = self.detail_display = self.display_spot_canvas(width=2*self.radius, height=2*self.radius)
        im_detail = self.image_detail = js_proxy.ProxyWidget()
        #self.assembly = self.label_selection # testing only
        info = self.info_area = widgets.Textarea(description="status")
        filename = self.file_name_text = widgets.Text(value=self.annotation_path)
        savebutton = self.save_button = widgets.Button(description="save")
        savebutton.onclick = self.save_click
        file_save = widgets.HBox(children=[filename, savebutton])
        work_area = widgets.HBox(children=[sd.target, ls, im_detail])
        self.assembly = widgets.VBox(children=[work_area, file_save, info])
        #self.assembly = widgets.HBox(children=[ls])
        display(self.assembly)
        self.configure_label_selector()
        self.configure_image_detail()
        self.label_selection.flush()
        self.draw_spots()
        #(width, height) = self.img_wh
        #(width, height) = detail.project(width, height)
        #add_image(detail.target, self.image_href, opacity=1.0, width=abs(width), height=abs(height))
        #self.label_handler(Undefined, None)

    def save_click(self, *args):
        filename = self.file_name_text.value
        self.info_area.value = "saving to " + filename
        self.dump_annotations(filename)
        self.info_area.value = "saved to " + filename

    def circle_name(self, index):
        return "circle_" + str(index)

    def draw_spots(self, image=True, dots=True):
        drawing = self.spot_display
        coords = self.coordinates
        labels = self.labels
        radius = self.radius
        #add_image(drawing.target, self.image_href)
        names, xs, ys, radii, clrs = [], [], [], [], []
        if dots:
            for (index, (x,y)) in enumerate(coords):
                (x,y) = self.adjust(x,y)
                label = labels[index]
                clr = COLORS[label]
                name = self.circle_name(index)
                #drawing.circle(name, x, y, radius, clr)
                names.append(name)
                xs.append(x)
                ys.append(y)
                radii.append(radius)
                clrs.append(clr)
            drawing.circles(names, xs, ys, radii, clr)
            #drawing.flush()
        if image:
            import time
            #time.sleep(0.2)
            self.img_wh = add_image(drawing.target, self.image_href, opacity=0.5)
            #drawing.flush()
        drawing.enable_events("click", self.spot_callback)
        self.info_area.value = "spots drawn " + repr(len(coords))
        print "drew", len(coords), "spots"
        drawing.flush()

    def spot_callback(self, info):
        [px, py] = info["point"]
        display = self.spot_display
        display.delete("click_circle")
        atts = {"opacity": 0.3}
        display.circle("click_circle", px, py, self.radius, "black", other_attributes=atts)
        chosen_index = chosen_xy = None
        radius = self.radius
        coords = self.coordinates
        for (index, (x,y)) in enumerate(coords):
            (x,y) = self.adjust(x,y)
            l1dist = max([abs(x-px), abs(y-py)])
            if l1dist < 2 * radius:
                chosen_index = index
                chosen_xy = (x, y)
                #print "chose", index, (x,y), (px, py)
        if chosen_index is not None:
            self.set_label(chosen_index)
            #self.set_detail(chosen_xy)
            detail = self.focus_image_detail(chosen_xy)
            self.info_area.value = "clicked %s: %s"  % (chosen_index, chosen_xy)
        else:
            chosen_xy = (px, py)
            detail = self.focus_image_detail((px,py))
            self.info_area.value = "no match for click " + repr((px, py))
        self.info_area.value = "clicked %s: %s %s"  % (chosen_index, chosen_xy, detail)
        self.highlight_detail(detail)

    def highlight_detail(self, detail):
        ((x,y), (w,h)) = detail
        (x,y) = self.radjust(x,y)
        display = self.spot_display
        display.delete("reference")
        atts = {"opacity": 0.1}
        y = y - h
        print "highlight", x, y, w, h
        display.rect("reference", x, y, w, h, "red", other_attributes=atts)
        display.flush()

    """def set_detail(self, xy):
        (x, y) = xy
        #detail = self.detail_display
        drawing = self.spot_display
        llx = x - self.radius
        lly = y - self.radius
        urx = llx + 2 * self.radius
        ury = lly + 2 * self.radius
        (llx_s, lly_s) = drawing.project(llx, lly)
        (urx_s, ury_s) = drawing.project(urx, ury)
        self.focus_image_detail(x, y)
        #detail.target.set_view_box(min(llx_s, urx_s), min(lly_s, ury_s), abs(llx_s-urx_s), abs(lly_s-ury_s))
        #print "detail view box", detail.target.viewBox"""
    
    def set_label(self, index):
        name = self.circle_name(index)
        label = self.selected_label
        if label is None:
            self.info_area.value = "no selected label " + repr(index)
            return
        self.labels[index] = label
        color = COLORS[label]
        self.spot_display.change(name, fill=color)
        self.info_area.value = "selected %s for %s" % (color, index)
        #print "set_label", index, name, label, color

    def adjust(self, x, y):
        (minx, miny) = self.mins
        (maxx, maxy) = self.maxes
        yoffset = y - miny
        newy = maxy - yoffset
        return (x, newy)

    radjust = adjust

    def display_spot_canvas(self, width=None, height=None):
        (minx, miny) = self.mins
        (maxx, maxy) = self.maxes
        if width is not None and height is not None:
            maxx = minx + width
            maxy = miny + height
        result = cartesian_svg.doodle(minx, miny, maxx, maxy, margin=0)
        result.buffered = True
        #add_image(result.target, self.image_href)
        return result
    
    def label_selector(self):
        w = js_proxy.ProxyWidget()
        return w
    
    def configure_label_selector(self):
        w = self.label_selection
        elt = w.element()
        jQuery = w.window().jQuery
        w(elt._set("top_div", jQuery("<div/>").width("100px").html("LABELS")))
        w(elt.append(elt.top_div))
        for label in LABELS:
            callback = w.callback(self.label_handler, data=label)
            clr = COLORS[label]
            w(elt._set(label,
                jQuery("<div/>")
                .html(label)
                .click(callback)
                .css('cursor', 'pointer')
                .css('color', clr)
                .appendTo(elt.top_div)))
        return w

    detail_size = 500  # TEMP!
    detail_side = 500

    def configure_image_detail(self):
        side = self.detail_side
        size_px = str(side) + "px"
        w = self.image_detail
        elt = w.element()
        jQuery = w.window().jQuery
        w(elt._set("canvas", jQuery("<canvas/>").attr({'width':side,'height':side}).width(size_px).height(size_px)))
        w(elt.append(elt.canvas))
        w(elt._set("ctx", elt.canvas[0].getContext("2d")))
        w(elt.ctx._set("fillStyle", "#ffffdd"))
        w(elt.ctx.fillRect(0, 0, 500, 500))
        w(elt._set("detail_image", jQuery("<img/>", {"src": self.image_href})))
        w(elt.detail_image.hide())
        w(elt.append(elt.detail_image))
        w.flush()

    def focus_image_detail(self, xy):
        # flip y 
        #print "input xy", xy
        offset_xy = np.array(xy) - self.mins
        #print "offset by mins", offset_xy, "extrema", self.mins, self.maxes, self.offset
        offset_xy[1] = self.offset[1] - offset_xy[1]# + self.radius
        #print "flipped y", offset_xy
        canvas_to_image = self.image_size / self.offset
        image_to_canvas = 1.0 / canvas_to_image
        #print "convert to canvas", image_to_canvas, "to image", canvas_to_image, "image size", self.image_size
        size = self.detail_size
        adjustment = np.array([0.5,0.5])
        size_shift = offset_xy - adjustment * size * image_to_canvas
        #print "shifted", size_shift, image_to_canvas
        size_shift = np.max([np.zeros((2,)), np.min([self.offset, size_shift], axis=0)], axis=0)
        #print "shifted", size_shift
        image_offset_xy = canvas_to_image * size_shift
        w = self.image_detail
        elt = w.element()
        w(elt.ctx.clearRect(0, 0, size, size))
        #sx = sy = self.image_width / 4 # temp
        (sx, sy) = image_offset_xy
        #print "image offset", sx, sy
        sw = sh = size
        dx = dy = 0
        dw = dh = self.detail_side
        w(w.function(["element"], "debugger;")(elt))
        #print "detail_image", sx, sy, sw, sh, dx, dy, dw, dh
        w(elt.ctx.drawImage(elt.detail_image[0], sx, sy, sw, sh, dx, dy, dw, dh))
        cx = cy = size * 0.5
        radius = self.radius * canvas_to_image[0]
        w(elt.ctx.beginPath())
        w(elt.ctx.arc(cx, cy, radius, 0, 6.28))
        w(elt.ctx.stroke())
        w.flush()
        corner = size_shift + self.mins
        dimensions = size * image_to_canvas
        #print "corner", corner, "dimensions", dimensions
        return (corner, dimensions)

    def label_handler(self, label, args):
        #print "label handler got", label
        w = self.label_selection
        elt = w.element()
        old_selected_label = self.selected_label
        self.selected_label = label
        if old_selected_label:
            w(elt[old_selected_label]
                .css("background-color", "white")
                .css("color", COLORS[old_selected_label]))
        w(elt[label]
            .css("background-color", COLORS[label])
            .css("color", "white"))
        w.flush()

    def get_coordinates(self):
        path = self.tsv_data_path
        (folder, filename) = os.path.split(path)
        (prefix, extension) = filename.split(".")
        (slide, array_name) = prefix.split("_")
        self.slide = slide
        self.array_name = array_name
        f = open(path)
        first_line = f.readline()
        coordinates = []
        labels = []
        for coords_str in first_line.split():
            xy = map(float, coords_str.strip().split("_"))
            assert len(xy) == 2
            coordinates.append(xy)
            labels.append(Undefined)
        self.coordinates = coordinates
        self.labels = labels

    def dump_annotations(self, to_path=None):
        if to_path is None:
            to_path = self.annotation_path
        f = open(to_path, "w")
        headers = ["Slide", "Array", "xPos", "yPos"] + LABELS
        f.write(("\t".join(headers) + "\n"))
        slide_array = [self.slide, self.array_name]
        for (i, xy) in enumerate(self.coordinates):
            label = self.labels[i]
            label_index = LABELS.index(label)
            indicators = [0] * len(LABELS)
            indicators[label_index] = 1
            row = slide_array + list(xy) + indicators
            srow = map(str, row)
            f.write(("\t".join(srow))+"\n")
        f.close()
