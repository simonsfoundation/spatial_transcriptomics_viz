
import numpy as np 
from numpy.linalg import norm
from jp_svg_canvas import canvas, cartesian_svg
from jp_gene_viz import proxy_html5_canvas
from jp_gene_viz import js_proxy
import contourist.lasso
import os
from PIL import Image
import shutil
from threading import Timer

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

NEARBY = "NEARBY"

LABELS = sorted([
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
])

# old colors
COLORS = {
    Vent_Med_White: '#aaaa66',
    Vent_Horn: '#ee45ba',
    Vent_Lat_White: '#66aaaa',
    Med_Grey: '#997744',
    Dors_Horn: '#447799',
    Dors_Edge: '#774499',
    Med_Lat_White: '#aa66aa',
    Vent_Edge: '#994477',
    Med_Edge: '#997744',
    Dors_Med_White: '#aaddee',
    Cent_Can: '#6666ff',
    Lat_Edge: '#ff6666',
    Undefined: '#dddddd',
}

# new colors
COLORS = {
    Vent_Med_White: '#66aaaa',
    Vent_Horn: '#ffaa00',
    Vent_Lat_White: '#99aaaa',
    Med_Grey: '#ff00aa',
    Dors_Horn: '#99aa00',
    Dors_Edge: '#66aa00',
    Med_Lat_White: '#ffaaaa',
    Vent_Edge: '#ff0000',
    Med_Edge: '#990000',
    Dors_Med_White: '#9900aa',
    Cent_Can: '#6600aa',
    Lat_Edge: '#660000',
    Undefined: '#dddddd',
}


ARROWS = [LEFT, UP, RIGHT, DOWN] = "LEFT UP RIGHT DOWN".split()

NAME_TO_ARROW_KEY_NUMBER = {
    LEFT: 37,
    UP: 38,
    RIGHT: 39,
    DOWN: 40,
}
NAME_TO_KEY_NUMBER = {}
LABEL_KEY_0 = ord('a')
for (index, label) in enumerate(LABELS):
    NAME_TO_KEY_NUMBER[label] = LABEL_KEY_0 + index
NAME_TO_KEY_NUMBER[NEARBY] = ord('x')

KEY_HANDLER_TEMPLATE = """function(event) {
        switch (event.which) {
            %(cases)s

            default: return; // allow defaults for other keys.
        }
        event.preventDefault();  // disallow defaults for handled keys
    };
"""

def key_handler(dictionary):
    cases = key_binding_cases(dictionary)
    return KEY_HANDLER_TEMPLATE % {"cases": cases}

KEY_BINDING_TEMPLATE = """
    // arguments: element, proxy_callback
    debugger;
    var keypress_handler = %(keypress_handler)s
    var keydown_handler = %(keydown_handler)s
    element.keydown(keydown_handler);
    element.keypress(keypress_handler);
"""

def key_binding_cases(dictionary):
    cases = []
    for (name, number) in sorted(dictionary.items()):
        case = "case %s: proxy_callback(%s); break;" % (number, repr(str(name)))
        cases.append(case)
    return "\n            ".join(cases)

def add_key_bindings_to_parent(w, callback, identifier=None):
    elt = w.element()
    parent = elt.parent();
    proxy_callback = w.callback(callback, identifier)
    D = {}
    D["keypress_handler"] = key_handler(NAME_TO_KEY_NUMBER)
    D["keydown_handler"] = key_handler(NAME_TO_ARROW_KEY_NUMBER)
    key_binding_fn_body = KEY_BINDING_TEMPLATE % D
    anon_fn = w.function(["element", "proxy_callback"], key_binding_fn_body)
    #w(anon_fn(elt, proxy_callback))
    w(anon_fn(parent, proxy_callback))
    w(parent.attr('tabindex', 0))

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
    selected_index = None

    def __init__(self, tsv_data_path, image_path=None, image_href=None, annotation_path=None, spots_path=None,
        mins=(1,1), maxes=None, radius=0.3, scratchdir="image_copies"):
        #self.radius = radius
        #self.mins = np.array(mins, dtype=np.float)
        #self.maxes = np.array(maxes, dtype=np.float)
        #self.offset = self.maxes - self.mins
        assert os.path.exists(tsv_data_path)
        spots_href = None
        if image_href is not None:
            assert image_path is None, "Only provide one of image_path or image_href, please. " + repr((image_path, image_href))
            image_path = image_href
        else:
            assert image_path is not None, "Image path or href is required."
            image_href = self.copy_image(image_path, scratchdir)
            if spots_path is not None:
                spots_href = self.copy_image(spots_path, scratchdir)
        assert os.path.exists(image_href), "image not found " + repr(image_href)
        self.image_href = image_href
        self.spots_href = spots_href
        self.tsv_data_path = tsv_data_path
        if annotation_path is None:
            annotation_path = tsv_data_path + ".annotations.tsv"
        self.annotation_path = annotation_path
        self.get_coordinates()
        # get image dimensions
        im = Image.open(image_href)
        self.image_size = np.array(im.size)
        self.mins = np.array(mins, dtype=np.float)
        if maxes is None:
            spot_to_pixels = self.image_size[0] * 194.0 / 6200.0
            self.maxes = self.mins + (self.image_size / spot_to_pixels)
        else:
            self.maxes = np.array(maxes, dtype=np.float)
        self.offset = self.maxes - self.mins
        self.radius = radius

    def copy_image(self, image_path, scratchdir):
        assert os.path.isfile(image_path), "No such file " + repr(image_path)
        # copy the image to scratch area
        if not os.path.isdir(scratchdir):
            os.mkdir(scratchdir)
        image_filename = os.path.split(image_path)[-1]
        image_href = scratchdir + "/" + image_filename
        shutil.copyfile(image_path, image_href)
        print(image_path + " copied to " + image_href)
        return image_href

    def show(self):
        ls = self.label_selection = self.label_selector()
        sd = self.spot_display = self.display_spot_canvas()
        #detail = self.detail_display = self.display_spot_canvas(width=2*self.radius, height=2*self.radius)
        im_detail = self.image_detail = js_proxy.ProxyWidget()
        #self.assembly = self.label_selection # testing only
        info = self.info_area = widgets.Textarea(description="status")
        filename = self.file_name_text = widgets.Text(value=self.annotation_path)
        savebutton = self.save_button = widgets.Button(description="save")
        savebutton.on_click(self.save_click)
        restorebutton = self.restore_button = widgets.Button(description="restore")
        restorebutton.on_click(self.restore_click)
        nearbybutton = self.nearby_button = widgets.Button(description="x: nearby")
        nearbybutton.on_click(self.nearby_click)
        file_save_widgets = [filename, savebutton, restorebutton]
        if self.spots_href is not None:
            spots_checkbox = self.spots_checkbox = widgets.Checkbox(description="spots")
            file_save_widgets.append(spots_checkbox)
            spots_checkbox.on_trait_change(self.redraw, "value")
        file_save_widgets.append(nearbybutton)
        file_save = widgets.HBox(children=file_save_widgets)
        work_area = widgets.HBox(children=[sd.target, ls, im_detail])
        self.assembly = widgets.VBox(children=[work_area, file_save, info])
        #self.assembly = widgets.HBox(children=[ls])
        display(self.assembly)
        self.configure_label_selector()
        self.configure_image_detail()
        self.add_key_bindings()
        self.label_selection.flush()
        if os.path.exists(self.annotation_path):
            self.restore_click()
        else:
            self.draw_spots()
        #(width, height) = self.img_wh
        #(width, height) = detail.project(width, height)
        #add_image(detail.target, self.image_href, opacity=1.0, width=abs(width), height=abs(height))
        #self.label_handler(Undefined, None)

    def restore_click(self, *args):
        filename = self.file_name_text.value
        self.info_area.value = "restoring from " + filename
        result = self.restore_annotations(filename)
        if result is True:
            self.draw_spots()
            self.info_area.value = "restored from " + filename
        else:
            self.info_area.value = "failed to restore " + repr((filename, result))

    def redraw(self, *args):
        t = Timer(0.1, self.draw_spots)
        t.start()

    def save_click(self, *args):
        filename = self.file_name_text.value
        self.info_area.value = "saving to " + filename
        self.dump_annotations(filename)
        self.info_area.value = "saved to " + filename

    def nearby_click(self, *args):
        changed = []
        last_xy = self.last_xy
        if last_xy is None:
            self.info_area.value = "No position selected: cannot find nearby spots."
            return
        (last_x, last_y) = last_xy
        label = self.selected_label
        if label is None:
            self.info_area.value = "No label selected: cannot label nearby spots."
            return
        selected_index = self.selected_index
        radius = self.radius * 4
        coords = self.coordinates
        for (index, (x, y)) in enumerate(coords):
            (x,y) = self.adjust(x,y)
            distance = max(abs(x - last_x), abs(y-last_y))
            if distance <= radius and index != selected_index:
                self.set_label(index)
                changed.append(index)
        if len(changed) == 0:
            self.info_area.value = "No spots found in radius " + repr(radius)
        else:
            self.spot_display.flush()
            self.info_area.value = "Set %s to %s within radius %s." % (changed, label, radius)

    def circle_name(self, index):
        return "circle_" + str(index)

    def draw_spots(self, image=True, dots=True):
        drawing = self.spot_display
        drawing.empty()
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
            drawing.circles(names, xs, ys, radii, clrs)
            #drawing.flush()
        if image:
            import time
            #time.sleep(0.2)
            self.img_wh = add_image(drawing.target, self.image_href, opacity=0.5)
        if self.spots_href is not None and self.spots_checkbox.value:
            add_image(drawing.target, self.spots_href, opacity=0.5)
        drawing.enable_events("click mousemove", self.spot_callback)
        self.info_area.value = "spots drawn " + repr(len(coords))
        #print "drew", len(coords), "spots"
        drawing.flush()

    last_xy = None

    def spot_callback(self, info, direction=None, epsilon=0.2):
        [px, py] = info["point"]
        ty = info.get("type")
        #self.info_area.value = repr((px, py, ty))
        display = self.spot_display
        display.delete("click_circle")
        atts = {"opacity": 0.1}
        if direction is None:
            display.circle("click_circle", px, py, self.radius, "black", other_attributes=atts)
        chosen_index = chosen_xy = None
        radius = self.radius
        coords = self.coordinates
        last_x = last_y = None
        if self.last_xy is not None:
            #print "last_xy", self.last_xy
            (last_x, last_y) = self.last_xy
        min_dist = None
        for (index, (x,y)) in enumerate(coords):
            (x,y) = self.adjust(x,y)
            if last_x is not None and direction is not None:
                if direction == UP and last_y + epsilon >= y:
                    continue
                elif direction == DOWN and last_y - epsilon <= y:
                    continue
                elif direction == LEFT and last_x - epsilon <= x:
                    continue
                elif direction == RIGHT and last_x + epsilon >= x:
                    continue
            #l1dist = max([abs(x-px), abs(y-py)])
            l1dist = abs(x-px) + abs(y-py)
            if direction is None and l1dist > 2 * radius:
                continue
            if min_dist is None or min_dist > l1dist:
                min_dist = l1dist
                chosen_index = index
                chosen_xy = (x, y)
                #print "chose", index, (x,y), (px, py), (last_x, last_y), repr(direction)
        if ty == "mousemove":
            if chosen_index is None:
                self.info_area.value = "No match at " + repr((x,y))
            else:
                label = self.labels[chosen_index]
                self.info_area.value = "%s: %s at %s" % (chosen_index, label, chosen_xy)
            return
        label = None
        if chosen_index is not None:
            label = self.set_label(chosen_index)
            #self.set_detail(chosen_xy)
            detail = self.focus_image_detail(chosen_xy)
            self.last_xy = chosen_xy
            #self.info_area.value = "clicked %s: %s"  % (chosen_index, map(int, chosen_xy))
            (px, py) = chosen_xy
            if direction is not None:
                display.circle("click_circle", px, py, self.radius, "black", other_attributes=atts)
        else:
            chosen_xy = (px, py)
            self.last_xy = chosen_xy
            detail = self.focus_image_detail((px,py))
            self.info_area.value = "no match for click " + repr((px, py))
        #self.info_area.value = "clicked %s: %s %s"  % (chosen_index, chosen_xy, detail)
        self.selected_index = chosen_index
        self.highlight_detail(detail)

    def highlight_detail(self, detail):
        ((x,y), (w,h)) = detail
        (x,y) = self.radjust(x,y)
        display = self.spot_display
        display.delete("reference")
        atts = {"opacity": 0.1}
        y = y - h
        #print "highlight", x, y, w, h
        display.rect("reference", x, y, w, h, "red", other_attributes=atts)
        display.flush()
    
    def set_label(self, index):
        name = self.circle_name(index)
        label = self.selected_label
        if label is None:
            self.info_area.value = "no selected label " + repr(index)
            return None
        self.labels[index] = label
        color = COLORS[label]
        self.spot_display.change(name, fill=color)
        xy = map(int, self.coordinates[index])
        self.info_area.value = "selected %s for %s at %s." % (label, index, xy)
        return label
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

    def key_binding_callback(self, info, args):
        arg = args["0"]
        #print "key binding callback", arg
        if arg in LABELS:
            later(self.label_handler, [arg])
        elif arg in ARROWS:
            #print "later moving", arg
            #later(self.move_spot, [arg])
            self.move_spot(arg)
        elif arg == NEARBY:
            self.nearby_click()
        else:
            self.info_area.value = "UNKNOWN KEY: " + repr(args)

    def move_spot(self, direction):
        #print "moving spot", direction, self.last_xy
        if self.last_xy is not None:
            info = {"point": self.last_xy}
            self.spot_callback(info, direction)

    def add_key_bindings(self):
        w = self.label_selection
        cb = w.callback(self.key_binding_callback, None)
        #w(w.element()["top_div"].keypress(cb))
        #w.flush()
        add_key_bindings_to_parent(w, self.key_binding_callback)
    
    def configure_label_selector(self):
        w = self.label_selection
        elt = w.element()
        jQuery = w.window().jQuery
        w(elt._set("top_div", jQuery("<div/>").width("100px").html("LABELS").attr('tabindex', 0).focus()))
        w(elt.append(elt.top_div))
        ch_ord = LABEL_KEY_0
        for label in LABELS:
            callback = w.callback(self.label_handler, data=label)
            clr = COLORS[label]
            ch = chr(ch_ord)
            w(elt._set(label,
                jQuery("<div/>")
                .html("&nbsp;"+ ch + ":" + label)
                .click(callback)
                .css('cursor', 'pointer')
                .css('color', clr)
                .appendTo(elt.top_div)))
            ch_ord += 1
        return w

    detail_size = 500
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
        if self.spots_href:
            w(elt._set("spots_image", jQuery("<img/>", {"src": self.spots_href})))
            w(elt.spots_image.hide())
            w(elt.append(elt.spots_image))
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
        #w(w.function(["element"], "debugger;")(elt))
        #print "detail_image", sx, sy, sw, sh, dx, dy, dw, dh
        w(elt.ctx.drawImage(elt.detail_image[0], sx, sy, sw, sh, dx, dy, dw, dh))
        cx = cy = self.detail_side * 0.5
        radius = self.radius * canvas_to_image[0]
        w(elt.ctx.beginPath())
        w(elt.ctx.arc(cx, cy, radius, 0, 6.28))
        w(elt.ctx.stroke())
        if self.spots_href:
            w(elt.ctx._set("globalAlpha", 0.1))
            w(elt.ctx.drawImage(elt.spots_image[0], sx, sy, sw, sh, dx, dy, dw, dh))
            w(elt.ctx._set("globalAlpha", 1.0))
        w.flush()
        corner = size_shift + self.mins
        dimensions = size * image_to_canvas
        #print "corner", corner, "dimensions", dimensions
        return (corner, dimensions)

    def label_handler(self, label, args=None):
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
        if self.selected_index is not None:
            self.set_label(self.selected_index)

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

    fixed_headers = ["Slide", "Array", "xPos", "yPos"]

    def dump_annotations(self, to_path=None):
        if to_path is None:
            to_path = self.annotation_path
        f = open(to_path, "w")
        headers = self.fixed_headers + LABELS
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

    def restore_annotations(self, from_path=None, epsilon=1e-4):
        if not os.path.isfile(from_path):
            return "No such file " + repr(from_path)
        f = open(from_path)
        headers = f.readline().strip().split("\t")
        if headers[:4] != self.fixed_headers:
            return "Fixed headers do not match " + repr((from_path, headers, self.fixed_headers))
        label_order = headers[4:]
        if not set(label_order).issubset(set(LABELS)):
            return "Unknown labels " + repr((set(label_order) - set(LABELS), from_path))
        line = f.readline().strip()
        while line:
            items = line.split("\t")
            [slide, array, xps, yps] = items[:4]
            xp = float(xps)
            yp = float(yps)
            if slide == self.slide and array == self.array_name:
                index = None
                for (i, (x, y)) in enumerate(self.coordinates):
                    n = norm([(x - xp), (y - yp)])
                    if n < epsilon:
                        if index is not None:
                            return "ambiguous position" + repr((i, x, y, from_path))
                        index = i
                label = None
                try:
                    indicators = map(int, items[4:])
                except ValueError:
                    return "bad indicator line " + repr(items[4:])
                label = None
                for (lindex, indicator) in enumerate(indicators):
                    if indicator != 0:
                        if label is not None:
                            return "position has more than one label " + (i, x, y, label_order[lindex], label)
                        label = label_order[lindex]
                self.labels[index] = label
            line = f.readline().strip()
        return True  # success


def later(operation, args=()):
    t = Timer(0.1, operation, args)
    t.start()
