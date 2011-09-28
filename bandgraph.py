
import os,sys

from scalespace import scalespace
from numpy import *

from enthought.enable.example_support import DemoFrame, demo_main

# Enthought imports
from enthought.enable.api import *
from enthought.traits.api import *
from enthought.traits.ui.api import Item, Group, View
from enthought.util.resource import find_resource

# Chaco imports
from enthought.chaco.api import *
from enthought.chaco.tools.api import RangeSelection

time_ticks = [
    ('us', 1e-6),
    ('ms', 1e-3),
    ('sec', 1),
    ('min', 60),
    ('hour', 3600),
    ('day', 86400),
    ('month', 86400.*365./12.),
    ('year', 86400*365),
    ('decade', 86400*365*10),
    ('century', 86400*365*100)]
class TimeIntervalTickGenerator(AbstractTickGenerator):
    def _mkticks(self, low, high):
        l = []
        t = []
        for name,v in time_ticks:
            if low <= v <= high:
                l.append(name)
                t.append(v)
        return l,t

    def get_ticks(self, data_low, data_high, bounds_low, bounds_high, interval, use_endpoints = False, scale = 'log'):
        if scale != 'log': raise Exception("Only log scales allowed")
        return self._mkticks(bounds_low, bounds_high)[1]

    def get_ticks_and_labels(self, data_low, data_high, bounds_low, bounds_high,
                             orientation = "h"):
        return self._mkticks(bounds_low, bounds_high)

        

class PickerOverlay(AbstractOverlay):
    source = Instance(BaseXYPlot)
    destination = Instance(Component)
    
    border_color = ColorTrait((0, 0, 0.7, 1))
    border_width = Int(1)
    fill_color = ColorTrait("dodgerblue")
    alpha = Float(0.3)
    
    def calculate_points(self, component):
        """
        Calculate the overlay polygon based on the selection and the location
        of the source and destination plots.
        """
        # find selection range on source plot
        start, end = self._get_selection_screencoords()
        if start > end:
            start, end = end, start
        
        if self.source.orientation == 'h':
            edge_end = self.source.y
            edge_start = self.source.y2
            coord = lambda x,y: array([x,y])
        else:
            edge_end = self.source.x
            edge_start = self.source.x2
            coord = lambda y,x: array([x,y])
        
        left_top = coord(start, edge_start)
        left_mid = coord(start, edge_end)
        right_top = coord(end, edge_start)
        right_mid = coord(end, edge_end)
        
        polygon = array((left_top, left_mid, 
                         right_mid, right_top))
        left_line = array((left_top, left_mid))
        right_line = array((right_top,right_mid))
        
        return left_line, right_line, polygon
    
    def overlay(self, component, gc, view_bounds=None, mode="normal"):
        """
        Draws this overlay onto 'component', rendering onto 'gc'.
        """

        tmp = self._get_selection_screencoords()
        if tmp is None:
            return
        
        left_line, right_line, polygon = self.calculate_points(component)
       
        gc.save_state()
        try:
            gc.translate_ctm(*component.position)
            gc.set_alpha(self.alpha)
            gc.set_fill_color(self.fill_color_)
            gc.set_line_width(self.border_width)
            gc.set_stroke_color(self.border_color_)
            gc.begin_path()
            gc.lines(polygon)
            gc.fill_path()
            
            gc.begin_path()
            gc.lines(left_line)
            gc.lines(right_line)
            gc.stroke_path()
        finally:
            gc.restore_state()
        return
    
    def _get_selection_screencoords(self):
        """
        Returns a tuple of (x1, x2) screen space coordinates of the start
        and end selection points.  If there is no current selection, then
        returns None.
        """
        selection = self.source.index.metadata["selections"]
        if selection is not None and len(selection) == 2:
            mapper = self.source.index_mapper
            return mapper.map_screen(array(selection))
        else:
            return None
    
    #------------------------------------------------------------------------
    # Trait event handlers
    #------------------------------------------------------------------------
    
    def _source_changed(self, old, new):
        if old is not None and old.controller is not None:
            old.controller.on_trait_change(self._selection_update_handler, "selection",
                                           remove=True)
        if new is not None and new.controller is not None:
            new.controller.on_trait_change(self._selection_update_handler, "selection")
        return

    def _selection_update_handler(self, value):
        if value is not None and self.destination is not None:
            self._selection_changed(amin(value), amax(value))
        self.source.request_redraw()
        self.destination.request_redraw()
    def _selection_changed(self, start, end):
        pass


class TimePickerOverlay(PickerOverlay):
    def _selection_changed(self, start, end):
        r = self.destination.index_mapper.range
        r.low = start
        r.high = end

class BandPickerOverlay(PickerOverlay):
    sample_dt = Float
    data = Any
#    raw_data = Any

    def _selection_changed(self, start, end):
        print start, end, end/self.sample_dt
        self.destination.value.set_data(self.data.get(end / self.sample_dt))

#        f1, f2 = self.sample_dt/start, self.sample_dt/end
#        if f2 < f1: f2,f1 = f1,f2
#        if f2 - f1 < 0.1: return
#        print f1, f2
#        b,a = scipy.signal.butter(*scipy.signal.buttord([f1,f2],
#                                                        [f1/2,f2*2],
#                                                        1, 10),
#                                   btype="bandpass")
#        b,a = scipy.signal.butter(*scipy.signal.buttord(f2,
#                                                        max(f2*2,1),
#                                                        1, 10),
#                                   btype="low")
#
#        newdata = scipy.signal.filtfilt(b,a, self.raw_data)
#        self.data[:] = newdata
#        self.destination.value.set_data(newdata)
#        print self.data[0:100]
        


def create_gridded_line_plot(x, y, orientation="h", color="red",
                             dash="solid", value_mapper_class=LinearMapper,
                             padding=30, **kwargs):

    assert len(x) == len(y)
    
    # If you know it is monotonically increasing, sort_order can
    # be set to 'ascending'
    index = ArrayDataSource(x,sort_order='none')
    value = ArrayDataSource(y, sort_order="none")
    
    index_range = DataRange1D(tight_bounds = True)
    index_range.add(index)
    index_mapper = LinearMapper(range=index_range)

    value_range = DataRange1D(tight_bounds = True)
    value_range.add(value)
    value_mapper = value_mapper_class(range=value_range)
    
    plot = LinePlot(index=index, value=value,
                    index_mapper = index_mapper,
                    value_mapper = value_mapper,
                    orientation = orientation,
                    color = color,
                    line_style = dash,
                    padding = [10, 10, 10, 10],   # left, right, top, bottom
                    border_visible = True,
                    border_width = 1,
                    bgcolor = "white",
                    use_backbuffer = True,
                    backbuffer_padding = False,
                    unified_draw = True,
                    draw_layer = "plot",
                    overlay_border = True,
                    **kwargs)

    vertical_grid = PlotGrid(component = plot,
                             mapper=index_mapper,
                             orientation="vertical",
                             line_color="gray",
                             line_style='dot',
                             use_draw_order = True)

    horizontal_grid = PlotGrid(component = plot,
                               mapper=value_mapper,
                               orientation="horizontal",
                               line_color="gray",
                               line_style='dot',
                               use_draw_order = True)
    
    vertical_axis = PlotAxis(orientation="left",
                             mapper=plot.value_mapper,
                             use_draw_order = True)
    
    horizontal_axis = PlotAxis(orientation="bottom",
                               title='Time (s)',
                               mapper=plot.index_mapper,
                               use_draw_order = True)
    
    plot.underlays.append(vertical_grid)
    plot.underlays.append(horizontal_grid)
    
    # Have to add axes to overlays because we are backbuffering the main plot,
    # and only overlays get to render in addition to the backbuffer.
    plot.overlays.append(vertical_axis)
    plot.overlays.append(horizontal_axis)
    return plot


def create_band_picker_plot(x, y, orientation="h", color="red",
                             dash="solid",
                             padding=30, **kwargs):

    assert len(x) == len(y)
    
    # If you know it is monotonically increasing, sort_order can
    # be set to 'ascending'
    index = ArrayDataSource(x,sort_order='none')
    value = ArrayDataSource(y, sort_order="none")
    
    index_range = DataRange1D(tight_bounds = True)
    index_range.add(index)
    index_mapper = LogMapper(range=index_range)

    value_range = DataRange1D(tight_bounds = True)
    value_range.add(value)
    value_mapper = LinearMapper(range=value_range)
    
    plot = LinePlot(index=index, value=value,
                    index_mapper = index_mapper,
                    value_mapper = value_mapper,
                    orientation = orientation,
                    color = color,
                    line_style = dash,
                    padding = [10, 10, 10, 10],   # left, right, top, bottom
                    border_visible = True,
                    border_width = 1,
                    bgcolor = "white",
                    use_backbuffer = True,
                    backbuffer_padding = False,
                    unified_draw = True,
                    draw_layer = "plot",
                    overlay_border = True,
                    **kwargs)


    value_grid = PlotGrid(component = plot,
                               mapper=value_mapper,
                               orientation="vertical",
                               line_color="gray",
                               line_style='dot',
                               use_draw_order = True)


    value_axis = PlotAxis(orientation="bottom",
                             mapper=plot.value_mapper,
                             use_draw_order = True)


    tgen = TimeIntervalTickGenerator()


    index_grid = PlotGrid(component = plot,
                             mapper=index_mapper,
                             orientation="horizontal",
                             line_color="gray",
                             line_style='dot',
                             tick_generator = tgen,
                             use_draw_order = True)
    
    index_axis = PlotAxis(orientation="left",
                               title='Time (s)',
                               mapper=plot.index_mapper,
                               tick_generator = tgen,
                               use_draw_order = True)
    
    plot.underlays.append(value_grid)
    plot.underlays.append(index_grid)
    
    # Have to add axes to overlays because we are backbuffering the main plot,
    # and only overlays get to render in addition to the backbuffer.
    plot.overlays.append(value_axis)
    plot.overlays.append(index_axis)
    return plot



import csv, time, datetime
print "Loading data...",
currency_data=array([[time.mktime(datetime.datetime(int(d[1][:4]),int(d[1][4:6]),int(d[1][6:8]),*map(int,d[2].split(":"))).timetuple()), float(d[-1])] for d in csv.reader(open("EURUSD_hour.csv")) if d[0] != '<TICKER>']).T
print "done"
currency_data_dt = 3600

def create_zoomed_plot():
    x = currency_data[0]
    y = array(currency_data[1])

    sz_main = [600,400]
    sz_pick = 70

    data_plot = create_gridded_line_plot(x,y,resizable="hv")
    time_picker_plot = create_gridded_line_plot(x,y, resizable="h", height=sz_pick)
    banddata_x = [3600, 86400, 86400*30]
    banddata_y = [20, 40, 30]
    band_picker_plot = create_band_picker_plot(banddata_x, banddata_y, resizable="v", width=sz_pick, orientation="v", value_mapper_class=LogMapper)

#    time_picker_plot.height = sz_pick
#    band_picker_plot.width = sz_pick

    outer_container = GridPlotContainer(padding=30,
                                        fill_padding=True,
                                        spacing=(10,10),
                                        bgcolor="lightgray",
                                        shape=(2,2),
                                        use_backbuffer=True)
#    outer_container = VPlotContainer(padding=30,
#                                     fill_padding=True,
#                                     spacing=50,
#                                     stack_order='top_to_bottom',
#                                     bgcolor='lightgray',
#                                     use_backbuffer=True)

    outer_container.add(data_plot)
    outer_container.add(band_picker_plot)
    outer_container.add(time_picker_plot)

    
    time_picker_plot.controller = RangeSelection(time_picker_plot, left_button_selects=True)
    time_picker_overlay = TimePickerOverlay(source=time_picker_plot, destination=data_plot)
    outer_container.overlays.append(time_picker_overlay)

    band_picker_plot.controller = RangeSelection(band_picker_plot, left_button_selects=True)
    band_picker_overlay = BandPickerOverlay(source=band_picker_plot, destination=data_plot, sample_dt=currency_data_dt, data=scalespace(currency_data[1],2))
    outer_container.overlays.append(band_picker_overlay)

    

#    data_plot.bounds = sz_main
#    time_picker_plot.bounds = [sz_main[0], sz_pick]
#    band_picker_plot.bounds = [sz_pick, sz_main[1]]
    


    return outer_container

size = (800, 600)
title = "EUR/USD"


#===============================================================================
# Stand-alone frame to display the plot.
#===============================================================================
class PlotFrame(DemoFrame):

    def _create_window(self):
        # Return a window containing our plots
        return Window(self, -1, component=create_zoomed_plot())
    
if __name__ == "__main__":
    demo_main(PlotFrame, size=size, title=title)

#--EOF---
