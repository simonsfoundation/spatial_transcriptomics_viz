#!/usr/bin/env python

import re
import glob
import pickle

import numpy
import pandas as pd

import bokeh.plotting
import bokeh.io
import bokeh.layouts
import bokeh.models

import colorcet

from PIL import Image

COVARIATES = {'Vent_Med_White': 'Ventral medial white','Vent_Horn': 'Ventral horn','Vent_Lat_White': 'Ventral lateral white','Med_Grey': 'Medial grey','Dors_Horn': 'Dorsal horn','Dors_Edge': 'Dorsal edge','Med_Lat_White': 'Medial lateral white','Vent_Edge': 'Ventral Edge','Dors_Med_White': 'Dorsal medial white','Cent_Can': 'Central canal','Lat_Edge': 'Lateral edge'}

def joy(category, data, scale=0.2):
  return list(zip([category]*len(data), scale*data))

def parse_P(P_df,count_files,metadata):
  P_pretty = []
  for count_file_idx in range(0,len(count_files)):
    filename = count_files[count_file_idx]

    coordinates = list(P_df[filename].columns)
 
    proportions = P_df[filename].as_matrix()
 
    tmp = metadata['Count file'] == filename.replace('Count_Tables/','')
    name,gt,tp,sex = metadata[tmp]['Name'].values[0],metadata[tmp]['Genotype'].values[0],metadata[tmp]['Timepoint'].values[0],metadata[tmp]['Sex'].values[0]
 
    if re.match(r".*_[0-9]+_stdata_aligned_counts_IDs.txt",filename):
      covariates = pd.read_table('Covariates/%s.tsv'%(re.sub(r"_[0-9]+_stdata_aligned_counts_IDs.txt",r"",filename.replace('Count_Tables/',''))),header=0,index_col=0)
    else:
      covariates = pd.read_table('Covariates/%s.tsv'%(filename.replace('_stdata_aligned_counts_IDs.txt','').replace('Count_Tables/','')),header=0,index_col=0)
 
    for idx,coordinate in enumerate(coordinates):
      P_pretty.append([name,coordinate,gt,tp,sex,covariates.index[covariates[coordinate] == 1][0]]+list(proportions[:,idx]))
  
  return pd.DataFrame(P_pretty,columns=['Slide','Coordinate','Genotype','Time point','Sex','Category']+list(P_df.index))

# "Expression on tissue sections" tab

# THESE ARE THE VARIABLES USER CAN CHANGE
#conditions = [['G93A','p100'],['G93A','p120'],['WT','p100'],['WT','p120']]
conditions = [['WT','p100']]
gene = 'Gfap'

# get the filenames
files = []
for condition in conditions:
  #print glob.glob('Estimates/CN_%s_%s_Y_norm.tsv'%(condition[0],condition[1])
  files.append('Estimates/CN_%s_%s_Y_norm.tsv'%(condition[0],condition[1]))

# read the metadata file
metadata = pd.read_table('Metadata/metadata_09252017.tsv',header=0,index_col=0)

# read result files
data_list = []
for filename in files:
  data_list.append(pd.read_table(filename,header=[0,1],index_col=0))

# get the minimum and maximum values
vmin = numpy.inf
vmax = -numpy.inf
for data in data_list:
  vmin = numpy.min([data[data.index == gene].as_matrix().min(),vmin])
  vmax = numpy.max([data[data.index == gene].as_matrix().max(),vmax])

# initialize the colormapper and ticker
log_scale = True
if log_scale:
  color_mapper = bokeh.models.mappers.LogColorMapper('Viridis256',low=vmin,high=vmax)
  ticker = bokeh.models.LogTicker(base=2)
else:
  color_mapper = bokeh.models.mappers.LinearColorMapper('Viridis256',low=vmin,high=vmax)
  ticker = bokeh.models.BasicTicker(base=2,mantissas=[1,5])

plots = []

# input thingies for user
textinput_gene = bokeh.models.widgets.TextInput(value='Gfap',title='Gene:')
checkbox_genotype = bokeh.models.widgets.CheckboxButtonGroup(labels=['WT','G93a'],name='Genotype:',active=[0,1])
checkbox_timepoint = bokeh.models.widgets.CheckboxButtonGroup(labels=['p30','p70','p100','p120'],name='Timepoint:',active=[2,3])
plots.append([textinput_gene])
plots.append([checkbox_genotype])
plots.append([checkbox_timepoint])

# loop over result files
for n in range(0,len(data_list)):
  subplots = []

  data = data_list[n]

  count_files = list(data.columns.levels[0])
  
  # loop over count files (arrays)
  for count_file_idx in range(0,len(count_files)):
    # extract coordinates
    coordinates = [map(float,foo.split('_')) for foo in data[count_files[count_file_idx]].columns]
    coordinates = numpy.array(coordinates)

    # prepare the array title
    tmp = metadata[metadata['Count file'] == count_files[count_file_idx].replace('Count_Tables/','')]
    gt,tp,sex = tmp['Genotype'].values[0],tmp['Timepoint'].values[0],tmp['Sex'].values[0]
    gt = gt[0:5]+(gt[5:] and '..')
    title = '%s %s %s (%s)'%(gt,tp,sex,count_files[count_file_idx].replace('_stdata_aligned_counts_IDs.txt','').replace('Count_Tables/',''))

    # read array annotations
    if re.match(r".*_[0-9]+_stdata_aligned_counts_IDs.txt",count_files[count_file_idx]):
      covariates = pd.read_table('Covariates/%s.tsv'%(re.sub(r"_[0-9]+_stdata_aligned_counts_IDs.txt",r"",count_files[count_file_idx].replace('Count_Tables/',''))),header=0,index_col=0)
    else:
      covariates = pd.read_table('Covariates/%s.tsv'%(count_files[count_file_idx].replace('_stdata_aligned_counts_IDs.txt','').replace('Count_Tables/','')),header=0,index_col=0)

    # prepare spot annotations for the hover thingie
    # this probably could be done without the for loop
    annotations = []
    for coordinate in data[count_files[count_file_idx]].columns:
      if coordinate in covariates.columns:
        annotations.append(covariates.index[covariates[coordinate] == 1][0])
      else:
        annotations.append('Undefined')

    # read the array he image
    image_filename = 'Images/'+count_files[count_file_idx].replace('_stdata_aligned_counts_IDs.txt','_small.jpg').replace('Count_Tables/','')
    tissue_image = Image.open('Images/'+count_files[count_file_idx].replace('_stdata_aligned_counts_IDs.txt','_small.jpg').replace('Count_Tables/','')).convert('RGBA')

    # mapping between pixels and spot coordinates
    xdim,ydim = tissue_image.size
    pixel_dim = 194.0/((6200.0)/(xdim/0.05))
    pixel_dim = pixel_dim * 0.05

    # tissue image will be reprented in the rgba space
    img = numpy.empty((ydim,xdim),dtype=numpy.uint32)
    view = img.view(dtype=numpy.uint8).reshape((ydim,xdim,4))
    view[:,:,:] = numpy.flipud(numpy.asarray(tissue_image))

    # spot expressions
    expression_values = data[data.index == gene][count_files[count_file_idx]].as_matrix()[0]

    # map the spot coordinates to pixels
    x_coordinates = pixel_dim*(coordinates[:,0]-1)
    y_coordinates = ydim-pixel_dim*(coordinates[:,1]-1)

    # relevant data
    data_spots = bokeh.models.ColumnDataSource({'x': x_coordinates, 'y': y_coordinates, 'z': expression_values, 'annotation': annotations})
    data_img = bokeh.models.ColumnDataSource({'image': [img]})

    # initialize hover thingie with expressions and annotations
    hover = bokeh.models.HoverTool(tooltips=[('Expression','@z'),('Annotation','@annotation')])

    # initialize a bokeh figure
    s = bokeh.plotting.figure(width=235,plot_height=250,title=title,x_range=(0,xdim),y_range=(0,ydim),match_aspect=True,aspect_scale=1,tools=[bokeh.models.tools.PanTool(),bokeh.models.tools.WheelZoomTool(),bokeh.models.tools.ResetTool(),hover])
    s.toolbar.logo = None
    s.toolbar_location = None
    s.axis.visible = False

    # plot the array he image
    #s.image_rgba(image='image',x=0,y=0,dw=xdim,dh=ydim,source=data_img)
    s.image_url(url=[image_filename],x=0,y=0,anchor='bottom_left',w=xdim,h=ydim)

    # plot the spots
    s.scatter(x='x',y='y',radius=5,fill_color={'field': 'z','transform': color_mapper},fill_alpha=0.8,line_color=None,source=data_spots)

    subplots.append(s)
  plots.append(subplots)

# initialize an empty bokeh figure for colorbar
s = bokeh.plotting.figure(width=250,plot_height=100,title=None,x_axis_location=None,y_axis_location=None,tools='pan,wheel_zoom,reset',min_border=0,outline_line_color=None)
color_bar = bokeh.models.ColorBar(color_mapper=color_mapper,ticker=ticker,label_standoff=6,border_line_color=None,location=(0,0),major_tick_line_color='black',title='Expression',orientation='horizontal')
s.add_layout(color_bar,'above')
s.toolbar.logo = None
s.toolbar_location = None

plots.append([s])
    
p = bokeh.layouts.gridplot(plots,merge_tools=True,toolbar_location='left')

# "Coefficients" tab

# read coefficents

# THESE ARE THE VARIABLES USER CAN CHANGE
gene = 'Gfap'

# read the file containing densities and other stuff
with open('Coefficients/data.pickle','r') as f:
  covariates,timepoints,genotypes,x,data = pickle.load(f)

# get (two) different colors for different (WT and G93a) genotypes
palette = [colorcet.rainbow[i*80] for i in range(2)]

plots = []

# an input thingie for user
textinput_gene = bokeh.models.widgets.TextInput(value='Gfap',title='Gene:')
plots.append([textinput_gene])

# get the maximum value (for setting y-axis ranges)
max_value = -numpy.inf
for covariate in covariates:
  for genotype in genotypes:
    for timepoint in timepoints:
      max_value = max([max_value,data[gene][covariate][genotype][timepoint].max()])

# loop over tissue categories
for idx,covariate in enumerate(covariates,start=0):
  # four subplots per row
  if idx % 4 == 0:
    subplots = []

  # initialize data objects for two genotypes
  source = {}
  source[genotypes[0]] = bokeh.models.ColumnDataSource({'x': x})
  source[genotypes[1]] = bokeh.models.ColumnDataSource({'x': x})

  # initialize a bokeh figure
  s = bokeh.plotting.figure(y_range=timepoints,x_range=(-5,5),width=200,plot_height=250,title=COVARIATES[covariate],tools='pan,wheel_zoom,reset')

  # loop over genotypes
  for genotype,color in zip(genotypes,palette):
    # loop over time points
    for timepoint in timepoints:
      # update the data object
      y = joy(timepoint,data[gene][covariate][genotype][timepoint],scale=1/max_value)
      source[genotype].add(y,timepoint)

      # plot the distribution
      # legend is only drawn in the first subplot
      if covariate == covariates[0] and timepoint == timepoints[0]:
        s.patch('x',timepoint,line_color='black',fill_color=color,alpha=0.6,source=source[genotype],line_width=0.5,legend=genotype)
      else:
        s.patch('x',timepoint,line_color='black',fill_color=color,alpha=0.6,source=source[genotype],line_width=0.5)

      # try to make the legend look nice
      s.legend.location = 'bottom_left'
      s.legend.label_text_font_size = '7pt'
      s.legend.background_fill_alpha = 0
      s.legend.border_line_alpha = 0
      s.legend.glyph_height = 10
      s.legend.glyph_width = 10
      s.legend.margin = 1
      s.legend.spacing = 2
      s.legend.padding = 2
      s.legend.orientation = 'horizontal'

  # some clipping will happen without this (p120)
  s.y_range.range_padding = 0.2

  subplots.append(s)

  # is the subplot row full or are running out of categories? 
  if (idx+1) % 4 == 0 or (idx+1) == len(covariates):
    plots.append(subplots)

p2 = bokeh.layouts.gridplot(plots,merge_tools=True,toolbar_location='left')

# "Cell types on tissue sections" tab

# THESE ARE THE VARIABLES USER CAN CHANGE
cell_type = 'Cell type #10' 
#conditions = [['G93A','p100'],['G93A','p120'],['WT','p100'],['WT','p120']]
conditions = [['WT','p100']]

plots = []

# input thingies for user
select_celltype = bokeh.models.widgets.Select(title='Cell type:', value='Cell type #1', options=['Cell type #%d'%(idx) for idx in range(1,11)])
checkbox_genotype = bokeh.models.widgets.CheckboxButtonGroup(labels=['WT','G93a'],name='Genotype:',active=[0,1])
checkbox_timepoint = bokeh.models.widgets.CheckboxButtonGroup(labels=['p30','p70','p100','p120'],name='Timepoint:',active=[2,3])
plots.append([select_celltype])
plots.append([checkbox_genotype])
plots.append([checkbox_timepoint])

# read the data
P_df = pd.read_table('Deconvolution/P_optimize_10.tsv',header=[0,1],index_col=0)
count_files = numpy.array(list(P_df.columns.levels[0]))

# get the count files that match our conditions
# loop over count files (arrays)
genotypes = {}
for count_file_idx in range(0,len(count_files)):
  tmp = metadata[metadata['Count file'] == count_files[count_file_idx].replace('Count_Tables/','')]
  gt,tp = tmp['Genotype'].values[0],tmp['Timepoint'].values[0]
  #gt = gt[0:5]+(gt[5:] and '..')
  if not genotypes.has_key(gt):
      genotypes[gt] = {}
  if not genotypes[gt].has_key(tp):
      genotypes[gt][tp] = []
  genotypes[gt][tp].append(count_file_idx)

# get the maximum value for colorbar
vmax = -numpy.inf
for condition in conditions:
  for count_file_idx in genotypes[condition[0]][condition[1]]:
    vmax = max([P_df[P_df.index == cell_type][count_files[count_file_idx]].as_matrix().max(),vmax])

color_mapper = bokeh.models.mappers.LinearColorMapper('Plasma256',low=0,high=vmax)
ticker = bokeh.models.BasicTicker(base=2,mantissas=[1,5])

for condition in conditions:
  subplots = []
  for count_file_idx in genotypes[condition[0]][condition[1]]:
    tmp = metadata[metadata['Count file'] == count_files[count_file_idx].replace('Count_Tables/','')]
    gt,tp,sex = tmp['Genotype'].values[0],tmp['Timepoint'].values[0],tmp['Sex'].values[0]
    gt = gt[0:5]+(gt[5:] and '..')
    title = '%s %s %s (%s)'%(gt,tp,sex,count_files[count_file_idx].replace('_stdata_aligned_counts_IDs.txt','').replace('Count_Tables/',''))
  
    # extract coordinates
    coordinates = [map(float,foo.split('_')) for foo in P_df[count_files[count_file_idx]].columns]
    coordinates = numpy.array(coordinates)
  
    # read the array he image
    image_filename = 'Images/'+count_files[count_file_idx].replace('_stdata_aligned_counts_IDs.txt','_small.jpg').replace('Count_Tables/','')
    tissue_image = Image.open('Images/'+count_files[count_file_idx].replace('_stdata_aligned_counts_IDs.txt','_small.jpg').replace('Count_Tables/','')).convert('RGBA')
  
    # mapping between pixels and spot coordinates
    xdim,ydim = tissue_image.size
    pixel_dim = 194.0/((6200.0)/(xdim/0.05))
    pixel_dim = pixel_dim * 0.05
  
    # tissue image will be reprented in the rgba space
    img = numpy.empty((ydim,xdim),dtype=numpy.uint32)
    view = img.view(dtype=numpy.uint8).reshape((ydim,xdim,4))
    view[:,:,:] = numpy.flipud(numpy.asarray(tissue_image))
  
    # cell type proportions
    proportion_values = P_df[P_df.index == cell_type][count_files[count_file_idx]].as_matrix()[0]
  
    # map the spot coordinates to pixels
    x_coordinates = pixel_dim*(coordinates[:,0]-1)
    y_coordinates = ydim-pixel_dim*(coordinates[:,1]-1)
  
    # relevant data
    data_spots = bokeh.models.ColumnDataSource({'x': x_coordinates, 'y': y_coordinates, 'z': proportion_values})
    data_img = bokeh.models.ColumnDataSource({'image': [img]})
  
    # initialize a bokeh figure
    s = bokeh.plotting.figure(width=235,plot_height=250,title=title,x_range=(0,xdim),y_range=(0,ydim),match_aspect=True,aspect_scale=1,tools=[bokeh.models.tools.PanTool(),bokeh.models.tools.WheelZoomTool(),bokeh.models.tools.ResetTool()])
    s.toolbar.logo = None
    s.toolbar_location = None
    s.axis.visible = False
  
    # plot the array he image
    #s.image_rgba(image='image',x=0,y=0,dw=xdim,dh=ydim,source=data_img)
    s.image_url(url=[image_filename],x=0,y=0,anchor='bottom_left',w=xdim,h=ydim)
  
    # plot the spots
    s.scatter(x='x',y='y',radius=5,fill_color={'field': 'z','transform': color_mapper},fill_alpha=0.8,line_color=None,source=data_spots)
  
    subplots.append(s)

  plots.append(subplots)

# initialize an empty bokeh figure for colorbar
s = bokeh.plotting.figure(width=250,plot_height=100,title=None,x_axis_location=None,y_axis_location=None,tools='pan,wheel_zoom,reset',min_border=0,outline_line_color=None)
color_bar = bokeh.models.ColorBar(color_mapper=color_mapper,ticker=ticker,label_standoff=6,border_line_color=None,location=(0,0),major_tick_line_color='black',title='Proportion',orientation='horizontal')
s.add_layout(color_bar,'above')
s.toolbar.logo = None
s.toolbar_location = None

plots.append([s])

p3 = bokeh.layouts.gridplot(plots,merge_tools=True,toolbar_location='left')

# "Cell type signatures" tab
# data for the table
M_df = pd.read_table('Deconvolution/M_optimize_10.tsv',header=0,index_col=0)
source_dict = {'Gene': list(M_df.index)}
for col in M_df.columns:
    source_dict[col] = M_df[col].values
    source_dict[col] = M_df[col].values

source = bokeh.models.ColumnDataSource(source_dict)

# columns of the table
formatter = formatter=bokeh.models.widgets.NumberFormatter(format='0.0000')
columns = [bokeh.models.widgets.TableColumn(field='Gene',title='Gene')]
for col in M_df.columns:
  columns.append(bokeh.models.widgets.TableColumn(field=col,title=col,formatter=formatter))

data_table = bokeh.models.widgets.DataTable(source=source,columns=columns,width=800,height=600,row_headers=False)

p4 = bokeh.layouts.widgetbox(data_table)

# "Cell types in anatomical sections" tab

# THESE ARE THE VARIABLES USER CAN CHANGE
anatomical_section = 'Vent_Horn'
cell_types = ['Cell type #5', 'Cell type #9']
genotypes = ['WT','G93A']
timepoints = ['p30','p70','p100','p120']

# read the data
P_df = pd.read_table('Deconvolution/P_optimize_10.tsv',header=[0,1],index_col=0)
count_files = numpy.array(list(P_df.columns.levels[0]))
covariates = list(pd.read_table('Covariates/CN68_E1.tsv',header=0,index_col=0).index)

# read the metadata file
metadata = pd.read_table('Metadata/metadata_09252017.tsv',header=0,index_col=0)

plots = []
select_as = bokeh.models.widgets.Select(title='Anatomical section:',value='Vent_Horn', options=covariates)
select_ct_1 = bokeh.models.widgets.Select(title='First cell type:',value='Cell type #1', options=list(P_df.index))
select_ct_2 = bokeh.models.widgets.Select(title='Second cell type:',value='Cell type #2', options=list(P_df.index))
plots.append([select_as,select_ct_1,select_ct_2])

# let parse our dataframe
data = parse_P(P_df,count_files,metadata)

subplots = []
for timepoint in timepoints:
  data_x_gt1 = data[(data['Genotype'] == genotypes[0]) & (data['Time point'] == timepoint) & (data['Category'] == anatomical_section)][cell_types[0]]
  data_y_gt1 = data[(data['Genotype'] == genotypes[0]) & (data['Time point'] == timepoint) & (data['Category'] == anatomical_section)][cell_types[1]]

  data_x_gt2 = data[(data['Genotype'] == genotypes[1]) & (data['Time point'] == timepoint) & (data['Category'] == anatomical_section)][cell_types[0]]
  data_y_gt2 = data[(data['Genotype'] == genotypes[1]) & (data['Time point'] == timepoint) & (data['Category'] == anatomical_section)][cell_types[1]]

  data_spots_gt1 = bokeh.models.ColumnDataSource({'data_x': data_x_gt1, 'data_y': data_y_gt1})
  data_spots_gt2 = bokeh.models.ColumnDataSource({'data_x': data_x_gt2, 'data_y': data_y_gt2})

  # initialize a bokeh figure
  s = bokeh.plotting.figure(width=235,plot_height=250,title=timepoint,x_range=(0,1),y_range=(0,1),x_axis_label=cell_types[0],y_axis_label=cell_types[1],match_aspect=True,aspect_scale=1,tools=[bokeh.models.tools.PanTool(),bokeh.models.tools.WheelZoomTool(),bokeh.models.tools.ResetTool()])
  s.toolbar.logo = None
  s.toolbar_location = None

  # plot the spots
  if timepoint == timepoints[0]:
    s.scatter(x='data_x',y='data_y',radius=0.005,fill_color='red',fill_alpha=0.5,line_color=None,source=data_spots_gt1,legend=genotypes[0])
    s.scatter(x='data_x',y='data_y',radius=0.005,fill_color='blue',fill_alpha=0.5,line_color=None,source=data_spots_gt2,legend=genotypes[1])
  else:
    s.scatter(x='data_x',y='data_y',radius=0.005,fill_color='red',fill_alpha=0.5,line_color=None,source=data_spots_gt1)
    s.scatter(x='data_x',y='data_y',radius=0.005,fill_color='blue',fill_alpha=0.5,line_color=None,source=data_spots_gt2)

  s.legend.location = 'top_left'
  s.legend.label_text_font_size = '7pt'
  s.legend.background_fill_alpha = 0
  s.legend.border_line_alpha = 0
  s.legend.glyph_height = 10
  s.legend.glyph_width = 10
  s.legend.margin = 1
  s.legend.spacing = 2
  s.legend.padding = 2
  s.legend.orientation = 'horizontal'

  subplots.append(s)

plots.append(subplots)

p5 = bokeh.layouts.gridplot(plots,merge_tools=True,toolbar_location='left')

# we have selectable tabs
tab1 = bokeh.models.widgets.Panel(child=p, title="Expression on tissue sections")
tab2 = bokeh.models.widgets.Panel(child=p2, title="Expression coefficients")
tab3 = bokeh.models.widgets.Panel(child=p3, title="Cell types on tissue sections")
tab4 = bokeh.models.widgets.Panel(child=p4, title="Cell type signatures")
tab5 = bokeh.models.widgets.Panel(child=p5, title="Cell types in anatomical sections")
tabs = bokeh.models.widgets.Tabs(tabs=[tab1,tab2,tab3,tab4,tab5])

# finally, let us display our tab pane object
bokeh.plotting.show(tabs)
