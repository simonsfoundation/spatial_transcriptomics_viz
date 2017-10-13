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

import argparse

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
      covariates = pd.read_table('../data/Covariates/%s.tsv'%(re.sub(r"_[0-9]+_stdata_aligned_counts_IDs.txt",r"",filename.replace('Count_Tables/',''))),header=0,index_col=0)
    else:
      covariates = pd.read_table('../data/Covariates/%s.tsv'%(filename.replace('_stdata_aligned_counts_IDs.txt','').replace('Count_Tables/','')),header=0,index_col=0)
 
    for idx,coordinate in enumerate(coordinates):
      P_pretty.append([name,coordinate,gt,tp,sex,covariates.index[covariates[coordinate] == 1][0]]+list(proportions[:,idx]))
  
  return pd.DataFrame(P_pretty,columns=['Slide','Coordinate','Genotype','Time point','Sex','Category']+list(P_df.index))

def get_expression_on_tissue_sections_data(gene,conditions,metadata_filename='../data/Metadata/metadata_09252017.tsv'):
  # get the filenames
  files = []
  for condition in conditions:
    #print glob.glob('Estimates/CN_%s_%s_Y_norm.tsv'%(condition[0],condition[1])
    files.append('../data/Estimates/CN_%s_%s_Y_norm.tsv'%(condition[0],condition[1]))
  
  # read the metadata file
  metadata = pd.read_table(metadata_filename,header=0,index_col=0)
  
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

  viz_data = []

  # loop over result files
  for n in range(0,len(data_list)):
    viz_data.append([])
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
        covariates = pd.read_table('../data/Covariates/%s.tsv'%(re.sub(r"_[0-9]+_stdata_aligned_counts_IDs.txt",r"",count_files[count_file_idx].replace('Count_Tables/',''))),header=0,index_col=0)
      else:
        covariates = pd.read_table('../data/Covariates/%s.tsv'%(count_files[count_file_idx].replace('_stdata_aligned_counts_IDs.txt','').replace('Count_Tables/','')),header=0,index_col=0)
  
      # prepare spot annotations for the hover thingie
      # this probably could be done without the for loop
      annotations = []
      for coordinate in data[count_files[count_file_idx]].columns:
        if coordinate in covariates.columns:
          annotations.append(covariates.index[covariates[coordinate] == 1][0])
        else:
          annotations.append('Undefined')
  
      # read the array he image
      image_filename = '../data/Images/'+count_files[count_file_idx].replace('_stdata_aligned_counts_IDs.txt','_small.jpg').replace('Count_Tables/','')
      tissue_image = Image.open('../data/Images/'+count_files[count_file_idx].replace('_stdata_aligned_counts_IDs.txt','_small.jpg').replace('Count_Tables/','')).convert('RGBA')
  
      # mapping between pixels and spot coordinates
      xdim,ydim = tissue_image.size
      pixel_dim = 194.0/((6200.0)/(xdim/0.05))
      pixel_dim = pixel_dim * 0.05
  
      # spot expressions
      expression_values = data[data.index == gene][count_files[count_file_idx]].as_matrix()[0]
  
      # map the spot coordinates to pixels
      x_coordinates = pixel_dim*(coordinates[:,0]-1)
      y_coordinates = ydim-pixel_dim*(coordinates[:,1]-1)
  
      # relevant data
      data_spots = bokeh.models.ColumnDataSource({'x': x_coordinates, 'y': y_coordinates, 'z': expression_values, 'annotation': annotations})
      data_img = bokeh.models.ColumnDataSource({'image': [image_filename]})
  
      viz_data[n].append({'title': title,'data_spots': data_spots,'data_img': data_img,'dims': [xdim,ydim]})

  return viz_data

def viz_expression_on_tissue_sections_data(viz_data,log_scale=True):
  # get the minimum and maximum values
  vmin = numpy.inf
  vmax = -numpy.inf
  # loop over conditions
  for condition_data in viz_data:
    # loop over count files (arrays)
    for array_data in condition_data:
      vmin = numpy.min([array_data['data_spots'].data['z'].min(),vmin])
      vmax = numpy.max([array_data['data_spots'].data['z'].max(),vmax])

  # initialize the colormapper and ticker
  if log_scale:
    color_mapper = bokeh.models.mappers.LogColorMapper('Viridis256',low=vmin,high=vmax)
    ticker = bokeh.models.LogTicker(base=2)
  else:
    color_mapper = bokeh.models.mappers.LinearColorMapper('Viridis256',low=vmin,high=vmax)
    ticker = bokeh.models.BasicTicker(base=2,mantissas=[1,5])
 
  plots = []

  # input thingies for user
  textinput_gene = bokeh.models.widgets.TextInput(value='Gfap',title='Gene:')
  #checkbox_genotype = bokeh.models.widgets.CheckboxButtonGroup(labels=['WT','G93a'],name='Genotype:',active=[0,1])
  #checkbox_timepoint = bokeh.models.widgets.CheckboxButtonGroup(labels=['p30','p70','p100','p120'],name='Timepoint:',active=[2,3])
  plots.append([textinput_gene])
  #plots.append([checkbox_genotype])
  #plots.append([checkbox_timepoint])

  # loop over result files
  for condition_data in viz_data:
    subplots = []
    
    # loop over count files (arrays)
    for array_data in condition_data:
  
      # initialize hover thingie with expressions and annotations
      hover = bokeh.models.HoverTool(tooltips=[('Expression','@z'),('Annotation','@annotation')])
  
      # initialize a bokeh figure
      s = bokeh.plotting.figure(width=235,plot_height=250,title=array_data['title'],x_range=(0,array_data['dims'][0]),y_range=(0,array_data['dims'][1]),match_aspect=True,aspect_scale=1,tools=[bokeh.models.tools.PanTool(),bokeh.models.tools.WheelZoomTool(),bokeh.models.tools.ResetTool(),hover])
      s.toolbar.logo = None
      s.toolbar_location = None
      s.axis.visible = False
  
      # plot the array he image
      #s.image_rgba(image='image',x=0,y=0,dw=xdim,dh=ydim,source=data_img)
      s.image_url(url='image',x=0,y=0,anchor='bottom_left',w=array_data['dims'][0],h=array_data['dims'][1],source=array_data['data_img'])
  
      # plot the spots
      s.scatter(x='x',y='y',radius=5,fill_color={'field': 'z','transform': color_mapper},fill_alpha=0.8,line_color=None,source=array_data['data_spots'])

      subplots.append(s)
    plots.append(subplots)
  
  # initialize an empty bokeh figure for colorbar
  s = bokeh.plotting.figure(width=250,plot_height=100,title=None,x_axis_location=None,y_axis_location=None,tools='pan,wheel_zoom,reset',min_border=0,outline_line_color=None)
  color_bar = bokeh.models.ColorBar(color_mapper=color_mapper,ticker=ticker,label_standoff=6,border_line_color=None,location=(0,0),major_tick_line_color='black',title='Expression',orientation='horizontal')
  s.add_layout(color_bar,'above')
  s.toolbar.logo = None
  s.toolbar_location = None
  
  plots.append([s])
  
  return plots

def get_expression_coefficients_data(filename='../data/Coefficients/data.pickle'):
  # read the file containing densities and other stuff
  with open(filename,'r') as f:
    covariates,timepoints,genotypes,x,data = pickle.load(f)
  return {'covariates': covariates,'timepoints': timepoints,'genotypes': genotypes,'x': x,'data': data}

  # MAYBE WE SHOULD RETURN COLUMNDATASOURCES

def viz_expression_coefficients_data(viz_data,gene,n_cols=4):
  # only works with two genotypes currently
  # get (two) different colors for different (WT and G93a) genotypes
  palette = [colorcet.rainbow[i*80] for i in range(2)]
  
  plots = []
  
  # an input thingie for user
  textinput_gene = bokeh.models.widgets.TextInput(value=gene,title='Gene:')
  plots.append([textinput_gene])
  
  # get the maximum value (for setting y-axis ranges)
  max_value = -numpy.inf
  for covariate in viz_data['covariates']:
    for genotype in viz_data['genotypes']:
      for timepoint in viz_data['timepoints']:
        max_value = max([max_value,viz_data['data'][gene][covariate][genotype][timepoint].max()])
  
  # loop over tissue categories
  for idx,covariate in enumerate(viz_data['covariates'],start=0):
    # four subplots per row
    if idx % n_cols == 0:
      subplots = []
  
    # initialize data objects for two genotypes
    source = {}
    source[viz_data['genotypes'][0]] = bokeh.models.ColumnDataSource({'x': viz_data['x']})
    source[viz_data['genotypes'][1]] = bokeh.models.ColumnDataSource({'x': viz_data['x']})
  
    # initialize a bokeh figure
    s = bokeh.plotting.figure(y_range=viz_data['timepoints'],x_range=(-5,5),width=200,plot_height=250,title=COVARIATES[covariate],tools='pan,wheel_zoom,reset')
  
    # loop over genotypes
    for genotype,color in zip(viz_data['genotypes'],palette):
      # loop over time points
      for timepoint in viz_data['timepoints']:
        # update the data object
        y = joy(timepoint,viz_data['data'][gene][covariate][genotype][timepoint],scale=1/max_value)
        source[genotype].add(y,timepoint)
  
        # plot the distribution
        # legend is only drawn in the first subplot
        if covariate == viz_data['covariates'][0] and timepoint == viz_data['timepoints'][0]:
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
    if (idx+1) % n_cols == 0 or (idx+1) == len(viz_data['covariates']):
      plots.append(subplots)

  return plots

def get_cell_types_on_tissue_sections_data(cell_type,conditions,proportions_filename='../data/Deconvolution/P_optimize_10.tsv',metadata_filename='../data/Metadata/metadata_09252017.tsv'):
  # read the data
  P_df = pd.read_table(proportions_filename,header=[0,1],index_col=0)
  count_files = numpy.array(list(P_df.columns.levels[0]))
  
  # read the metadata file
  metadata = pd.read_table(metadata_filename,header=0,index_col=0)
  
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
  
  viz_data = []

  for n,condition in enumerate(conditions,start=0):

    viz_data.append([])

    for count_file_idx in genotypes[condition[0]][condition[1]]:
      tmp = metadata[metadata['Count file'] == count_files[count_file_idx].replace('Count_Tables/','')]
      gt,tp,sex = tmp['Genotype'].values[0],tmp['Timepoint'].values[0],tmp['Sex'].values[0]
      gt = gt[0:5]+(gt[5:] and '..')
      title = '%s %s %s (%s)'%(gt,tp,sex,count_files[count_file_idx].replace('_stdata_aligned_counts_IDs.txt','').replace('Count_Tables/',''))
    
      # extract coordinates
      coordinates = [map(float,foo.split('_')) for foo in P_df[count_files[count_file_idx]].columns]
      coordinates = numpy.array(coordinates)
    
      # read the array he image
      image_filename = '../data/Images/'+count_files[count_file_idx].replace('_stdata_aligned_counts_IDs.txt','_small.jpg').replace('Count_Tables/','')
      tissue_image = Image.open('../data/Images/'+count_files[count_file_idx].replace('_stdata_aligned_counts_IDs.txt','_small.jpg').replace('Count_Tables/','')).convert('RGBA')
    
      # mapping between pixels and spot coordinates
      xdim,ydim = tissue_image.size
      pixel_dim = 194.0/((6200.0)/(xdim/0.05))
      pixel_dim = pixel_dim * 0.05
    
      # cell type proportions
      proportion_values = P_df[P_df.index == cell_type][count_files[count_file_idx]].as_matrix()[0]
    
      # map the spot coordinates to pixels
      x_coordinates = pixel_dim*(coordinates[:,0]-1)
      y_coordinates = ydim-pixel_dim*(coordinates[:,1]-1)
    
      # relevant data
      data_spots = bokeh.models.ColumnDataSource({'x': x_coordinates, 'y': y_coordinates, 'z': proportion_values})
      data_img = bokeh.models.ColumnDataSource({'image': [image_filename]})

      viz_data[n].append({'title': title,'data_spots': data_spots,'data_img': data_img,'dims': [xdim,ydim]})
    
  return viz_data 

def viz_cell_types_on_tissue_sections_data(viz_data):
  # get the maximum value for colorbar
  vmax = -numpy.inf
  # loop over conditions
  for condition_data in viz_data:
    # loop over count files (arrays)
    for array_data in condition_data:
      vmax = numpy.max([array_data['data_spots'].data['z'].max(),vmax])
  
  color_mapper = bokeh.models.mappers.LinearColorMapper('Plasma256',low=0,high=vmax)
  ticker = bokeh.models.BasicTicker(base=2,mantissas=[1,5])

  plots = []

  # input thingies for user
  select_celltype = bokeh.models.widgets.Select(title='Cell type:', value='Cell type #1', options=['Cell type #%d'%(idx) for idx in range(1,11)])
  #checkbox_genotype = bokeh.models.widgets.CheckboxButtonGroup(labels=['WT','G93a'],name='Genotype:',active=[0,1])
  #checkbox_timepoint = bokeh.models.widgets.CheckboxButtonGroup(labels=['p30','p70','p100','p120'],name='Timepoint:',active=[2,3])
  plots.append([select_celltype])
  #plots.append([checkbox_genotype])
  #plots.append([checkbox_timepoint])
  
  for condition_data in viz_data:
    subplots = []
    for array_data in condition_data:
      # initialize a bokeh figure
      s = bokeh.plotting.figure(width=235,plot_height=250,title=array_data['title'],x_range=(0,array_data['dims'][0]),y_range=(0,array_data['dims'][1]),match_aspect=True,aspect_scale=1,tools=[bokeh.models.tools.PanTool(),bokeh.models.tools.WheelZoomTool(),bokeh.models.tools.ResetTool()])
      s.toolbar.logo = None
      s.toolbar_location = None
      s.axis.visible = False
    
      # plot the array he image
      s.image_url(url='image',x=0,y=0,anchor='bottom_left',w=array_data['dims'][0],h=array_data['dims'][1],source=array_data['data_img'])
    
      # plot the spots
      s.scatter(x='x',y='y',radius=5,fill_color={'field': 'z','transform': color_mapper},fill_alpha=0.8,line_color=None,source=array_data['data_spots'])
    
      subplots.append(s)
  
    plots.append(subplots)

  # initialize an empty bokeh figure for colorbar
  s = bokeh.plotting.figure(width=250,plot_height=100,title=None,x_axis_location=None,y_axis_location=None,tools='pan,wheel_zoom,reset',min_border=0,outline_line_color=None)
  color_bar = bokeh.models.ColorBar(color_mapper=color_mapper,ticker=ticker,label_standoff=6,border_line_color=None,location=(0,0),major_tick_line_color='black',title='Proportion',orientation='horizontal')
  s.add_layout(color_bar,'above')
  s.toolbar.logo = None
  s.toolbar_location = None
  
  plots.append([s])

  return plots

def get_cell_type_signatures_data(filename='../data/Deconvolution/M_optimize_10.tsv'):
  # data for the table
  M_df = pd.read_table(filename,header=0,index_col=0)
  source_dict = {'Gene': list(M_df.index)}
  for col in M_df.columns:
      source_dict[col] = M_df[col].values
      source_dict[col] = M_df[col].values
  
  return {'columns': list(M_df.columns),'data': bokeh.models.ColumnDataSource(source_dict)}

def viz_cell_type_signatures_data(viz_data):
  # columns of the table
  formatter = formatter=bokeh.models.widgets.NumberFormatter(format='0.0000')
  columns = [bokeh.models.widgets.TableColumn(field='Gene',title='Gene')]
  for col in viz_data['columns']:
    columns.append(bokeh.models.widgets.TableColumn(field=col,title=col,formatter=formatter))
  
  return bokeh.models.widgets.DataTable(source=viz_data['data'],columns=columns,width=800,height=600,row_headers=False)

def get_cell_types_in_anatomical_sections_data(anatomical_section,cell_types,timepoints,genotypes,proportions_filename='../data/Deconvolution/P_optimize_10.tsv',metadata_filename='../data/Metadata/metadata_09252017.tsv'):
  # read the data
  P_df = pd.read_table(proportions_filename,header=[0,1],index_col=0)
  count_files = numpy.array(list(P_df.columns.levels[0]))
  covariates = list(pd.read_table('../data/Covariates/CN68_E1.tsv',header=0,index_col=0).index)
  
  # read the metadata file
  metadata = pd.read_table(metadata_filename,header=0,index_col=0)

  # let parse our dataframe
  data = parse_P(P_df,count_files,metadata)

  viz_data = {'anatomical_section': anatomical_section,'cell_types': list(P_df.index),'covariates': covariates, 'genotypes': genotypes, 'timepoints': timepoints, 'data': []}

  for timepoint in timepoints:
    data_x_gt1 = data[(data['Genotype'] == genotypes[0]) & (data['Time point'] == timepoint) & (data['Category'] == anatomical_section)][cell_types[0]]
    data_y_gt1 = data[(data['Genotype'] == genotypes[0]) & (data['Time point'] == timepoint) & (data['Category'] == anatomical_section)][cell_types[1]]
  
    data_x_gt2 = data[(data['Genotype'] == genotypes[1]) & (data['Time point'] == timepoint) & (data['Category'] == anatomical_section)][cell_types[0]]
    data_y_gt2 = data[(data['Genotype'] == genotypes[1]) & (data['Time point'] == timepoint) & (data['Category'] == anatomical_section)][cell_types[1]]
  
    data_spots_gt1 = bokeh.models.ColumnDataSource({'data_x': data_x_gt1, 'data_y': data_y_gt1})
    data_spots_gt2 = bokeh.models.ColumnDataSource({'data_x': data_x_gt2, 'data_y': data_y_gt2})

    viz_data['data'].append({'data_spots': [data_spots_gt1, data_spots_gt2],'title': timepoint})
  
  return viz_data

def viz_cell_types_in_anatomical_sections_data(viz_data):
  plots = []

  select_as = bokeh.models.widgets.Select(title='Anatomical section:',value=viz_data['anatomical_section'], options=viz_data['covariates'])
  select_ct_1 = bokeh.models.widgets.Select(title='First cell type:',value=viz_data['cell_types'][0], options=viz_data['cell_types'])
  select_ct_2 = bokeh.models.widgets.Select(title='Second cell type:',value=viz_data['cell_types'][1], options=viz_data['cell_types'])

  plots.append([select_as,select_ct_1,select_ct_2])
  
  subplots = []
  for timepoint,col_data in zip(viz_data['timepoints'],viz_data['data']):
  
    # initialize a bokeh figure
    s = bokeh.plotting.figure(width=235,plot_height=250,title=col_data['title'],x_range=(0,1),y_range=(0,1),x_axis_label=viz_data['cell_types'][0],y_axis_label=viz_data['cell_types'][1],match_aspect=True,aspect_scale=1,tools=[bokeh.models.tools.PanTool(),bokeh.models.tools.WheelZoomTool(),bokeh.models.tools.ResetTool()])
    s.toolbar.logo = None
    s.toolbar_location = None
  
    # plot the spots
    if timepoint == timepoints[0]:
      s.scatter(x='data_x',y='data_y',radius=0.005,fill_color='red',fill_alpha=0.5,line_color=None,source=col_data['data_spots'][0],legend=viz_data['genotypes'][0])
      s.scatter(x='data_x',y='data_y',radius=0.005,fill_color='blue',fill_alpha=0.5,line_color=None,source=col_data['data_spots'][1],legend=viz_data['genotypes'][1])
    else:
      s.scatter(x='data_x',y='data_y',radius=0.005,fill_color='red',fill_alpha=0.5,line_color=None,source=col_data['data_spots'][0])
      s.scatter(x='data_x',y='data_y',radius=0.005,fill_color='blue',fill_alpha=0.5,line_color=None,source=col_data['data_spots'][1])
  
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

  return plots

if __name__ == '__main__':
  parser = argparse.ArgumentParser(description='Slice report (by default all the tabs are drawn)')
  parser.add_argument('-a','--expression_on_tissues',action='store_true',dest='expression_on_tissues',help='Draw the expression on tissue sections tab')
  parser.add_argument('-b','--expression_coefficients',action='store_true',dest='expression_coefficients',help='Draw the expression coefficients tab')
  parser.add_argument('-c','--cell_types_on_tissues',action='store_true',dest='cell_types_on_tissues',help='Draw the cell types on tissue sections tab')
  parser.add_argument('-d','--cell_type_signatures',action='store_true',dest='cell_type_signatures',help='Draw the cell type signatures tab')
  parser.add_argument('-e','--cell_types_in_anatomical_sections',action='store_true',dest='cell_types_in_anatomical_sections',help='Draw the cell types in anatomical sections tab')
  parser.add_argument('-v','--version',action='version',version='%(prog)s 0.666')
  options = parser.parse_args()

  # by default all the tabs are drawn
  if not options.expression_on_tissues and not options.expression_coefficients and not options.cell_types_on_tissues and not options.cell_type_signatures and not options.cell_types_in_anatomical_sections:
    options.cell_types_on_tissues = True
    options.expression_coefficients = True
    options.cell_types_on_tissues = True
    options.cell_type_signatures = True
    options.cell_types_in_anatomical_sections = True

  tab_list = []

  if options.expression_on_tissues:
    conditions = []
    for genotype in ['WT','G93A']:
      for timepoint in ['p030','p070','p100','p120']:
        conditions.append([genotype,timepoint])
    #conditions = [['WT','p100']]
    #conditions = [['WT','p100']]
    gene = 'Gfap'
  
    viz_data = get_expression_on_tissue_sections_data(gene,conditions)
    plots = viz_expression_on_tissue_sections_data(viz_data)
    p = bokeh.layouts.gridplot(plots,merge_tools=True,toolbar_location='left')
    tab_list.append(bokeh.models.widgets.Panel(child=p,title='Expression on tissue sections'))
  
  if options.expression_coefficients:
    gene = 'Gfap'
  
    viz_data = get_expression_coefficients_data()
    plots = viz_expression_coefficients_data(viz_data,gene)
    p2 = bokeh.layouts.gridplot(plots,merge_tools=True,toolbar_location='left')
    tab_list.append(bokeh.models.widgets.Panel(child=p2,title='Expression coefficients'))
  
  if options.cell_types_on_tissues:
    cell_type = 'Cell type #10' 
    conditions = []
    for genotype in ['WT','G93A']:
      for timepoint in ['p30','p70','p100','p120']:
        conditions.append([genotype,timepoint])
    #conditions = [['G93A','p100']]
  
    viz_data = get_cell_types_on_tissue_sections_data(cell_type,conditions)
    plots = viz_cell_types_on_tissue_sections_data(viz_data)
    p3 = bokeh.layouts.gridplot(plots,merge_tools=True,toolbar_location='left')
    tab_list.append(bokeh.models.widgets.Panel(child=p3,title='Cell types on tissue sections'))
  
  if options.cell_type_signatures:
    viz_data = get_cell_type_signatures_data()
    data_table = viz_cell_type_signatures_data(viz_data)
    p4 = bokeh.layouts.widgetbox(data_table)
    tab_list.append(bokeh.models.widgets.Panel(child=p4,title='Cell type signatures'))
  
  if options.cell_types_in_anatomical_sections:
    anatomical_section = 'Vent_Horn'
    cell_types = ['Cell type #5', 'Cell type #9']
    genotypes = ['WT','G93A']
    timepoints = ['p30','p70','p100','p120']
  
    viz_data = get_cell_types_in_anatomical_sections_data(anatomical_section,cell_types,timepoints,genotypes)
    plots = viz_cell_types_in_anatomical_sections_data(viz_data)
    p5 = bokeh.layouts.gridplot(plots,merge_tools=True,toolbar_location='left')
    tab_list.append(bokeh.models.widgets.Panel(child=p5,title='Cell types in anatomical sections'))
  
  # create our tab pane object from our tabs
  tabs = bokeh.models.widgets.Tabs(tabs=tab_list)
  
  # finally, let us display our tab pane object
  bokeh.plotting.show(tabs)
