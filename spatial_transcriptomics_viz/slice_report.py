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

class ExpressionOnTissueSections:
  def __init__(self,gene,conditions,metadata_filename='../data/Metadata/metadata_09252017.tsv',log_scale=True):
    self.genes,self.data = self.__read_data(conditions,metadata_filename)
    self.source_spots,self.vmin,self.vmax = self.__create_source_spots(gene)
    self.source_image = self.__create_source_image()

    # initialize the colormapper and ticker
    if log_scale:
      self.color_mapper = bokeh.models.mappers.LogColorMapper('Viridis256',low=self.vmin,high=self.vmax)
      self.ticker = bokeh.models.LogTicker(base=2)
    else:
      self.color_mapper = bokeh.models.mappers.LinearColorMapper('Viridis256',low=self.vmin,high=self.vmax)
      self.ticker = bokeh.models.BasicTicker(base=2,mantissas=[1,5])

    self.error_pretext = bokeh.models.widgets.PreText(text='',width=100,height=20)

  def __read_data(self,conditions,metadata_filename):
    # names of the results files of the given conditions
    filenames = []
    for condition in conditions:
      filenames.append('../data/Estimates/CN_%s_%s_Y_norm.tsv'%(condition[0],condition[1]))

    # read the results files
    results = []
    for filename in filenames:
      results.append(pd.read_table(filename,header=[0,1],index_col=0))

    # read the metadata file
    metadata = pd.read_table(metadata_filename,header=0,index_col=0)

    # get the gene names
    genes = list(results[0].index)

    data = []

    # loop over the result files
    for n,result in enumerate(results,start=0):
      data.append([])

      # get the count filenames of the condition
      count_files = list(result.columns.levels[0])
      
      # loop over the count files (arrays) included in the result file
      for count_file_idx in range(0,len(count_files)):

        # extract the coordinates of the spots on the array
        coordinates = [map(float,foo.split('_')) for foo in result[count_files[count_file_idx]].columns]
        coordinates = numpy.array(coordinates)
    
        # prepare the array title
        tmp = metadata[metadata['Count file'] == count_files[count_file_idx].replace('Count_Tables/','')]
        gt,tp,sex = tmp['Genotype'].values[0],tmp['Timepoint'].values[0],tmp['Sex'].values[0]
        # truncate long genotype names
        gt = gt[0:5]+(gt[5:] and '..')
        title = '%s %s %s (%s)'%(gt,tp,sex,count_files[count_file_idx].replace('_stdata_aligned_counts_IDs.txt','').replace('Count_Tables/',''))
    
        # read the annotations of the spots on the array
        if re.match(r".*_[0-9]+_stdata_aligned_counts_IDs.txt",count_files[count_file_idx]):
          covariates = pd.read_table('../data/Covariates/%s.tsv'%(re.sub(r"_[0-9]+_stdata_aligned_counts_IDs.txt",r"",count_files[count_file_idx].replace('Count_Tables/',''))),header=0,index_col=0)
        else:
          covariates = pd.read_table('../data/Covariates/%s.tsv'%(count_files[count_file_idx].replace('_stdata_aligned_counts_IDs.txt','').replace('Count_Tables/','')),header=0,index_col=0)
    
        # prepare the spot annotations for the hover thingie
        annotations = []
        for coordinate in result[count_files[count_file_idx]].columns:
          if coordinate in covariates.columns:
            annotations.append(covariates.index[covariates[coordinate] == 1][0])
          # this should not happen
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
        #expression_values = result[result.index == gene][count_files[count_file_idx]].as_matrix()[0]
        expressions = result[count_files[count_file_idx]].as_matrix()
    
        # map the spot coordinates to pixels
        x_coordinates = pixel_dim*(coordinates[:,0]-1)
        y_coordinates = ydim-pixel_dim*(coordinates[:,1]-1)
    
        # relevant data
        data[n].append({'title': title,'image_filename': image_filename,'annotations': annotations,'expressions': expressions,'coordinates': [x_coordinates,y_coordinates],'dims': [xdim,ydim]})

    return genes,data

  def __create_source_spots(self,gene):
    # get the minimum and maximum values
    vmin = numpy.inf
    vmax = -numpy.inf

    # loop over conditions
    for condition_data in self.data:
      # loop over count files (arrays)
      for array_data in condition_data:
        vmin = numpy.min([array_data['expressions'][self.genes.index(gene),:].min(),vmin])
        vmax = numpy.max([array_data['expressions'][self.genes.index(gene),:].max(),vmax])
  
    source_spots = []
    for n,condition_data in enumerate(self.data,start=0):
      source_spots.append([])
      for array_data in condition_data:
        source_spots[n].append(bokeh.models.ColumnDataSource({'x': array_data['coordinates'][0], 'y': array_data['coordinates'][1], 'z': array_data['expressions'][self.genes.index(gene),:], 'annotation': array_data['annotations']}))
    return source_spots, vmin, vmax

  def __create_source_image(self):
    source_image = []
    for n,condition_data in enumerate(self.data,start=0):
      source_image.append([])
      for array_data in condition_data:
        source_image[n].append(bokeh.models.ColumnDataSource({'image': [array_data['image_filename']]}))
    return source_image

  def __update_plot(self,attr,old,new):
    if new not in self.genes:
      self.error_pretext.text = 'gene not found!'
      return
    self.error_pretext.text = ''

    source_spots,vmin,vmax = self.__create_source_spots(new)

    for n in range(0,len(source_spots)):
      for m in range(0,len(source_spots[n])):
        self.source_spots[n][m].data = source_spots[n][m].data

    self.color_mapper.low = vmin
    self.color_mapper.high = vmax

    self.vmax_value = vmax
    self.vmin_value = vmin
  
  def plot(self,gene='Gfap'):
    plots = []
  
    # input thingies for user
    textinput_gene = bokeh.models.widgets.TextInput(value='Gfap',title='Gene:')
    textinput_gene.on_change('value',self.__update_plot)
    plots.append([textinput_gene,self.error_pretext])
  
    # loop over result files
    for n,condition_data in enumerate(self.data,start=0):
      subplots = []
      
      # loop over count files (arrays)
      for m,array_data in enumerate(condition_data,start=0):
    
        # initialize hover thingie with expressions and annotations
        hover = bokeh.models.HoverTool(tooltips=[('Expression','@z'),('Annotation','@annotation')])
    
        # initialize a bokeh figure
        s = bokeh.plotting.figure(width=235,plot_height=250,title=array_data['title'],x_range=(0,array_data['dims'][0]),y_range=(0,array_data['dims'][1]),match_aspect=True,aspect_scale=1,tools=[bokeh.models.tools.PanTool(),bokeh.models.tools.WheelZoomTool(),bokeh.models.tools.ResetTool(),hover])
        s.toolbar.logo = None
        s.toolbar_location = None
        s.axis.visible = False
    
        # plot the h&e stained image of the array
        s.image_url(url='image',x=0,y=0,anchor='bottom_left',w=array_data['dims'][0],h=array_data['dims'][1],source=self.source_image[n][m])
    
        # plot the spots
        s.scatter(x='x',y='y',radius=5,fill_color={'field': 'z','transform': self.color_mapper},fill_alpha=0.8,line_color=None,source=self.source_spots[n][m])
  
        subplots.append(s)

      plots.append(subplots)
    
    # initialize an empty bokeh figure for colorbar
    s = bokeh.plotting.figure(width=250,plot_height=100,title=None,x_axis_location=None,y_axis_location=None,tools='pan,wheel_zoom,reset',min_border=0,outline_line_color=None)
    color_bar = bokeh.models.ColorBar(color_mapper=self.color_mapper,ticker=self.ticker,label_standoff=6,border_line_color=None,location=(0,0),major_tick_line_color='black',title='Expression',orientation='horizontal')
    s.add_layout(color_bar,'above')
    s.toolbar.logo = None
    s.toolbar_location = None
    
    plots.append([s])
    
    return plots

class ExpressionCoefficients:
  def __init__(self,gene='Gfap',filename='../data/Coefficients/data.pickle'):
    self.covariates,self.timepoints,self.genotypes,self.x,self.data = self.__read_data(filename)
    self.source,self.max_value = self.__create_source(gene)

    self.gene = gene

    self.error_pretext = bokeh.models.widgets.PreText(text='',width=100,height=20)

  def __joy(self,category,data,scale=0.2):
    return list(zip([category]*len(data), scale*data))

  def __read_data(self,filename):
    # read the file containing densities and other stuff
    with open(filename,'r') as f:
      covariates,timepoints,genotypes,x,data = pickle.load(f)
    return covariates, timepoints, genotypes, x, data

  def __create_source(self,gene):
    source = {}
    for genotype in self.genotypes:
      source[genotype] = {}
      for covariate in self.covariates:
        source[genotype][covariate] = bokeh.models.ColumnDataSource({'x': self.x})

    # get the maximum value (for setting y-axis ranges)
    max_value = -numpy.inf
    for covariate in self.covariates:
      for genotype in self.genotypes:
        for timepoint in self.timepoints:
          max_value = max([max_value,self.data[gene][covariate][genotype][timepoint].max()])

    # loop over anatomical sections
    for covariate in self.covariates:
      # loop over genotypes
      for genotype in self.genotypes:
        # loop over time points
        for timepoint in self.timepoints:
          # update the data object
          y = self.__joy(timepoint,self.data[gene][covariate][genotype][timepoint],scale=1/max_value)
          source[genotype][covariate].add(y,timepoint)

    return source,max_value

  def __update_plot(self,attr,old,new):
    if new not in self.data.keys():
      self.error_pretext.text = 'gene not found!'
      return
    self.error_pretext.text = ''

    source,max_value = self.__create_source(new)

    # TODO: update the y-axis range
    for covariate in self.covariates:
      for genotype in self.genotypes:
        self.source[genotype][covariate].data = source[genotype][covariate].data
    self.max_value = max_value

  def plot(self,n_cols=4):
    palette = [colorcet.rainbow[idx] for idx in numpy.linspace(0,255,len(self.genotypes)).astype(int)]

    plots = []

    # an input thingie for user
    textinput_gene = bokeh.models.widgets.TextInput(value=self.gene,title='Gene:')
    textinput_gene.on_change('value',self.__update_plot)
    plots.append([textinput_gene,self.error_pretext])

    # loop over tissue categories
    for idx,covariate in enumerate(self.covariates,start=0):
      # four subplots per row
      if idx % n_cols == 0:
        subplots = []

      # initialize a bokeh figure
      s = bokeh.plotting.figure(y_range=self.timepoints,x_range=(-5,5),width=200,plot_height=250,title=COVARIATES[covariate],tools='pan,wheel_zoom,reset')

      # loop over genotypes
      for genotype,color in zip(self.genotypes,palette):
        # loop over time points
        for timepoint in self.timepoints:
          # plot the distribution
          # legend is only drawn in the first subplot
          if covariate == self.covariates[0] and timepoint == self.timepoints[0]:
            s.patch('x',timepoint,line_color='black',fill_color=color,alpha=0.6,source=self.source[genotype][covariate],line_width=0.5,legend=genotype)
          else:
            s.patch('x',timepoint,line_color='black',fill_color=color,alpha=0.6,source=self.source[genotype][covariate],line_width=0.5)
    
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
      if (idx+1) % n_cols == 0 or (idx+1) == len(self.covariates):
        plots.append(subplots)
  
    return plots

class CellTypesOnTissueSections:
  def __init__(self,cell_type,conditions,proportions_filename='../data/Deconvolution/P_optimize_10.tsv',metadata_filename='../data/Metadata/metadata_09252017.tsv'):
    self.cell_types,self.data = self.__read_data(conditions,proportions_filename,metadata_filename)
    self.source_spots,self.vmax = self.__create_source_spots(cell_type)
    self.source_image = self.__create_source_image()

    # initialize the colormapper and ticker
    color_mapper = bokeh.models.mappers.LinearColorMapper('Plasma256',low=0,high=self.vmax)
    ticker = bokeh.models.BasicTicker(base=2,mantissas=[1,5])

  def __read_data(self,conditions,proportions_filename,metadata_filename):
    # read the data
    results = pd.read_table(proportions_filename,header=[0,1],index_col=0)
    count_files = numpy.array(list(results.columns.levels[0]))

    cell_types = list(results.index)
    
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
    
    data = []

    for n,condition in enumerate(conditions,start=0):

      data.append([])

      for count_file_idx in genotypes[condition[0]][condition[1]]:
        tmp = metadata[metadata['Count file'] == count_files[count_file_idx].replace('Count_Tables/','')]
        gt,tp,sex = tmp['Genotype'].values[0],tmp['Timepoint'].values[0],tmp['Sex'].values[0]
        gt = gt[0:5]+(gt[5:] and '..')
        title = '%s %s %s (%s)'%(gt,tp,sex,count_files[count_file_idx].replace('_stdata_aligned_counts_IDs.txt','').replace('Count_Tables/',''))
      
        # extract coordinates
        coordinates = [map(float,foo.split('_')) for foo in results[count_files[count_file_idx]].columns]
        coordinates = numpy.array(coordinates)

        # read the annotations of the spots on the array
        if re.match(r".*_[0-9]+_stdata_aligned_counts_IDs.txt",count_files[count_file_idx]):
          covariates = pd.read_table('../data/Covariates/%s.tsv'%(re.sub(r"_[0-9]+_stdata_aligned_counts_IDs.txt",r"",count_files[count_file_idx].replace('Count_Tables/',''))),header=0,index_col=0)
        else:
          covariates = pd.read_table('../data/Covariates/%s.tsv'%(count_files[count_file_idx].replace('_stdata_aligned_counts_IDs.txt','').replace('Count_Tables/','')),header=0,index_col=0)
    
        # prepare the spot annotations for the hover thingie
        annotations = []
        for coordinate in results[count_files[count_file_idx]].columns:
          if coordinate in covariates.columns:
            annotations.append(covariates.index[covariates[coordinate] == 1][0])
          # this should not happen
          else:
            annotations.append('Undefined')
      
        # read the array he image
        image_filename = '../data/Images/'+count_files[count_file_idx].replace('_stdata_aligned_counts_IDs.txt','_small.jpg').replace('Count_Tables/','')
        tissue_image = Image.open('../data/Images/'+count_files[count_file_idx].replace('_stdata_aligned_counts_IDs.txt','_small.jpg').replace('Count_Tables/','')).convert('RGBA')
      
        # mapping between pixels and spot coordinates
        xdim,ydim = tissue_image.size
        pixel_dim = 194.0/((6200.0)/(xdim/0.05))
        pixel_dim = pixel_dim * 0.05
      
        # cell type proportions
        proportions = results[count_files[count_file_idx]].as_matrix()
      
        # map the spot coordinates to pixels
        x_coordinates = pixel_dim*(coordinates[:,0]-1)
        y_coordinates = ydim-pixel_dim*(coordinates[:,1]-1)
      
        # relevant data
        data[n].append({'title': title,'image_filename': image_filename,'annotations': annotations,'proportions': proportions,'coordinates': [x_coordinates,y_coordinates],'dims': [xdim,ydim]})
      
    return cell_types, data 

  def __create_source_spots(self,cell_type):
    # get the maximum values
    vmax = -numpy.inf

    # loop over conditions
    for condition_data in self.data:
      # loop over count files (arrays)
      for array_data in condition_data:
        vmax = numpy.max([array_data['proportions'][self.cell_types.index(cell_type),:].max(),vmax])
  
    source_spots = []
    for n,condition_data in enumerate(self.data,start=0):
      source_spots.append([])
      for array_data in condition_data:
        source_spots[n].append(bokeh.models.ColumnDataSource({'x': array_data['coordinates'][0], 'y': array_data['coordinates'][1], 'z': array_data['proportions'][self.cell_types.index(cell_type),:], 'annotation': array_data['annotations']}))
    return source_spots, vmax

  def __create_source_image(self):
    source_image = []
    for n,condition_data in enumerate(self.data,start=0):
      source_image.append([])
      for array_data in condition_data:
        source_image[n].append(bokeh.models.ColumnDataSource({'image': [array_data['image_filename']]}))
    return source_image

  def __update_plot(self,attr,old,new):
    source_spots,vmax = self.__create_source_spots(new)

    for n in range(0,len(source_spots)):
      for m in range(0,len(source_spots[n])):
        self.source_spots[n][m].data = source_spots[n][m].data

    self.color_mapper.high = vmax 

    self.vmax_value = vmax

  def plot(self,cell_type):
    plots = []
  
    # input thingies for user
    select_cell_type = bokeh.models.widgets.Select(title='Cell type:', value=cell_type, options=self.cell_types)
    select_cell_type.on_change('value',self.__update_plot)
    plots.append([select_cell_type])
  
    # loop over result files
    for n,condition_data in enumerate(self.data,start=0):
      subplots = []
      
      # loop over count files (arrays)
      for m,array_data in enumerate(condition_data,start=0):
    
        # initialize hover thingie with expressions and annotations
        hover = bokeh.models.HoverTool(tooltips=[('Proportion','@z'),('Annotation','@annotation')])
    
        # initialize a bokeh figure
        s = bokeh.plotting.figure(width=235,plot_height=250,title=array_data['title'],x_range=(0,array_data['dims'][0]),y_range=(0,array_data['dims'][1]),match_aspect=True,aspect_scale=1,tools=[bokeh.models.tools.PanTool(),bokeh.models.tools.WheelZoomTool(),bokeh.models.tools.ResetTool(),hover])
        s.toolbar.logo = None
        s.toolbar_location = None
        s.axis.visible = False
    
        # plot the h&e stained image of the array
        s.image_url(url='image',x=0,y=0,anchor='bottom_left',w=array_data['dims'][0],h=array_data['dims'][1],source=self.source_image[n][m])
    
        # plot the spots
        s.scatter(x='x',y='y',radius=5,fill_color={'field': 'z','transform': self.color_mapper},fill_alpha=0.8,line_color=None,source=self.source_spots[n][m])
  
        subplots.append(s)

      plots.append(subplots)
    
    # initialize an empty bokeh figure for colorbar
    s = bokeh.plotting.figure(width=250,plot_height=100,title=None,x_axis_location=None,y_axis_location=None,tools='pan,wheel_zoom,reset',min_border=0,outline_line_color=None)
    color_bar = bokeh.models.ColorBar(color_mapper=self.color_mapper,ticker=self.ticker,label_standoff=6,border_line_color=None,location=(0,0),major_tick_line_color='black',title='Expression',orientation='horizontal')
    s.add_layout(color_bar,'above')
    s.toolbar.logo = None
    s.toolbar_location = None
    
    plots.append([s])
    
    return plots

class CellTypeSignatures:
  def __init__(self,filename='../data/Deconvolution/M_optimize_10.tsv'):
    self.data = pd.read_table(filename,header=0,index_col=0)
    self.source = self.__create_source()

  def __create_source(self):
    source_dict = {'Gene': list(self.data.index)}
    for col in self.data.columns:
      source_dict[col] = self.data[col].values

    return bokeh.models.ColumnDataSource(source_dict)

  def plot(self):
    # columns of the table
    formatter = formatter=bokeh.models.widgets.NumberFormatter(format='0.0000')
    columns = [bokeh.models.widgets.TableColumn(field='Gene',title='Gene')]
    for col in list(self.data.columns):
      columns.append(bokeh.models.widgets.TableColumn(field=col,title=col,formatter=formatter))
  
    return bokeh.models.widgets.DataTable(source=self.source,columns=columns,width=800,height=600,row_headers=False)

class CellTypesinAnatomicalSections:
  def __init__(self,anatomical_section,cell_types,timepoints=['p30','p70','p100','p120'],genotypes=['WT','G93A'],proportions_filename='../data/Deconvolution/P_optimize_10.tsv',metadata_filename='../data/Metadata/metadata_09252017.tsv'):
    self.covariates,self.cell_types,self.data = self.__read_data(proportions_filename,metadata_filename)
    self.timepoints = timepoints
    self.genotypes = genotypes
    self.source = self.__create_source(anatomical_section,cell_types)

    self.select_as = bokeh.models.widgets.Select(title='Anatomical section:',value=anatomical_section, options=self.covariates)
    self.select_ct_1 = bokeh.models.widgets.Select(title='First cell type:',value=cell_types[0], options=self.cell_types)
    self.select_ct_2 = bokeh.models.widgets.Select(title='Second cell type:',value=cell_types[1], options=self.cell_types)

  def __parse_proportions_df(self,proportions_df,count_files,metadata):
    data = []
    for count_file_idx in range(0,len(count_files)):
      filename = count_files[count_file_idx]
  
      coordinates = list(proportions_df[filename].columns)
   
      proportions = proportions_df[filename].as_matrix()
   
      tmp = metadata['Count file'] == filename.replace('Count_Tables/','')
      name,gt,tp,sex = metadata[tmp]['Name'].values[0],metadata[tmp]['Genotype'].values[0],metadata[tmp]['Timepoint'].values[0],metadata[tmp]['Sex'].values[0]
   
      if re.match(r".*_[0-9]+_stdata_aligned_counts_IDs.txt",filename):
        covariates = pd.read_table('../data/Covariates/%s.tsv'%(re.sub(r"_[0-9]+_stdata_aligned_counts_IDs.txt",r"",filename.replace('Count_Tables/',''))),header=0,index_col=0)
      else:
        covariates = pd.read_table('../data/Covariates/%s.tsv'%(filename.replace('_stdata_aligned_counts_IDs.txt','').replace('Count_Tables/','')),header=0,index_col=0)
   
      for idx,coordinate in enumerate(coordinates):
        data.append([name,coordinate,gt,tp,sex,covariates.index[covariates[coordinate] == 1][0]]+list(proportions[:,idx]))
    
    return pd.DataFrame(data,columns=['Slide','Coordinate','Genotype','Time point','Sex','Category']+list(proportions_df.index))


  def __read_data(self,proportions_filename,metadata_filename):
    # read the data
    results = pd.read_table(proportions_filename,header=[0,1],index_col=0)
    count_files = numpy.array(list(results.columns.levels[0]))
    covariates = list(pd.read_table('../data/Covariates/CN68_E1.tsv',header=0,index_col=0).index)

    cell_types = list(results.index)
    
    # read the metadata file
    metadata = pd.read_table(metadata_filename,header=0,index_col=0)
  
    # let parse our dataframe
    data = self.__parse_proportions_df(results,count_files,metadata)

    return covariates,cell_types,data

  def __create_source(self,anatomical_section,cell_types):
    source = {}
  
    for timepoint in self.timepoints:
      source[timepoint] = {}
      for genotype in self.genotypes:
        tmp_x = self.data[(self.data['Genotype'] == genotype) & (self.data['Time point'] == timepoint) & (self.data['Category'] == anatomical_section)][cell_types[0]]
        tmp_y = self.data[(self.data['Genotype'] == genotype) & (self.data['Time point'] == timepoint) & (self.data['Category'] == anatomical_section)][cell_types[1]]

        source[timepoint][genotype] = bokeh.models.ColumnDataSource({'x': tmp_x, 'y': tmp_y})
    
    return source 

  def __update_plot(self,attr,old,new):
    source = self.__create_source(self.select_as.value,[self.select_ct_1.value,self.select_ct_2.value])

    for timepoint in self.timepoints:
      for genotype in self.genotypes:
        self.source[timepoint][genotype].data = source[timepoint][genotype].data
  
  def plot(self,anatomical_section,cell_types):
    palette = [colorcet.rainbow[idx] for idx in numpy.linspace(0,255,len(self.genotypes)).astype(int)]
    plots = []
  
    self.select_as.on_change('value',self.__update_plot)
    self.select_ct_1.on_change('value',self.__update_plot)
    self.select_ct_2.on_change('value',self.__update_plot)
  
    plots.append([self.select_as,self.select_ct_1,self.select_ct_2])

    subplots = []
    for timepoint in self.timepoints:
      # initialize a bokeh figure
      s = bokeh.plotting.figure(width=235,plot_height=250,title=timepoint,x_range=(0,1),y_range=(0,1),x_axis_label='First cell type',y_axis_label='Second cell type',match_aspect=True,aspect_scale=1,tools=[bokeh.models.tools.PanTool(),bokeh.models.tools.WheelZoomTool(),bokeh.models.tools.ResetTool()])
      s.toolbar.logo = None
      s.toolbar_location = None
    
      # plot the spots
      if timepoint == self.timepoints[0]:
        for color,genotype in zip(palette,self.genotypes):
          s.scatter(x='x',y='y',radius=0.005,fill_color=color,fill_alpha=0.5,line_color=None,source=self.source[timepoint][genotype],legend=genotype)
      else:
        for color,genotype in zip(palette,self.genotypes):
          s.scatter(x='x',y='y',radius=0.005,fill_color=color,fill_alpha=0.5,line_color=None,source=self.source[timepoint][genotype])
    
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
    gene = 'Gfap'
    expression_on_tissue_sections = ExpressionOnTissueSections(gene,conditions)
  
    #viz_data = get_expression_on_tissue_sections_data(gene,conditions)
    plots = expression_on_tissue_sections.plot()
    p = bokeh.layouts.gridplot(plots,merge_tools=True,toolbar_location='left')
    tab_list.append(bokeh.models.widgets.Panel(child=p,title='Expression on tissue sections'))
  
  if options.expression_coefficients:
    gene = 'Gfap'
  
    #viz_data = get_expression_coefficients_data()
    #plots = viz_expression_coefficients_data(viz_data,gene)
    expression_coefficients = ExpressionCoefficients(gene=gene)
    plots = expression_coefficients.plot()
    p2 = bokeh.layouts.gridplot(plots,merge_tools=True,toolbar_location='left')
    tab_list.append(bokeh.models.widgets.Panel(child=p2,title='Expression coefficients'))
  
  if options.cell_types_on_tissues:
    cell_type = 'Cell type #6' 
    conditions = []
    for genotype in ['WT','G93A']:
      for timepoint in ['p30','p70','p100','p120']:
        conditions.append([genotype,timepoint])
    #conditions = [['G93A','p100']]
  
    #viz_data = get_cell_types_on_tissue_sections_data(cell_type,conditions)
    #plots = viz_cell_types_on_tissue_sections_data(viz_data)
    cell_types_on_tissue_sections = CellTypesOnTissueSections(cell_type,conditions)
    plots = cell_types_on_tissue_sections.plot(cell_type)

    p3 = bokeh.layouts.gridplot(plots,merge_tools=True,toolbar_location='left')
    tab_list.append(bokeh.models.widgets.Panel(child=p3,title='Cell types on tissue sections'))
  
  if options.cell_type_signatures:
    #viz_data = get_cell_type_signatures_data()
    #data_table = viz_cell_type_signatures_data(viz_data)
    cell_type_signatures = CellTypeSignatures()
    plots = cell_type_signatures.plot()


    p4 = bokeh.layouts.widgetbox(plots)
    tab_list.append(bokeh.models.widgets.Panel(child=p4,title='Cell type signatures'))
  
  if options.cell_types_in_anatomical_sections:
    anatomical_section = 'Vent_Horn'
    cell_types = ['Cell type #5', 'Cell type #9']

    cell_types_in_anatomical_section = CellTypesinAnatomicalSections(anatomical_section,cell_types)
    plots = cell_types_in_anatomical_section.plot(anatomical_section,cell_types)
  
    #viz_data = get_cell_types_in_anatomical_sections_data(anatomical_section,cell_types,timepoints,genotypes)
    #plots = viz_cell_types_in_anatomical_sections_data(viz_data)
    p5 = bokeh.layouts.gridplot(plots,merge_tools=True,toolbar_location='left')
    tab_list.append(bokeh.models.widgets.Panel(child=p5,title='Cell types in anatomical sections'))
  
  # create our tab pane object from our tabs
  tabs = bokeh.models.widgets.Tabs(tabs=tab_list)
  
  # finally, let us display our tab pane object
  bokeh.plotting.show(tabs)
