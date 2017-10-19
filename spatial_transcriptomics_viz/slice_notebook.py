"""
Helpers for displaying slice report data in jupyter notebooks.
"""

from spatial_transcriptomics_viz import slice_report
from bokeh.io import output_notebook
import bokeh.layouts
import bokeh.plotting
import bokeh.io
import ipywidgets as widgets
from IPython.display import display

output_notebook()

class GeneExpressionView(object):

    def __init__(self, gene='Gfap', conditions=[("WT", "p030"), ('G93A', 'p100')]):
        self.gene = gene
        self.conditions = conditions

    def show(self, E=None):
        if E is None:
            E = slice_report.ExpressionOnTissueSections(self.gene, self.conditions)
        plots = E.plot()
        p = bokeh.layouts.gridplot(plots,merge_tools=True,toolbar_location='left')
        bokeh.io.show(p, notebook_handle=True)

def showGeneExpressions(gene='Gfap', conditions=[("WT", "p030"), ('G93A', 'p100')]):
    G = GeneExpressionView(gene, conditions)
    G.show()

class ExpressionSelector(object):

    "Interactive Controls for selecting conditions and expressions in jupyter"

    def __init__(self):
        genes = self.get_genes()
        timepoints = self.get_timepoints()
        genotypes = self.get_genotypes()
        self.select_gene = sg = widgets.Select(
            options=genes,
            value=genes[0],
            description="genes",
            disabled=False
        )
        self.select_timepoints = st = widgets.SelectMultiple(
            options=timepoints,
            value=timepoints[:1],
            description="time points",
            disabled=False
        )
        self.select_genotypes = sgt = widgets.SelectMultiple(
            options=genotypes,
            value=genotypes[:1],
            description="genotypes",
            disabled=False
        )
        showbutton = self.show_button = widgets.Button(description="show")
        showbutton.on_click(self.show_click)
        self.assembly = widgets.HBox(children=(sg, st, sgt, showbutton))
        display(self.assembly)

    def show_click(self, *ignored):
        gene = self.select_gene.value
        timepoints = self.select_timepoints.value
        genotypes = self.select_genotypes.value
        conditions = []
        for gt in genotypes:
            for tp in timepoints:
                conditions.append((gt, tp))
        showGeneExpressions(gene, conditions)

    def get_genes(self):
        # XXXX temp
        return ["Gfap"]

    def get_timepoints(self):
        # XXXX temp
        return ['p030','p070','p100','p120']

    def get_genotypes(self):
        # XXXX temp
        return ['WT','G93A']