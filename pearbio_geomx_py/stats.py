
import pandas as pd
import anndata as ad
import seaborn as sns
import scanpy as sc
import statsmodels as sm
import patsy

from matplotlib.figure import Figure
from statsmodels.genmod.generalized_linear_model import GLM
from statsmodels.genmod.families import Poisson

adata_path = "/home/nourdine/dev/pearbio-geomx-nf/output/00_nanopore_data/GX0000.h5ad"
width = 8
height = 8
img_basepath = "test"

adata = ad.read_h5ad(adata_path)

#adata.obs["n_genes"] = np.apply_along_axis(
#a = np.apply_along_axis(
#    func1d = lambda x: np.nonzero(x),
#    axis = 0,
#    arr = adata.X
#    )

#sc.pp.calculate_qc_metrics(adata, inplace=True)

factors = [
    "Type",
    "Infiltration",
    "Patient"
]

variables = [
    "AlignedReads",
    "AOISurfaceArea",
    "AOINucleiCount"
]

#fig = Figure(figsize=(width,height))
#ax = fig.add_subplot(111)
#sns.violinplot(data=adata.to_df().T, ax=ax)
#fig.savefig(f"{img_basepath}.pdf")


#"""
#Visualizing mean-variance relationship and overdispersion.
#"""
#m_bar = np.apply_along_axis(
#    func1d = np.mean,
#    axis = 0,
#    arr = adata.X
#)
#v_bar = np.apply_along_axis(
#    func1d = np.var,
#    axis = 0,
#    arr = adata.X
#)
#z = np.log(np.linspace(0, np.max(m_bar), 1000))
#fig = Figure(figsize = (width,height))
#ax = fig.add_subplot(111)
#sns.scatterplot(x = np.log(m_bar), y = np.log(v_bar), ax = ax)
#ax.plot(z, z)
#fig.savefig(f"{img_basepath}.pdf")


"""
Visualizing count distribution.
"""
fig = Figure(figsize=(width, height))
ax = fig.add_subplot(111)
sns.histplot(data = adata.X[0], log_scale = True, ax = ax)
fig.savefig(f"{img_basepath}.pdf")


"""
Violinplot
"""
fig = Figure(figsize=(width, height))
ax = fig.add_subplot(111)
sns.violinplot(data = np.log(adata.to_df().T), ax = ax)
fig.savefig(f"{img_basepath}.pdf")

predictors = [
    {
        "name": "Type",
        "type": "categorical",
        "levels": ["Normal", "Tumor"]
    },
    {
        "name": "Infiltration",
        "type": "categorical",
        "levels": ["Non-infiltrated", "Infiltrated"]
    }
]


gene = "A1BG"
expr = adata[:,gene]

y = np.apply_along_axis(
    func1d = lambda x: np.int64(x[0]),
    axis = 1,
    arr = expr.X
)

formula = "y ~ Type * Infiltration"
names = ["Intercept", "Type", "Infiltration", "Type:Infiltration"]


predictors_names = list(map(lambda x: x["name"], predictors))
#X = adata.obs[predictors_names].reset_index(drop = True)

data = {"y": y}
for predictor in predictors:

    name = predictor["name"]
    levels = predictor["levels"]

    x = adata.obs[name]
    data[name] = pd.Categorical(x, categories=levels)

data = pd.DataFrame(data)
design_matrix = patsy.dmatrices(formula, data)
X = np.asarray(design_matrix[1])
X = pd.DataFrame(data = X, columns = names)

model = GLM(y, X, family = Poisson())
model = model.fit()





#pair_grid = sns.pairplot(adata.obs[variables])
#pair_grid.savefig(f"{img_basepath}.pdf")
#
#for factor in factors:
#    pair_grid = sns.pairplot(
#        adata.obs[ variables + [factor] ],
#        hue=factor
#    )
#    pair_grid.savefig(f"{img_basepath}.{factor}.pdf")
#    pair_grid.savefig(f"{img_basepath}.{factor}.png")

