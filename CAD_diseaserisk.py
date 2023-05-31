import pandas as pd 
import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt

CADdata = pd.read_csv('CAD-logisticpred-manual-plot.csv',names=['percentiles','disease_risk']).sort_values(by='percentiles')
percentile = list(CADdata['percentiles'])
percentile = [round(x)*0.01 for x in percentile]
disease_risk = list(CADdata['disease_risk'])
disease_risk = [x*0.01 for x in disease_risk]
