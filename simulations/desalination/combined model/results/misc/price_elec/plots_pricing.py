#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec 13 23:46:12 2023

@author: elizabethkeith
"""

import seaborn as sns
import pandas as pd
import plotly.io as pio

dat = pd.read_csv('price_schedule.csv')
fig = sns.lineplot(data=dat, x='Hour', y='LMP', hue='Season')


pio.write_image(fig, "price_schedule.svg", width=1.5*300, height=0.75*300, scale=1)