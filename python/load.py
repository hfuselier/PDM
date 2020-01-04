#!/usr/bin/env python
# coding: utf-8

# In[1]:


import numpy as np


# In[4]:


file_name_list = ['DDsand.txt',
             'Yuubarishale.txt',
             'RedWildmoorsand.txt',
             'Bereasand.txt',
             'Apuliancalcarenite.txt',
             'Indianalime.txt',
             'Mizuhotrachyte.txt',
             'Shirahamasand.txt',
             'Taiwansilt.txt',
             'Coconinosand.txt',
             'Bentheimsand.txt',
             'Dunnvillesand.txt',
             'Dunnvillesand_wTT.txt']
name_list = ['Darley Dale Sandstone',
        'Yuubari Shale',
        'Red Wildmoor Sandstone',
        'Berea Sandstone',
        'Apulian Calcarenite',
        'Indiana Limestone',
        'Mizuho Trachyte',
        'Shirahama Sandstone',
        'Taiwan Siltstone',
        'Coconino Sandstone',
        'Bentheim Sandstone',
        'Dunnville Sandstone Lit',
        'Dunnville Sandstone']


# In[12]:


def load_data(idx):
    fname = file_name_list[idx-1]
    name = name_list[idx-1]
    data = np.loadtxt('data/{}'.format(fname))
    return data




