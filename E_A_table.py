#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 10 08:58:53 

  What are some realistic expections regarding shales? Shales are the finest
in size minerials, with some organics, floated andblown offshore into seas and 
lakes. They are deposited in mechanically quiet areas. Since river mouths,
shorelines, tectonics, and weather change rapidly over million-year time
spans the rates of deposition and composition change rapidly in vertical
sections. Shales are largely clays (40-60%) and quartz=SiO2,(20-50%); see 
Wikipedia. They can be interbedded with sands and carbonates.
  Silica is very important beyond quartz, as the SiO4(-4) group is a major
part of clays. Furthermore, when melted crustal materal cools as it comes
to the surface,quartz has the lowest melting temperature, about 600 deg C,
and is thus the 'pore fluid' and coats other minerals. See Youtube,C V 
Shorey, "Bowen's Reaction Series". Si is the closest analog of carbon.
  The above are some reasons I favor quartz pressure solution as being the
factor behind 'E' in the compaction model. See Miyakawa and Kawabe,
Pressure Solution of quartz aggregates under low effective stress.. etc.,
2014.
  Two more points: The lower limit on forces necessary to cause pressure
solution aren't known. Mechanisms are speculation I think. Moon and sun
forces produce extremely tiny displacements of atomic dimensions, but 
are continuous and changing 360 deg every day everywhere over all time,
and the moon was much closer 300 million years ago. Not usually mentioned.   
   
E_A_table

The forward model of shale deposition requires and defines parameters
E, A, m and n.
Derived parameters E and A for deposition rates between 0.2 km/My and
10 km/My for all six "wells" are put in the 'E_A_table' here. The parameters
m and n are the same as in the original paper. Parameter A encompasses
reaction frequency, but also a well's mineral composition, mineral 
morphology, mineral relative distributions in the section,geometrical
relationship between pores and minerals, etc.. Six varieties of quartz
have different solubilities at 25 deg C, and perhaps different E and A for
pressure solution Source: http://www.quartzpage.de/gen_chem.html. Google 
"The Quartz Page".      

@author: ed
"""

from pandas import DataFrame
#import pandas as pd

title = 'Parameters E and A for different matrix deposition rates'
data = {'Example':['Akita','5.-146 My','m=.95','n=1.','',
                   'Macran1','2.6-66 My','m=.85','n=1.','',  \
                   'Makran2','2.6-66 My','m=.85','n=.95','', \
                   '"SuluSea"','0.-23 My','m=.9','n=.8','',\
                   'Oklahoma','254-323 My','m=.85','n=1.','',\
                   'Maracaibo','2.6-66 My','m=.9','n=.9',''    ], 
        
        'Km/My':[.2, .3, 1., 4., 10,\
                 .2, .3, 1., 4., 10,\
                 .2, .3, 1., 4., 10,\
                 .2, .3, 1., 4., 10,\
                 .2, .3, 1., 4., 10,\
                 .2, .3, 1., 4., 10   ],
            
        'E /KJ/m':[28.8,29.2,29.7,29.8,29.9,\
                   16.6,17.5,19.2,19.8,19.9,\
                   12.2,13.2,14.8,15.5,15.6,\
                   11.9,13.4,16.9,18.6,19.0,
                   32.8,32.9,33.1,33.1,33.1,\
                   19.8,20.4,21.4,21.7,21.8],
            
        'A /s *10^-12': [18, 21, 25, 27, 27,\
                         1.2,1.7,3.4,4.4,4.7,\
                         .4, .61,1.2,1.6, 1.7,\
                         .49,.96,4.0,8.0,9.3,\
                         69.,71.,75.,76.,77.,\
                         3.7,4.8,7.2,8.3,8.5 ],
        }
frame = DataFrame(data)

frame.to_csv(r'/home/ed/Desktop/Z_Code0/'+ 'E_A_table '+'.csv', index=False)
print(frame)
pass
#print(' THE END')
