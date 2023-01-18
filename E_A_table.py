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
  When combining pressure solution and low-permeability models, the 
montmorillionite semi-permeable membrane phenomena should be mentioned as
impeding compaction. 
 
E_A_table

The forward model of shale deposition requires and defines parameters
E, A, m and n. Minimum t = 25 deg C.  
Derived parameters E and A for deposition rates between 0.2 km/My and
10 km/My for all six "wells" are put in the 'E_A_table' here. The parameters
m and n are the same as in the original paper. Parameter A encompasses
reaction frequency, but also a well's mineral composition, mineral 
morphology, mineral relative distributions in the section,geometrical
relationship between pores and minerals, etc.. Six varieties of quartz
have different solubilities at 25 deg C, and perhaps different E and A for
pressure solution Source: http://www.quartzpage.de/gen_chem.html. Google 
"The Quartz Page.  ""A Colmpilation of Rate Parameters.." J L Palandri,Y K 
Kharaka, USGSopen file report 2004-1068.   
    

@author: ed
"""

from pandas import DataFrame
#import pandas as pd
print( 'Deposition Rates With Parameters m,n,A and E: Forward Model ' )
title = 'Parameters E and A for different matrix deposition rates'
data = {'Example':['Akita','5.-146 My','m,n=.95,1.','t>=25C','','',     \
                   'Macran1','2.6-66 My','m,n=.85,1.','t>=25C','','',   \
                   'Makran2','2.6-66 My','m,n=.85,.95','t>=25C','','',  \
                   '"SuluSea"','0.-23 My','m,n=.9,.8','t>=25C','','',   \
                   'Oklahoma','254-323 My','m,n=.85,1.','t>=25C','','', \
                   'Maracaibo','2.6-66 My','m,n=.9,.9','t>=25C','',''   ], 
        
        '   Km/My':[.2, .3, 1., 4, 20,'',\
                 .2, .3, 1., 4, 20,'', \
                 .2, .3, 1., 4, 20,'', \
                 .2, .3, 1., 4, 20,'', \
                 .2, .3, 1., 4, 20,'',\
                 .2, .3, 1., 4, 20 ,'' ],
            
        '   E kJ/m':[30.8,31.1,31.6,31.8,31.8,'',\
                   17.2,18.1,19.7,20.3,20.4,'',\
                   14.4,15.4,17.0,17.5,17.7,'',\
                   12.3,13.9,17.4,19.0,19.5,'',\
                   31.7,31.8,31.9,31.9,31.9,  '',\
                   18.4,19.1,20.0,20.3,20.4,'' ],
            
        '    A/s*10^-12': [43,49,59,62,64,      '',\
                         1.4,2.1,4.0,5.1,5.4, '',\
                         .86,1.3,2.5,3.2,3.4, '',\
                         .55,1.1,4.6,8.9,11.0,'',\
                         48.,49.,52.,52.,53., '',\
                         2.4,3.1,4.5,5.2,5.3, '' ],
        }
frame = DataFrame(data)

frame.to_csv(r'/home/ed/Desktop/Z_Code0/'+ 'E_A_table '+'.csv', index=False)
print(frame)
pass
#print(' THE END')
