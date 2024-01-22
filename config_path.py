#    Data analysis in OfficeEVparkingLot
#    Filter, process and analyze EV data collected at Dutch office building parking lot
#    Statistical analysis of EV data at ASR facilities - GridShield project - developed by 
#    Leoni Winschermann, University of Twente, l.winschermann@utwente.nl
#    Nataly Bañol Arias, University of Twente, m.n.banolarias@utwente.nl
#    
#    Copyright (C) 2022 CAES and MOR Groups, University of Twente, Enschede, The Netherlands
#
#    This library is free software; you can redistribute it and/or
#    modify it under the terms of the GNU Lesser General Public
#    License as published by the Free Software Foundation; either
#    version 2.1 of the License, or (at your option) any later version.
#
#    This library is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#    Lesser General Public License for more details.
#
#    You should have received a copy of the GNU Lesser General Public
#    License along with this library; if not, write to the Free Software
#    Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301
#    USA

import os, sys

class Config:
    def initializePath():
        #specify os path
        #os.chdir("C:/Users/WinschermannL/OneDrive - University of Twente/Documenten/Gridshield/PriorityAggregatedSched")
        #os.chdir("C:/Users/WinschermannL/OneDrive - University of Twente/Documenten/Gridshield/ProvideHomeCommute")
        os.chdir("C:/Users/WinschermannL/OneDrive - University of Twente/Documenten/Gridshield/ValueOfInformation")
        #os.chdir('C:/Users/BanolAriasMN/OneDrive - University of Twente/mnba/projects/project-Griedshield/EV-ASR')
    
    def path(dir):
        if dir == "VoI":
            os.chdir("C:/Users/WinschermannL/OneDrive - University of Twente/Documenten/Gridshield/ValueOfInformation")
        elif dir == "VoI_code":
            sys.path.insert(0, 'C:/Users/WinschermannL/OneDrive - University of Twente/Documenten/Gridshield/ValueOfInformation/commercial-building-ev-profiles-main')
        elif dir == "FOCS_code": 
            sys.path.insert(0, 'C:/Users/WinschermannL/OneDrive - University of Twente/Documenten/Gridshield/Criticalintervals/FOCS_code')
        elif dir == "FOCS_data": 
            sys.path.insert(0, 'C:/Users/WinschermannL/OneDrive - University of Twente/Documenten/Gridshield/Criticalintervals/FOCS_code/data')



