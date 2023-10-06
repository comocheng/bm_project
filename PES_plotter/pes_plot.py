import os
import sys
from energydiagram import ED
import logging
import cantera as ct
print(ct.__file__)
import numpy as np
import re
############################################
#
#   Plots a potential energy surface 
#   (enthalpy vs rxn coordinate) for a 
#   given cti file mechanism
#   
#   uses https://github.com/giacomomarchioro/PyEnergyDiagrams
#
############################################

class pes_reaction_combine():
    """
    feed in a cantera reaction, get an object out of it that we can use for making a chart

    arguments: 
    reaction - a ct reaction object
    phase_gas - gas phase in mechanism (for looking up species)
    phase_surf - solid phase in mechanism (for looking up species)
    reverse - bool, True if we swap reactants and products, and change Ea to Ea from products

    properties:
    reactants - dict, reactant string as key (e.g. "CO2+H2O"), combined Hf as value
    products - dict, reactant string as key (e.g. "CO2+H2O"), combined Hf as value
    barrier - float, Ea for reaction as value
    equation - string, cantera reaction equation
    links - list of ids for connecting the diagram
            [reactant id, Ea id, product id]
    positions - int, list of positions for reactants, reactions, and products. 
            [reactant pos, Reaction Ea pos, product pos]
    """
    def __init__(
        self,
        reaction,
        phase_gas,
        phase_surf,
        reverse=False,
        ):

        self.equation = reaction.equation
        self.reactants  = {}
        self.products  = {}
        self.links = [-1,-1,-1]
        self.positions = [-1, -1, -1]


        if reverse:
            reactants_obj = reaction.products
            products_obj = reaction.reactants

            # flip reaction equation, e.g. A <=> B is now B <=> A
            str_orig = reaction.equation
            split_list = str_orig.split("<=>")
            str1 = split_list[1] + " <=> " + split_list[0]
            str1 = str1.strip() 
            self.equation = str1
            print("flipped equation: ", str_orig, self.equation)

        else:
            reactants_obj = reaction.reactants
            products_obj = reaction.products
            self.equation = reaction.equation

        # lookup each reactant, put in dictionary as 
        # {species_name : enthalpy at phase temp (in Ev) * stoich coeff}
        total_reac_enth = 0
        reac_str = ""
        for i in reactants_obj: 
            if i in phase_gas.species_names:
                phase = phase_gas
            elif i in phase_surf.species_names:
                phase = phase_surf
            else: 
                logging.error(f"species {i} cannot be found")

            reac_str += f"{i} " * int(reactants_obj[i])
            total_reac_enth += reactants_obj[i] * (phase.species(i).thermo.h(phase.T)/1000**2)/96.4915666370759

        reac_str = reac_str.strip()
        self.reactants[reac_str] = total_reac_enth 
        self.reactants[reac_str] = round(self.reactants[reac_str],3)

        total_prod_enth = 0
        prod_str = ""
        for i in products_obj: 
            if i in phase_gas.species_names:
                phase = phase_gas
            elif i in phase_surf.species_names:
                phase = phase_surf
            else: 
                logging.error(f"species {i} cannot be found")

            prod_str += f"{i} " * int(products_obj[i])
            total_prod_enth += products_obj[i] * (phase.species(i).thermo.h(phase.T)/1000**2)/96.4915666370759
        
        prod_str = prod_str.strip()
        self.products[prod_str] = total_prod_enth 
        self.products[prod_str] = round(self.products[prod_str],3)

        # reaction activation energy. need to add to 
        # reactant enthalpy to get barrier
        # if reversed, need to get barrier from the products
        # duplicate surface Blowers Masel reactions are not considered
        if reaction.rate.type == 'interface-Blowers-Masel' or reaction.rate.type == 'sticking-Blowers-Masel':
            for i in range(phase_surf.n_reactions):
                if phase_surf.reaction(i).equation == reaction.equation:
                    delta_H = phase_surf.delta_enthalpy[i]
                    activation_energy = get_Ea_from_E0_dH(reaction.rate.activation_energy, delta_H, w0=1e9) / 96.4915666370759 / 1e6
                    if reverse: 
                        self.barrier = activation_energy + total_prod_enth
                    else: 
                        self.barrier = activation_energy + total_reac_enth
        else:   
            if reverse: 
                self.barrier = (reaction.rate.activation_energy/1000**2)/96.4915666370759 + total_prod_enth
            else: 
                self.barrier = (reaction.rate.activation_energy/1000**2)/96.4915666370759 + total_reac_enth

        self.barrier = round(self.barrier, 3)

class pes_reaction():
    """
    feed in a cantera reaction, get an object out of it that we can use for making a chart

    arguments: 
    reaction - a ct reaction object
    phase_gas - gas phase in mechanism (for looking up species)
    phase_surf - solid phase in mechanism (for looking up species)
    reverse - bool, True if we swap reactants and products, and change Ea to Ea from products

    properties:
    reactants - dict, species names as keys, Hf as value
    products - dict, species names as keys, Hf as value
    barrier - float, Ea for reaction as value
    equation - string, cantera reaction equation
    links - list of ids for connecting the diagram
            [[reactant ids], Ea id, [product ids]]
    """
    def __init__(
        self,
        reaction,
        phase_gas,
        phase_surf,
        ):

        self.equation = reaction.equation
        self.reactants  = {}
        self.products  = {}

        # lookup each reactant, put in dictionary as 
        # {species_name : enthalpy at phase temp (in Ev) * stoich coeff}
        total_reac_enth = 0
        for i in reaction.reactants: 
            if i in phase_gas.species_names:
                phase = phase_gas
            elif i in phase_surf.species_names:
                phase = phase_surf
            else: 
                logging.error(f"species {i} cannot be found")
            
            self.reactants[i] = reaction.reactants[i] * (phase.species(i).thermo.h(phase.T)/1000**2)/96
            total_reac_enth += self.reactants[i]
            self.reactants[i] = round(self.reactants[i],3)

        for i in reaction.products: 
            if i in phase_gas.species_names:
                phase = phase_gas
            elif i in phase_surf.species_names:
                phase = phase_surf
            else: 
                logging.error(f"species {i} cannot be found")

            self.products[i] = reaction.products[i] * (phase.species(i).thermo.h(phase.T)/1000**2)/96
            self.products[i] = round(self.products[i],3)

        # reaction activation energy. need to add to 
        # reactant enthalpy to get barrier
        self.barrier = (reaction.rate.activation_energy/1000**2)/96 + total_reac_enth
        self.barrier = round(self.barrier, 3)

class pes_plot():
    """
    Plots a potential energy surface 
    (enthalpy vs rxn coordinate) for a 
    given cti file mechanism
    """
    def __init__(
        self,
        yaml_file,
        temp=700,
        press=1,
        ):
        """
        initialize model
        yaml_file = cti or yaml file for mechanism
        temp = temperature (K)
        press = pressure (atm)
        """

        # set initial temps & pressures
        self.temp = temp # kelvin
        self.pressure =  press * ct.one_atm  # Pascals

        # create thermo phases
        self.yaml_file = yaml_file
        self.gas = ct.Solution(yaml_file, "gas")
        self.surf = ct.Interface(yaml_file,"surface1", [self.gas])

        # initialize T and P for each phase
        self.gas.TP = self.temp, self.pressure
        self.surf.TP = self.temp, self.pressure

        # create reaction diagram object
        self.diagram = ED()

        # initialize the reaction object dictionary
        self.pes_rxn_dict = {}

    def get_h_ev(self, species, temp):
        """
        gets species enthalpy in eV. 
        species is a cantera Species object
        """
        h_eV = (species.thermo.h(temp)/1000**2)/96
        print(f'{species.name} enthalpy = {h_eV} eV')
        return h_eV

    def get_ea_ev(self, reaction):
        """
        gets reaction Ea in eV. 
        reaction is a cantera Reaction object
        """
        Ea_eV = (reaction.rate.activation_energy/1000**2)/96
        print(f'{reaction.equation} enthalpy = {Ea_eV} eV')
        return Ea_eV

    def find_reactions(self, species, temp):
        """
        find all reactions that involve a certain species or set of species.
        species is a species object
        pes_rxn_dict is a dictionary, reaction equation is the key, PES reaction object is the value
        """
        pes_rxn_dict = {}
        species_names = [i.name for i in species]
        print(species_names)
        # get combined H for species as the starting point for Ea
        
        for i,j in enumerate(self.gas.reactions()):
            if all(x in j.reactants.keys() for x in species_names):
                pes_obj = pes_reaction_combine(j, self.gas, self.surf)
                pes_rxn_dict[pes_obj.equation] = pes_obj

            # if we want to show the reverse reaction, specify that 
            # in call to pes_reaction_combine
            if all(x in j.products.keys() for x in species_names):
                pes_obj = pes_reaction_combine(j, self.gas, self.surf, reverse=True)
                pes_rxn_dict[pes_obj.equation] = pes_obj
                
        for i,j in enumerate(self.surf.reactions()):
            if all(x in j.reactants.keys() for x in species_names):
                pes_obj = pes_reaction_combine(j, self.gas, self.surf)
                pes_rxn_dict[pes_obj.equation] = pes_obj

            # if we want to show the reverse reaction, specify that 
            # in call to pes_reaction_combine
            if all(x in j.products.keys() for x in species_names):
                pes_obj = pes_reaction_combine(j, self.gas, self.surf, reverse=True)
                pes_rxn_dict[pes_obj.equation] = pes_obj
        
        # if no reactions are found, throw an error
        if len(pes_rxn_dict)==0:
            raise Exception(f"no reactions found with reactants {species}")
                
        return pes_rxn_dict
        

    def plot_pes_diagram(
        self, 
        species, 
        width=20, 
        height=40, 
        offset=None,
        dimension=None,
        space=None,
        combined=True,
        ):
        """
        plots a potential energy surface given an input for species.
        the "species" are the starting point for the mechanism. each successive 
        run will take an input of the species and get all reactions for that pair. 

        inputs:
        species - (str or [strs]) matching starting reactant species name in cantera mechanism.
        width - (float) matplotlib plot width in inches
        height - (float) matplotlib plot height in inches
        offset - (float) vertical distance that energy level and upper/lower labels are spaced
        dimension - (float) width of platform used for energy level 
        space - (float) distance between bars for energy levels
        combined - (bool) if true combine all reactants to a single energy level. do the same for products.
        """

        # get a list of all reactions containing the two species identified
        species_obj = []
        for i in species:
            if i in self.gas.species_names:
                species_obj.append(self.gas.species(i))
            elif i in self.surf.species_names:
                species_obj.append(self.surf.species(i))
            else:
                print(f'species {i} not found!')

        rxns = self.find_reactions(species_obj, self.temp)
        self.pes_rxn_dict.update(rxns)

        link_num = 0
        for i,j in self.pes_rxn_dict.items():
            # generate a pes plot for each pes_reaction reactant
            for k,l in j.reactants.items():
                reac = k
                H_r = l

                # make a new energy level
                self.diagram.add_level(H_r, k, 0)
                j.positions[0] = 0
                    
            j.links[0] = link_num
            link_num+=1

        for i,j in self.pes_rxn_dict.items():
            # plot rxn Ea. for it to show up between species, should be here
            rxn_eq = j.equation
            rxn_Ea = j.barrier

            # make a new energy level
            self.diagram.add_level(rxn_Ea, rxn_eq, 1)
            j.positions[1] = 1
            
            j.links[1] = link_num
            link_num+=1

        for i,j in self.pes_rxn_dict.items():
            # generate a pes plot for each pes_reaction product

            for k,l in j.products.items():
                prod = k
                H_p = l

                # make a new energy level
                self.diagram.add_level(H_p, prod,2)
                j.positions[2] = 2

            # add link id for line drawing
            j.links[2] = link_num
            link_num+=1

        
        for i in self.pes_rxn_dict.values():
            # get connections between each reac - Ea and each Ea - product
            self.diagram.add_link(i.links[0],i.links[1])
            self.diagram.add_link(i.links[1],i.links[2])

        # optional arguements 
        if space: 
            self.diagram.space = space
        if offset:
            self.diagram.offset = offset
        if dimension:
            self.diagram.dimension = dimension
        
        self.diagram.ax.set_yticks(labelsize=36)
        self.diagram.plot(show_IDs=True, ylabel="Energy, $\mathregular{eV}$", width=width, height=height)

    def add_next_reaction(
        self, 
        species, 
        width, 
        height, 
        offset=None,
        dimension=None,
        space=None,
        ):
        """
        adds the next reaction to the plot. 
        
        species is the product species selected for the next step
        """

        # get a list of all reactions containing the species identified
        species_obj = []
        for i in species:
            if i in self.gas.species_names:
                species_obj.append(self.gas.species(i))
            elif i in self.surf.species_names:
                species_obj.append(self.surf.species(i))

        rxns = self.find_reactions(species_obj, self.temp)
        self.pes_rxn_dict.update(rxns)

        # ED.position can be assigned to an integer (1,2,3,4, etc) 
        # so we do not need to use "l"
        # get starting position (whatever the last energy level position was)
        starting_pos = max(self.diagram.positions)

        link_num = len(self.diagram.data)

        for i,j in rxns.items():
            # generate a pes plot for each pes_reaction
            for k,l in self.pes_rxn_dict[i].reactants.items():
                reac = k
                H_r = l

                # make a new energy level
                self.diagram.add_level(H_r, k, starting_pos)
                self.pes_rxn_dict[i].positions[0] = starting_pos
                    
            self.pes_rxn_dict[i].links[0] = link_num
            link_num+=1

        for i,j in rxns.items():   
            # plot rxn Ea. for it to show up between species, should be here
            rxn_eq = self.pes_rxn_dict[i].equation
            rxn_Ea = self.pes_rxn_dict[i].barrier

            # make a new energy level
            self.diagram.add_level(rxn_Ea, rxn_eq, starting_pos + 1)
            self.pes_rxn_dict[i].positions[1] = starting_pos + 1
            
            self.pes_rxn_dict[i].links[1] = link_num
            link_num+=1

        for i,j in rxns.items():  
            # generate a pes plot for each pes_reaction 
            for k,l in self.pes_rxn_dict[i].products.items():
                prod = k
                H_p = l

                # make a new energy level
                self.diagram.add_level(H_p, prod,starting_pos + 2)
                self.pes_rxn_dict[i].positions[2] = starting_pos + 2

            # add link id for line drawing
            self.pes_rxn_dict[i].links[2] = link_num
            link_num+=1

        
        for i,j in rxns.items():
            # get connections between each reac - Ea and each Ea-product
            self.diagram.add_link(self.pes_rxn_dict[i].links[0],self.pes_rxn_dict[i].links[1])
            self.diagram.add_link(self.pes_rxn_dict[i].links[1],self.pes_rxn_dict[i].links[2])

        # optional arguements 
        if space: 
            self.diagram.space = space
        if offset:
            self.diagram.offset = offset
        if dimension:
            self.diagram.dimension = dimension
        
        self.diagram.plot(show_IDs=True, ylabel="Energy / $eV$", width=width, height=height)



    def _redraw(
        self,
        pes_rxn_dict,
        width=20, 
        height=40,
        offset=None,
        dimension=None,
        space=None,
        ):

        """ redraw after a trim"""
        
        # create new reaction diagram object 
        # (is there a more efficient way to erase the old one?)
        self.diagram = ED()

        link_num = 0
        for i,j in pes_rxn_dict.items():

            # generate a pes plot for each pes_reaction reactant
            for k,l in j.reactants.items():
                reac = k
                H_r = l

                # make a new energy level
                self.diagram.add_level(H_r, k, j.positions[0])
                
            rxn_eq = j.equation
            rxn_Ea = j.barrier

            # make a new energy level
            self.diagram.add_level(rxn_Ea, rxn_eq, j.positions[1])

            # generate a pes plot for each pes_reaction product
            for k,l in j.products.items():
                prod = k
                H_p = l

                # make a new energy level
                self.diagram.add_level(H_p, prod, j.positions[2])

        self.diagram.create_data()

        # go through reaction dictionary and match reactant and product entries
        # match only adjacent ones
        for i,j in pes_rxn_dict.items():
            for m,n in enumerate(self.diagram.data):

                # there has to be a better way to do this ".keys())[0]"" nonsense
                if str(list(j.reactants.keys())[0]) == n[2] and j.positions[0] == n[1]:
                    j.links[0] = m
                    position = n[1]
                
                if j.equation == n[2] and j.positions[1] == n[1]:
                    j.links[1] = m

                if str(list(j.products.keys())[0]) == n[2] and j.positions[2] == n[1]:
                    j.links[2] = m


        # need to figure out how to add links
        for i in pes_rxn_dict.values():
            # get connections between each reac - Ea and each Ea - product
            self.diagram.add_link(i.links[0],i.links[1])
            self.diagram.add_link(i.links[1],i.links[2])

        # optional arguements 
        if space: 
            self.diagram.space = space
        if offset:
            self.diagram.offset = offset
        if dimension:
            self.diagram.dimension = dimension

        self.diagram.plot(show_IDs=True, ylabel="Energy / $eV$", width=width, height=height)

    def trim(
        self, 
        reacs, 
        width=20, 
        height=40,
        offset=None,
        dimension=None,
        space=None,
        ):
        """
        trims the specified reactions and their species from plot. 

        updates the pes_rxn_object to remove reactions we do not want. 
        runs through diagram.data and pes_rxn_object to update all of the links
        
        reac is the reaction to be trimmed (will remove reactants + products)
        """
        
        for key in list(self.pes_rxn_dict.keys()):
            if key in reacs: 
                del self.pes_rxn_dict[key]
        
        self._redraw(self.pes_rxn_dict)
        
    def plot_rxn_path(
        self, 
        rxns,
        width=20, 
        height=40,
        offset=None,
        dimension=None,
        space=None,
        ):
        """
        plots the energy path of a given surface reaction list
        make sure the reactions are in the order you want
        
        Inputs:
        rxns - a list of tuples cosists of reaction index and bool value
        """
        for i in rxns:
            rxn = pes_reaction_combine(self.surf.reaction(i[0]), self.gas, self.surf, i[1])
            
            # add level for reactants
            self.diagram.add_level(list(rxn.reactants.values())[0], list(rxn.reactants.keys())[0])
            
            # add level for barrier    
            rxn_eq = rxn.equation
            rxn_Ea = rxn.barrier
            self.diagram.add_level(rxn_Ea, rxn_eq)    
            
            # add level for products
            self.diagram.add_level(list(rxn.products.values())[0], list(rxn.products.keys())[0])
            
            if space: 
                self.diagram.space = space
            if offset:
                self.diagram.offset = offset
            if dimension:
                self.diagram.dimension = dimension
        
        for i in range(3*len(rxns)):
            if i != 3 * len(rxns) - 1:
                self.diagram.add_link(i, i+1)

        self.diagram.plot(show_IDs=True, ylabel="Energy / $eV$", width=width, height=height, fontsize=36)
        
def plot_multi_mech(
    mechs, 
    rxns,
    width=20, 
    height=40,
    offset=None,
    dimension=None,
    space=None,
    level_linewidth=1,
    link_linewidth=1
    ):
    """
    plot the same reaction path for multiple mechanisms
    
    Inputs:
    mechs - a list of tuples consists with Cantera gas and surface objects
    rxns - a list of tuples consists of reaction index and bool value
    """
    diagram = ED(level_linewidth=level_linewidth)
    diagram.color_bottom_text = 'black'
    level_ids = np.zeros((len(rxns), len(mechs), 3), dtype=np.int)
    start_id = -1
    colors = 'ybgrk'
    for r, i in enumerate(rxns):
        for j, mech in enumerate(mechs):
            # start_id += 1
            # level_ids[r][j][0] = start_id
            rxn = pes_reaction_combine(mech[1].reaction(i[0]), mech[0], mech[1], i[1])
            reactants = eliminate_index(list(rxn.reactants.keys())[0].split(' '))
            reactants_on_pt = eliminate_index(list(rxn.reactants.keys())[0].split(' '), Pt=True)
            reactants = '\n'.join(reactants)
            reactants_on_pt = '\n'.join(reactants_on_pt)
            if r != 0:
                rxn_last = pes_reaction_combine(mech[1].reaction(rxns[r-1][0]), mech[0], mech[1], rxns[r-1][1]) 
                # if list(rxn_last.products.keys())[0] != list(rxn.reactants.keys())[0]:
                #     # add level for reactants
                #     start_id += 1
                #     level_ids[r][j][0] = start_id
                #     if j != 0:
                #         diagram.add_level(list(rxn.reactants.values())[0], reactants, 'l', color=colors[j], top_text='')
                #     else:
                #         diagram.add_level(list(rxn.reactants.values())[0], reactants, color=colors[j],top_text='')
                if list(rxn_last.products.keys())[0] == list(rxn.reactants.keys())[0]:
                    level_ids[r][j][0] = start_id
                else:
                    # add level for reactants
                    start_id += 1
                    level_ids[r][j][0] = start_id
                    if j != 0:
                        diagram.add_level(list(rxn.reactants.values())[0], reactants, 'l', color=colors[j], top_text='')
                    else:
                        # diagram.add_level(list(rxn.reactants.values())[0], '', color=colors[j],top_text=reactants_on_pt)
                        diagram.add_level(list(rxn.reactants.values())[0], '', color=colors[j],top_text=' ')
            else:
                start_id += 1
                level_ids[r][j][0] = start_id
                if j != 0:
                    diagram.add_level(list(rxn.reactants.values())[0], reactants, 'l', color=colors[j], top_text='')
                else:
                    # diagram.add_level(list(rxn.reactants.values())[0], '', color=colors[j], top_text=reactants_on_pt)
                    diagram.add_level(list(rxn.reactants.values())[0], '', color=colors[j], top_text=' ')      
        for j, mech in enumerate(mechs):
            start_id += 1
            level_ids[r][j][1] = start_id
            # add level for barrier  
            rxn = pes_reaction_combine(mech[1].reaction(i[0]), mech[0], mech[1], i[1])
            rxn_eq = rxn.equation
            rxn_Ea = rxn.barrier
            if j!= 0:
                # diagram.add_level(rxn_Ea, rxn_eq, 'l')
                diagram.add_level(rxn_Ea, '', 'l', color=colors[j],top_text='')
            else:
                # diagram.add_level(rxn_Ea, rxn_eq)
                diagram.add_level(rxn_Ea, '', color=colors[j], top_text='')
        for j, mech in enumerate(mechs):
            start_id += 1
            level_ids[r][j][2] = start_id
            rxn = pes_reaction_combine(mech[1].reaction(i[0]), mech[0], mech[1], i[1])
            print(rxn.products)
            products = eliminate_index(list(rxn.products.keys())[0].split(' '))
            products_on_pt = eliminate_index(list(rxn.products.keys())[0].split(' '), Pt=True)
            products = '\n'.join(products)
            products_on_pt = '\n'.join(products_on_pt)
            # add level for products
            if j!=0:
                diagram.add_level(list(rxn.products.values())[0], products, 'l', color=colors[j], top_text='')
            else:
                # diagram.add_level(list(rxn.products.values())[0], '', color=colors[j],top_text=products_on_pt)
                diagram.add_level(list(rxn.products.values())[0], '', color=colors[j],top_text=' ')
    # add links from reactants to products for each reaction
    # For CMPO
    level_ids = [[[ 0,  3,  6],
                  [ 1,  4,  7],
                  [ 2,  5,  8],],

                 [[ 9, 12, 15],
                  [10, 13, 16],
                  [11, 14, 17],],

                 [[15, 18, 21],
                  [16, 19, 22],
                  [17, 20, 23],]]
    # For CMPO-BMA
#     level_ids = [[[ 0,  3,  6],
#                   [ 1,  4,  7],
#                   [ 2,  5,  8],],

#                  [[9,  12, 15],
#                   [10, 13, 16],
#                   [11, 14, 17],],

#                  [[18, 21, 24],
#                   [19, 22, 25],
#                   [20, 23, 26],],

#                  [[27, 30, 33],
#                   [28, 31, 34],
#                   [29, 32, 35],],

#                  [[36, 39, 42],
#                   [37, 40, 43],
#                   [38, 41, 44],]]

    for id_comb in level_ids:
        diagram.add_link(id_comb[1][0], id_comb[1][1], color=colors[1], linewidth=link_linewidth)
        diagram.add_link(id_comb[1][1], id_comb[1][2], color=colors[1], linewidth=link_linewidth)
        diagram.add_link(id_comb[0][0], id_comb[0][1], color=colors[0], linewidth=link_linewidth)
        diagram.add_link(id_comb[0][1], id_comb[0][2], color=colors[0], linewidth=link_linewidth)
        diagram.add_link(id_comb[2][0], id_comb[2][1], color=colors[2], linewidth=link_linewidth)
        diagram.add_link(id_comb[2][1], id_comb[2][2], color=colors[2], linewidth=link_linewidth)
        for rxn_ids in id_comb:
            diagram.add_link(rxn_ids[0], rxn_ids[1])
            diagram.add_link(rxn_ids[1], rxn_ids[2])
    
    # add links between reactions
    i = 0
    while i < len(rxns):
        if i != len(rxns) - 1:
            for num in range(len(mechs)):
                diagram.add_link(level_ids[i][num][2], level_ids[i+1][num][0], color=colors[num], linewidth=link_linewidth)
        i += 1
            
    if space: 
        diagram.space = space
    if offset:
        diagram.offset = offset
    if dimension:
        diagram.dimension = dimension
        
    diagram.plot(show_IDs=False, ylabel="Energy / $eV$", width=width, height=height)

# This function calculates Ea with given E0 and w0
def get_Ea_from_E0_dH(E0, dHrxn, w0=1e9):
        """
        Return the activation energy in J/kmol corresponding to the given
        E0, enthalpy of reaction `dHrxn`, and w0, all in J/kmol.
        """
        E0 = float(E0)
        w0 = max(w0, 2*E0)
        assert w0 >= 2*E0, f"seems to work best when w0 > 2*Eo = {2*E0/1e3:.1f} kJ/mol, but w0={w0/1e3:.1f} kJ/mol"
        
        if E0 == dHrxn == 0:
            return 0
        if dHrxn < -4 * E0:
            return 0.0
        elif dHrxn > 4 * E0:
            return dHrxn
        else:
            Vp = 2 * w0 * (2 * w0 + 2 * E0) / (2 * w0 - 2 * E0)
            return (w0 + dHrxn / 2.0) * (Vp - 2 * w0 + dHrxn) ** 2 / (Vp ** 2 - (2 * w0) ** 2 + dHrxn ** 2)

def eliminate_index(species, Pt=False):
    """
    This function eliminate the indices after species 
    species:a list of strings
    """
    new_sp_list = []
    for sp in species:
        new_r = []
        element_num = []
        for num in range(10):
            element_num.append(str(num))
        sp = sp.replace('_', '')
        sp = re.sub(r'\(\d+\)','', sp)
        sp = re.sub(r'(\d+)',r'_\1', sp)
        if Pt:
            sp = re.sub(r"X", "Pt", sp)
        else:
            sp = re.sub(r"X", "*", sp)
        sp = '$\mathregular{'+ sp + '}$'
        new_sp_list.append(sp) 
    return new_sp_list