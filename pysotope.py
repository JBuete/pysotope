# Pysotope v0.5
# This is a basic implementation of a Pythonic class for isotopes.
# This file contains the main outline of the isotopes definition and the
# functions to parse the AME-2016 table as well as the NUBASE table.
#
# New in v0.5 LaTeX formatting output and half-life uncertainties
# 
# New in v0.4 ISOMERS
# 
# 
# Jacob Buete, January 2019
import logging
import os
import re

# define the version number
__version__ = '0.5.0'

# let's get the real path to this file so we can always find the mass table
_filepath = os.path.realpath(__file__)[:-11]

# and set up the logging information
logging.basicConfig(format="[%(levelname)s]\x1b[0;33;40m %(message)s \x1b[0m")

_mu_amu = 931.77 * 1e-6

_si_prefixes = {'Y': 24,  # yotta
                'Z': 21,  # zetta
                'E': 18,  # exa
                'P': 15,  # peta
                'T': 12,  # tera
                'G': 9,  # giga
                'M': 6,  # mega
                'k': 3,  # kilo
                'm': -3,  # milli
                'u': -6,  # micro
                'n': -9,  # nano
                'p': -12,  # pico
                'f': -15,  # femto
                'a': -18,  # atto
                'z': -21,  # zepto
                'y': -24,  # yocto
                }

_time_units = {'s': 1,
               'm': 60,
               'h': 3600,
               'd': 84600,
               'y': 31557600,
               }

def _correct_formatting(isotope):
    """Corrects the formatting of a provided isotope."""
    pattern = r"([0-9]*)?-?([a-z]*)-?([0-9]*)?"

    values = re.findall(pattern, isotope.lower())[0]

    # since we don't know the formatting that people use
    # we need to account for A-Symb and Symb-A
    if values[2] == "":
        corrected = values[0] + values[1]
    else:
        corrected = values[2] + values[1]

    # there are conventions that we need to account for
    if corrected == 'a':
        corrected = '4he'
    elif corrected == 't':
        corrected = '3h'
    elif corrected == 'd':
        corrected = '2h'
    elif corrected == 'p':
        corrected = '1h'
    elif corrected == 'n':
        corrected = '1n'

    return corrected



def _parse_mass_table():
    """Parses the isotope information from the Atomic Mass Evaluation table."""
    
    # read in the file
    with open(_filepath + "files/mass16.txt", "r") as f:
        raw_data = f.read()
        
    # now we can parse it with the regex
    pattern = r"^0?\s+\-?\d+\s+\d+\s+(\d+)\s+(\d+)\s+(\w+)"  # Z, A, name
    pattern += r"\s+([\+\-epnxaITt\d]*)?"  # decay processes
    pattern += r"\s+(\-?[\d\.#]+)\s+(\-?[\d\.#]+)"  # mass excess + unc.
    pattern += r"\s+(\-?[\d\.#]+)\s+(\-?[\d\.#]+)"  # binding energy + unc
    pattern += r"\s+B\-\s+(\-?[\d\.#\*]+)\s+(\-?[\d\.#]+)?"  # beta decay energy + unc
    pattern += r"\s+(\d+\s+[\d\.#]+)\s+([\d\.#]+)"  # atomic mass + unc

    data = re.findall(pattern, raw_data, re.MULTILINE)
    
    return data


def _parse_nubase_isotopes():
    """Parses the nubase text file to extract the stability and g.s. spin level."""
    # first open the file
    with open(_filepath + "files/nubase.txt", "r") as f:
        raw_data = f.read()
    # define the pattern
    pattern = r"^0*(\d+)\s+0*(\d+)0\s+(\w+)"  # A, Z, name
    pattern += r"\s+(\-?[\d\.#]+)\s+(\-?[\d\.#]+)"  # mass excess and uncertainty
    pattern += r"\s+(?:(?:2p|EU|RQ|MD|IT|\*&|AD|&|RN|BD|Nm|\*)*)?\s+(stbl|p\-unst#?|[<>~]?\d+\.?\d*#?)?"  # halflife or stability
    pattern += r"\s*([YZPETGMkmunpfazy]?[mhdsy])?\s*([\d\.]+#?|[<>]?[\d\.]+[YZPETGMkmunpfazy]?[mhdsy])?"
    pattern += r"\s+([\d/\+\-\(\),#\*]+#?)?\s+"

    data = re.findall(pattern, raw_data, re.MULTILINE)

    return data


def _parse_nubase_isomers():
    """Parses the nubase text file to extract the isomeric states."""
    # first open the file
    with open(_filepath + "files/nubase.txt", "r") as f:
        raw_data = f.read()
    
    # define the pattern
    pattern = r"^0*(\d+)\s+0*(\d+)\dW?\s+\d+(\D{2})([mnpqr])"  # A, Z, name
    pattern += r"\s+(\-?[\d\.#]+)\s+(\-?[\d\.#]+)"  # mass excess and uncertainty
    pattern += r"\s{1,9}(\-?[\d\.#]+)?\s{1,7}(\-?[\d\.#]+)?"  # other values
    pattern += r"\s+(?:(?:2p|EU|RQ|MD|IT|\*&|AD|&|RN|BD|Nm|\*)*)?\s+(stbl|p\-unst#?|[<>~]?\d+\.?\d*#?)?"# the decay path? and stability
    pattern += r"\s+(Yy|Zy|Ey|Py|Ty|Gy|My|ky|y|d|m|s|ms|ns|us|ps|fs|as|zs|ys)?\s*([\d\.]+#?|[<>]?[\d\.]+[YZPGTGMkmunpfazy]?[mhdsy])?"
    pattern += r"\s+([\d/\+\-\(\),#\*]+#?)?\s+"

    data = re.findall(pattern, raw_data, re.MULTILINE)
    
    return data

def _parse_tabmas():
    """Parses the Tabmas2012 file to extract deformation parameters."""
    
    with open(_filepath + "files/tabmas2012.dat", "r") as f:
        raw_data = f.read()
        
        
    # define the pattern
    pattern = r"^" + r"\s+(\d+)"*3 + r"\s+(\-?\d\.\d+)"*8

    data = re.findall(pattern, raw_data, re.MULTILINE)
    
    return data

def _gen_dicts():
    """Generates the isotope dictionary from the parsed mass table."""
    # load the data sources
    ame_data = _parse_mass_table()
    nu_data = _parse_nubase_isotopes()
    nu_isomers = _parse_nubase_isomers()
    nu_deformation = _parse_tabmas()
    
    # now load elements into the dictionary
    # using the formatted name as the key
    _iso_dict = {}
    _name_dict = {}
    _isomer_name_dict = {}
    for isotope in ame_data:
        if isotope[0] == '113':
            _name_dict[(int(isotope[0]), int(isotope[1]))] = isotope[1] + "nh"
        elif isotope[0] == '115':
            _name_dict[(int(isotope[0]), int(isotope[1]))] = isotope[1] + "mc"
        elif isotope[0] == '117':
            _name_dict[(int(isotope[0]), int(isotope[1]))] = isotope[1] + "ts"
        elif isotope[0] == '118':
            _name_dict[(int(isotope[0]), int(isotope[1]))] = isotope[1] + "og"
        else:
            _name_dict[(int(isotope[0]), int(isotope[1]))] = isotope[1] + isotope[2].lower()


        # now we should add the properties for the isotope
        properties = {}
        properties['Z'] = int(isotope[0])
        properties['A'] = int(isotope[1])
        properties['decay_proc'] = isotope[3]
        properties['mass excess'] = isotope[4].replace(" ", "")
        properties['mass excess uncertainty'] = isotope[5].replace(" ", "")
        properties['binding energy'] = isotope[6].replace(" ", "")
        properties['binding energy uncertainty'] = isotope[7].replace(" ", "")
        properties['beta decay energy'] = isotope[8].replace(" ", "")
        properties['beta decay energy uncertainty'] = isotope[9].replace(" ", "")
        properties['mass'] = isotope[10].replace(" ", "")
        properties['mass uncertainty'] = isotope[11].replace(" ", "")
            
        _iso_dict[_name_dict[(int(isotope[0]), int(isotope[1]))]] = properties

    # now we can use _name_dict to add in the stabilty and half_lifes
    for isotope in nu_data:
        name = _name_dict[(int(isotope[1]), int(isotope[0]))]

        # make a new dictionary for these properties
        nu_properties = {}
        nu_properties['half life'] = isotope[5]
        nu_properties['half life units'] = isotope[6]
        nu_properties['half life uncertainty'] = isotope[7]
        nu_properties['spin'] = isotope[8]
        
        
        _iso_dict[name].update(nu_properties)

    # now add the isomers
    for isomer in nu_isomers:
        base_name = _name_dict[(int(isomer[1]), int(isomer[0]))]
        
        if 'isomers' not in _iso_dict[base_name].keys():
            _iso_dict[base_name]['isomers'] = {}
        
        iso_properties = {}
        iso_properties['Z'] = int(isomer[1])
        iso_properties['A'] = int(isomer[0])
        iso_properties['mass excess'] = isomer[4]
        iso_properties['mass excess uncertainty'] = isomer[5]
        iso_properties['excitation energy'] = isomer[6]
        iso_properties['excitation energy uncertainty'] = isomer[7]
        iso_properties['half life'] = isomer[8]
        iso_properties['half life units'] = isomer[9]
        iso_properties['half life uncertainty'] = isomer[10]
        iso_properties['spin'] = isomer[11]
        
        
        _iso_dict[base_name]['isomers'][isomer[3]] = iso_properties
        # also define a mapping from the isomeric names to the base isotopes
        _isomer_name_dict[isomer[0] + isomer[2].lower() + isomer[3]] = base_name
        
    # now the deformation parameters
    for isotope in nu_deformation:
        try:
            name = _name_dict[(int(isotope[0]), int(isotope[2]))]
        except KeyError:
            continue
            
        def_properties = {}
        
        def_properties['epsilon2'] = isotope[3]
        def_properties['epsilon3'] = isotope[4]
        def_properties['epsilon4'] = isotope[5]
        def_properties['epsilon6'] = isotope[6]
        def_properties['beta2'] = isotope[7]
        def_properties['beta3'] = isotope[8]
        def_properties['beta4'] = isotope[9]
        def_properties['beta6'] = isotope[10]
        
        _iso_dict[name].update(def_properties)

    return _iso_dict, _name_dict, _isomer_name_dict

class Cluster(object):
    
    def __init__(self, iso, num: int):
        """Initialises a cluster of the same kind of isotopes."""
        self.type = iso
        self.num = num
        
        
    def __add__(self, b):
        """Defines the addition of Clusters with Isotopes."""
        total = b
        
        for i in range(self.num):
            total += self.type
            
        return total
    
    def __sub__(self, b):
        """Defines the subtraction of Clusters with Isotopes."""
        total = b
        for i in range(self.num):
            total -= self.type
            
        return total
    
    def __repr__(self):
        """Defines the representation of Clusters."""
        return "{} x {}".format(self.num, self.type)
    
    def get_z(self):
        """Returns the charge of the kind of Isotope in the Cluster."""
        return self.type.get_z()
    
    def get_a(self):
        """Returns the mass number of the kind of Isotope in the Cluster."""
        return self.type.get_a()
    
    def get_mass(self):
        """Returns the mass of the kind of Isotope in the Cluster in amu."""
        return self.type.get_mass()
    

class Isotope(object):
    """A class for performing basic linear operations on isotopes.
    
    All the information for the Isotopes has been taken from the 2016 Atomic Mass
    Evaluation table, which is parsed whenever the module is imported.
    
    Isotopes can be initialised by either their symbol and mass number as,

    >>> target = Isotope("208pb")
    >>> target = Isotope("pb208")

    or by giving their charge and mass number,

    >>> target = Isotope(82, 208)
    >>> target = Isotope(208, 82)

    Isotope also provides the following internal functions:

    >>> target.get_mass()  # returns the mass in micro-amu
    207976651.918  
    >>> target.get_mass_uncertainty()  # returns the uncertainty in the mass in micro-amu
    1.231
    >>> target.get_binding_energy()  # returns the binding energy in keV
    7867.453
    >>> target.get_be_uncertainty()  # returns the uncertainty in the binding energy in keV
    0.006
    >>> target.get_mass_excess()  # returns the mass excess in keV
    -21748.598
    >>> target.get_me_uncertainty()  # returns the uncertainty in the mass excess in keV
    1.148
    >>> target.get_z()  # returns the charge number
    82
    >>> target.get_a()  # returns the mass number
    208
    >>> target.get_info()  # prints some basic information about the isotope
    Isotope: 208Pb
    Mass = 207.9766519180 +/- 0.0000012310 u
    M_ex = -21.748598 +/- 0.001148 MeV
    B/A  = 7.867453 +/- 0.000006 MeV
    S_p  = 8.005489 MeV
    S_2p = 15.385295 MeV
    S_n  = 7.370052 MeV
    S_2n = 14.109828 MeV
    """
 
    def __init__(self, z, a=0):
        """Initialises the isotope with its name and properties."""
        if isinstance(z, str) and a == 0:  # case where we have the symbol and a
            self.name = _correct_formatting(z)
        elif isinstance(z, int) and a != 0:  # case where we have the charge and mass
            if a >= z:  # the good way
                self.name = _name_dict[(z, a)]
            else:  # to satiate ECS
                self.name = _name_dict[(a, z)]
        else:
            raise ValueError("Invalue Arguments.")
        # see if we can get the properties out
        try:
            self.properties = _iso_dict[self.name]
            self._isomer = None
        except KeyError:
            # this could mean that it's an isomeric state
            try:
                # get the base isotope name from the isomer list
                base_name = _isomer_name_dict[self.name]
                self._isomer = self.name[-1]  # get the isomer idenifier
                self.name = base_name  # bring the default name back
                self._base_properties = _iso_dict[self.name]
                self.properties = _iso_dict[self.name]['isomers'][self._isomer]  # and get the properties
            except KeyError:
                raise ValueError("Invalid Isotope Name")
            
    def __str__(self):
        """Defines the string representation for an isotope."""
        pattern = "(\d+)(\D+)"
        
        bits = re.findall(pattern, self.name)[0]

        result = bits[0] + bits[1].capitalize()
        
        if self._isomer is not None:
            result += self._isomer
            
        return result
    
    
    def __repr__(self):
        """Defines the representation of an isotope."""
        return "Isotope: " + self.__str__()
       
    
    def __add__(self, b):
        """Defines the addition (fusion) of isotopes."""
        
        if type(b) == Isotope:
            target_z = self.get_z() + b.get_z()
            target_a = self.get_a() + b.get_a()

            try:
                product = _name_dict[(target_z, target_a)]
            except KeyError:
                raise ValueError("Fusion of {} and {} is not known.".format(a, b))
        elif type(b) == Cluster:
            return b.__add__(self)
        else:
            raise ValueError("Addition only defined for Isotopes and Clusters.")
        
        return Isotope(product)
 
    def __radd__(self, b):
        """Defines reverse addition os isotopes so that the sum operator works."""
        return self

    def __sub__(self, b):
        """Defines the subtraction (fission?) of isotopes."""
        
        if type(b) == Isotope:
            target_z = abs(self.get_z() - b.get_z())
            target_a = abs(self.get_a() - b.get_a())

            try:
                product = _name_dict[(target_z, target_a)]
            except KeyError:
                raise ValueError("Breakup of " + str(self) + " and " + str(b) + " is not known.")
        elif type(b) == Cluster:
            return b.__sub__(self)
        else:
            raise ValueError("Subtraction only defined for Isotopes and Clusters.")
                
        return Isotope(product) 
    
    def __rmul__(self, b: int):
        """Defines 'groups' of similar nuclei"""
  
        return Cluster(self, b)
    
    def __mul__(self, b: int):
        """Defines 'groups' of similar nuclei"""
        
        return Cluster(self, b)
    
    def __eq__(self, b):
        """Defines the equality of isotopes."""
        return self.get_z() == b.get_z() and self.get_a() == b.get_a()
    
    def __eq__(self, b):
        """Defines the equality of isotopes."""
        return self.get_z() == b.get_z() and self.get_a() == b.get_a()

    def __lt__(self, b):
        """Defines an ordering for the isotopes."""
        return self.get_z() < b.get_z() or (self.get_z() == b.get_z and self.get_a() < b.get_z())

    def __hash__(self):
        """Defines the hash process for isotopes."""
        return hash(repr(self))
    
    def get_z(self):
        """Returns the Z of the isotope."""
        return int(self.properties['Z'])
    
    def get_a(self):
        """Returns the A of the isotope."""
        return int(self.properties['A'])
    
    def get_r(self, r_0=1.16):
        """Returns the radius of the isotope."""

        return r_0 * float(self.properties['A'])**(1/3)

    def get_mass(self):
        """Return the mass of the isotope."""
        if self._isomer is not None:
            mass = self._base_properties['mass']
            if "#" in mass:
                # logging.warning("Atomic mass of " + str(self) + " is estimated.")
                mass = mass[:-1]
            mass = float(mass) + 1e-3 * float(self.properties['excitation energy']) / _mu_amu
        else:
            mass = self.properties['mass']
            if "#" in mass:
                # logging.warning("Atomic mass of " + str(self) + " is estimated.")
                mass = mass[:-1]
            mass = float(mass)
        
        return mass
    
    def get_mass_uncertainty(self):
        """Return the mass of the isotope."""
        if self._isomer is not None:
            mass = self._base_properties['mass uncertainty']
            if "#" in mass:
                logging.warning("Atomic mass of " + str(self) + " is estimated.")
                mass = mass[:-1]
            mass = (float(mass)**2 + (1e-3 * float(self.properties['excitation energy uncertainty']) / _mu_amu)**2)**0.5
        else:
            mass = self.properties['mass uncertainty']
            if "#" in mass:
                logging.warning("Atomic mass of " + str(self) + " is estimated.")
                mass = mass[:-1]
            mass = float(mass)
        
        return mass
    
    def get_excitation_energy(self):
        """Returns the excitation energy of the isomer or zero for the groundstate."""
        if self._isomer is None:
            return 0
        else:
            energy = self.properties['excitation energy']
            if "#" in energy:
                logging.warning("Excitation energy of " + str(self) + " is estimated.")
                energy = energy[:-1]
            
            return float(energy)
        
    def get_excitation_energy_uncertainty(self):
        """Returns the excitation energy uncertainty of the isomer or zero for the groundstate."""
        if self._isomer is None:
            return 0
        else:
            energy = self.properties['excitation energy uncertainty']
            if "#" in energy:
                logging.warning("Excitation energy uncertainty of " + str(self) + " is estimated.")
                energy = energy[:-1]

            return float(energy)
    
    def get_binding_energy(self):
        """Returns the binding energy."""
        if self._isomer is None:
            binding = self.properties['binding energy']
        else:
            binding = self._base_properties['binding energy']
            logging.warning("The binding energy for an isomer defaults to the groundstate value")
        if "#" in binding:
            logging.warning("Binding energy for " + str(self) + " is estimated.")
            binding = binding[:-1]
        
        return float(binding)
    
    def get_beta(self, n, epsilon=False):
        """Returns the beta_n deformation parameter for the nucleus."""

        #  we don't have deformation values for Z < 8 nuclides
        if self.get_z() < 8:
            return 0

        if epsilon:
            pad = "epsilon"
        else:
            pad = "beta"

        if n not in [2, 3, 4, 6]:
            raise KeyError(pad + str(n) + " deformation not known.")

        return float(self.properties[pad + str(n)])

    def get_mass_excess(self):
        """Returns the mass excess."""
        excess = self.properties['mass excess']
        
        if "#" in excess:
            logging.warning("Mass excess for " + str(self) + " is estimated.")
            excess = excess[:-1]
        
        return float(excess)
    
    def get_be_uncertainty(self):
        """Returns the uncertainty in the binding energy."""
        if self._isomer is None:
            binding = self.properties['binding energy uncertainty']
        else:
            binding = self._base_properties['binding energy uncertainty']
            logging.warning("The binding energy for an isomer defaults to the groundstate value")
        if "#" in binding:
            binding = binding[:-1]
        
        return float(binding)
    
    def get_me_uncertainty(self):
        """Returns the uncertainty in the mass excess."""
        excess = self.properties['mass excess uncertainty']
        if "#" in excess:
            excess = excess[:-1]
        
        return float(excess)

    def _parse_time_units(self, units):
        """Parses the units of the half_life and returns the appropriate number of seconds."""
        # first get the scale of the unit
        if len(units) == 2:
            scale = _si_prefixes[units[0]]
        else:
            scale = 0
        
        # and now the base unit
        base_unit = _time_units[units[-1]]

        return base_unit * 10**scale
        
    def get_half_life(self, h=False):
        """Returns the half_life of the isotope in seconds or 1e99 if the isotope is stable."""
        half_life = self.properties['half life']

        # first we should determine if the isotope is stable or not
        if half_life in ['stbl', 'p-stbl']:
            return 1e99  # the isotope is stable, give back a huge number

        if half_life in ['p-unst', 'p-unst#']:
            return 0  # the isotope is kinda unstable

        if half_life == '':
            logging.warning("{} has no available half_life information, returning 0.".format(self))
            return 0

        # before moving on we should see if the 'h' argument has been used,
        # this implements a 'human readable' half_life and returns a string instead
        # this should only be used when someone makes something and wants to see the 
        # half_life instead of using it for a calculation
        # I don't know how I feel about including this
        if h:
            return half_life + " " + self.properties['half life units']

        # by this point we should only have isotopes for which we actually have
        # numbers for the half_lifes?
        # we will first remove the < and > signs because in general the limit of the
        # half_life should be enough
        half_life = re.findall(r"\d+\.?\d*#?", half_life)[0]

        # also check if the half_life has been estimated
        if '#' in half_life:
            logging.warning("Warning: The half_life of {} is estimated.".format(self))
            half_life = half_life.replace('#', '')

        # now get the units that we should scale by
        units = self._parse_time_units(self.properties['half life units'])

        return float(half_life) * units

    def get_half_life_uncertainty(self, h=False):
        """Returns the uncertainty in the half_life of the isotope."""
        # first check that there is an associated uncertainty for the half-life
        if self.properties['half life uncertainty'] == "":
            return 0
        
        # now we should check for units being present
        unit_pattern = r"[YZPGTGMkmunpfazy]?[mhdsy]"
        contained_units = re.findall(unit_pattern, self.properties['half life uncertainty'])
        
        if contained_units == []:  # then there are no units in the number
            units = self.properties['half life units']
        else:
            units = contained_units[0]
            
        scale = self._parse_time_units(units)

        half_life = re.findall(r"\d+\.?\d*", self.properties['half life uncertainty'])[0]
        print(half_life) 

        
        # before moving on we should see if the 'h' argument has been used,
        # this implements a 'human readable' half_life and returns a string instead
        # this should only be used when someone makes something and wants to see the 
        # half_life instead of using it for a calculation
        # I don't know how I feel about including this
        if h:
            return half_life + " " + units
        
        return float(half_life) * scale

    def is_stable(self, threshold=None):
        """Returns true if the isotope is stable."""
        if self.properties['half life'] == 'stbl' or self.properties['half life'] == 'p-stbl':
            return True

        if threshold is not None:
            # now get the half_life
            half_life = self.get_half_life()

            # and convert the threshold if we need to
            if isinstance(threshold, str):
                # break the value into the the number and the time units
                number, units = re.findall(r"([\d\.]+)(\D+)", threshold.replace(" ", ""))[0]
                threshold = float(number) * self._parse_time_units(units)

            if half_life > threshold:
                return True
        
        return False

    def has_isomers(self):
        """Returns true if the isotope has any isomeric states."""
        return 'isomers' in self.properties.keys()
    
    def select_isomer(self, isomer=None):
        """Selects the isomeric state for the isotope if it exists."""
        if isomer is None:
            self._isomer = None
        elif self._isomer is None and isomer in self.properties['isomers'].keys():
            self._isomer = isomer
            self._base_properties = _iso_dict[self.name]
            self.properties = self._base_properties['isomers'][isomer]
        elif isomer in self._base_properties['isomers'].keys():
            self._isomer = isomer
            self._base_properties = _iso_dict[self.name]
            self.properties = self._base_properties['isomers'][isomer]
        else:
            raise KeyError("{} is not a known isomeric state of {}".format(isomer, self))
    
        return self
    
    def get_isomers(self, list=False):
        """Returns the known isomers of the isotope."""
        if self.has_isomers():
            return self.properties['isomers'].keys()
            
        return []
        
    def get_spin(self):
        """Returns the spin of the isotope."""
        if self.properties['spin'] != "":
            return self.properties['spin']
        else:
            logging.warning("Spin not defined for {}".format(self))
            return None
        
    def get_info(self):
        """Provides a print-out of the known information about the isotope."""
        # basically we want to print out the name, mass, binding energy, and first few separations
        print(self.__repr__())
        print("Mass = {:.10f} +/- {:.10f} u".format(1e-6 * self.get_mass(),
                                                  1e-6 * self.get_mass_uncertainty()))
        print("M_ex = {:.6f} +/- {:.6f} MeV".format(1e-3 * self.get_mass_excess(),
                                                   1e-3 * self.get_me_uncertainty()))
        if self._isomer is not None:
            print("E_ex = {:.6f} +/- {:.6f} keV".format(self.get_excitation_energy(),
                                                        self.get_excitation_energy_uncertainty()))
        print("B/A  = {:.6f} +/- {:.6f} MeV".format(1e-3 * self.get_binding_energy(),
                                                   1e-3 * self.get_be_uncertainty()))
        if self.is_stable():
            print("t_1/2 = Stable")
        else:
            print("t_1/2 = {} +/- {} {}".format(self.properties['half life'],
                                                self.properties['half life uncertainty'],
                                                self.properties['half life units']))
        if self.get_spin is not None:
            print("Spin = " + self.get_spin())
        
        try:
            p_loss = self - Isotope('p')
            print("S_p  = {:.6f} MeV".format(- Reaction("{} -> {} + p".format(str(self),
                                                                              str(p_loss))).get_qvalue()))
        except ValueError:
            pass
        
        try:
            p2_loss = self - 2*Isotope('p')
            print("S_2p = {:.6f} MeV".format(- Reaction("{} -> {} + p + p".format(str(self),
                                                                             str(p2_loss))).get_qvalue()))
        except ValueError:
            pass
        
        try:
            n_loss = self - Isotope('n')
            print("S_n  = {:.6f} MeV".format(- Reaction("{} -> {} + n".format(str(self),
                                                                              str(n_loss))).get_qvalue()))
        except ValueError:
            pass
        
        try:
            n2_loss = self - 2*Isotope('n')
            print("S_2n = {:.6f} MeV".format(- Reaction("{} -> {} + n + n".format(str(self),
                                                                             str(n2_loss))).get_qvalue()))
        except ValueError:
            pass
       
    def latex(self):
        """Returns the LaTeX markup for the given isotope."""
        # get the symbol
        symb = re.findall(r"\D+", self.name)[0]
        
        return r"${{}}^{{{:}}}${:}".format(self.get_a(), symb.capitalize())
            
class Reaction(object):
    
    def _parse_reaction(self, reaction):
        """Parses a given reaction into the products and the reactants."""
        # first split into products and reactants
        (reactants, products) = reaction.split("->")

        # parse the list of reactants and products removing spaces
        reactants = reactants.lower().replace(" ", "").split("+")
        products = products.lower().replace(" ", "").split("+")

        
        # convert everything into isotope objects
        r_list = []
        p_list = []

        for reactant in reactants:
            if "*" in reactant:  # this just allows for some syntactic sugar in the reaction definition
                num, iso = reactant.split("*")
                for i in range(int(num)):
                    r_list.append(Isotope(iso))
            else:
                r_list.append(Isotope(reactant))
                
        for product in products:
            if "*" in product:
                num, iso = product.split("*")
                for i in range(int(num)):
                    p_list.append(Isotope(iso))
            else:
                p_list.append(Isotope(product))

        return r_list, p_list
        
    def _determine_missing_component(self, reactants, products):
        """Determines the missing isotope from a given set of reactants and products.
        
        Returns the missing isotope, or None if the isotope doesn't exist.
        """
        # first we should determine the charge and nucleon number of the reactants
        charge = 0
        nucleon = 0
        for iso in reactants:
            charge += iso.get_z()
            nucleon += iso.get_a()

        # now that we we can determine the missing components
        for iso in products:
            charge -= iso.get_z()
            nucleon -= iso.get_a()

        try:
            link = _name_dict[(charge, nucleon)]
            return Isotope(link)
        except KeyError:
            logging.error("Reaction cannot be balanced.")
            return None
    
    def __init__(self, reaction, ebeam=None, e_cm=None):
        """Initalises a reaction to contain the reaction definition, and the list of reactants
        and isotopes.
        """
        self.reaction = reaction
        self.reactants, self.products = self._parse_reaction(reaction)
        self.ebeam = ebeam
        self.e_cm = e_cm

    def get_barrier(self):
        """Calculates the barrier using the Swiatecki empirical equation."""

        if len(self.reactants) < 2:
            raise ValueError("Reaction needs two reactants to have a barrier!")

        # get the first reactant
        r1 = self.reactants[0]
        r2 = self.reactants[1]
        
        # define the coulomb parameter
        coulomb_param = r1.get_z() * r2.get_z() / (r1.get_a()**(1/3) + r2.get_a()**(1/3))

        # now calculate the barrier
        barrier = 0.85247 * coulomb_param 
        barrier += 0.001361 * coulomb_param**2 
        barrier -= 0.00000223 * coulomb_param**3

        return barrier

    def get_excitation_energy(self):
        """Returns the excitation energy for the compound nucleus."""

        if len(self.products) > 1:
            raise ValueError("No single compound nucleus formed.")

        return self.get_qvalue() + self.get_e_cm()

    def get_qvalue(self):
        """Calculates the q-value for the reaction."""
        
        reactant_energy = 0
        for r in self.reactants:
            reactant_energy += r.get_mass()
            
        product_energy = 0
        for p in self.products:
            product_energy += p.get_mass()
            
        return (reactant_energy - product_energy) * _mu_amu
        
    def get_q_opt(self, ebeam=None):
        """Calculates the optimal q-value for the reaction at the given beam energy."""
        # first let's calculate the cenre of mass energy
        e_cm = self.get_e_cm(ebeam)

        # now let's work on the charge fraction bit of q_opt
        r_charge = 1
        for r in self.reactants:
            r_charge *= r.get_z()

        # the products will be difficult as we could have someone putting in a break up reaction
        # that someone probably being me... so we need to put all of the smaller break up
        # fragments into a single variable
        heavy_product = max(self.products, key=(lambda x: x.get_mass()))
        mid_isotope = []
        for p in self.products:
            if p != heavy_product:
                mid_isotope.append(p)

        mid_isotope = sum(mid_isotope)

        p_charge = mid_isotope.get_z() * heavy_product.get_z()

        return e_cm * (p_charge / r_charge - 1)
         
    def set_ebeam(self, ebeam):
        """Sets the beam energy for the reaction."""
        self.ebeam = ebeam
        self.e_cm = self.get_e_cm(ebeam)

    def set_e_cm(self, e_cm):
        """Sets the centre of mass energy for the reaction."""
        self.e_cm = e_cm

        # get the target and projectile
        target = max(self.reactants, key=(lambda x: x.get_mass()))
        projectile = min(self.reactants, key=(lambda x: x.get_mass()))
        

        self.ebeam = self.e_cm * (target.get_mass() + projectile.get_mass()) / target.get_mass()

    def get_ebeam(self):
        """Returns the beam energy (if set)."""

        if self.ebeam is not None:
            return self.ebeam

        raise ValueError("Reaction hasn't been specified with a beam energy.")

    def get_e_cm(self, ebeam=None):
        """Calculates the centre of mass energy for the reactants."""
        if ebeam is None and self.ebeam is None:
            raise ValueError("Reaction hasn't been specified with a beam energy.")

        target = max(self.reactants, key=(lambda x: x.get_mass()))
        projectile = min(self.reactants, key=(lambda x: x.get_mass()))

        return self.ebeam * target.get_mass() / (target.get_mass() + projectile.get_mass())

    def is_valid(self):
        """Determines if a given reaction is valid by comparing charge and nucleon conservation."""

        # first determine the charge and nucleon numbers for the reactants
        charge = 0
        nucleon = 0
        for r in self.reactants:
            charge += r.get_z()
            nucleon += r.get_a()

        # now subtract the same quantities for the products
        for p in self.products:
            charge -= p.get_z()
            nucleon -= p.get_a()

        return charge == 0 and nucleon == 0


    def latex(self, math_env=False):
        """Returns the LaTeX formatted reaction."""
        
        
        r_side = " + ".join([r.latex() for r in self.reactants])
        p_side = " + ".join([p.latex() for p in self.products])
        
        # check if the user wants it inside a maths environment
        if math_env:

            r_side = r_side.replace("$", "")
            r_side = re.sub(r"([A-Z][a-z]?)", r"\text{\1}", r_side)

            p_side = p_side.replace("$", "")
            p_side = re.sub(r"([A-Z][a-z]?)", r"\text{\1}", p_side)

            return r_side + " \\rightarrow " + p_side

        return r_side + " $\\rightarrow$ " + p_side

    def __str__(self):
        """Defines the string representation of the reaction."""
        
        r_side = " + ".join([str(r) for r in self.reactants])
        p_side = " + ".join([str(p) for p in self.products])
        

        return r_side + " -> " + p_side
    
    def __repr__(self):
        """Defines the representation of the reaction object."""
        return "Reaction : " + str(self)
    
    def __hash__(self):
        """Defines the hash process for reactions."""
        return hash(repr(self))
 
_iso_dict, _name_dict, _isomer_name_dict = _gen_dicts()

# define the range of mass numbers for each of the isotopes by charge
mass_range = {}
symbol_dict = {}

# define the pattern for getting the symbol out of the name
_pattern = r"(\D+)"

for _z, _a in sorted(_name_dict.keys()):
    _name_string = _name_dict[(_z, _a)]
    _symb = re.findall(_pattern, _name_string)[0]
    if _z not in mass_range.keys():
        mass_range[_z] = []
        mass_range[_symb] = []
        mass_range[_symb.capitalize()] = []
    
    mass_range[_z].append(_a)
    mass_range[_symb].append(_a)
    mass_range[_symb.capitalize()].append(_a)

    if _z not in symbol_dict.keys():
        symbol_dict[_z] = _symb.capitalize()
