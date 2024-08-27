import xml.etree.ElementTree as ET
import numpy as np
import subprocess
import copy
import os

class SkiGen(object):

    def __init__(self, cone_tau, cone_width=10, cone_dust_type='mrn77', cone_dust_minsize=None, cone_dust_maxsize=None, cone_dust_exponent=None):

        #Save the input parameters. 
        self.cone_tau = cone_tau
        self.cone_width = cone_width
        self.cone_dust_type = cone_dust_type
        self.cone_dust_minsize=cone_dust_minsize
        self.cone_dust_maxsize=cone_dust_maxsize
        self.cone_dust_exponent=cone_dust_exponent

        #The directory where this script lives. 
        skg_folder = os.path.dirname(os.path.realpath(__file__))

        #Start by reading the template.
        parser = ET.XMLParser(target=ET.TreeBuilder(insert_comments=True))
        self.tree = ET.parse(skg_folder+"/template.ski", parser)
        self.root = self.tree.getroot()

        #Separate the main components. 
        self.tor, self.cone = self.root.findall("./MonteCarloSimulation/mediumSystem/MediumSystem/media/GeometricMedium")
        self.incs = self.root.findall("./MonteCarloSimulation/instrumentSystem/InstrumentSystem/instruments")[0]

        #Set the cone dust properties.
        self.cone_dust = self.cone.find("./materialMix/ConfigurableDustMix/populations/GrainPopulation/sizeDistribution/PowerLawGrainSizeDistribution")

        if self.cone_dust_type is None: 
            if self.cone_dust_minsize is None or self.cone_dust_maxsize is None or self.cone_dust_exponent:
                print("You need to provide dust properties or a dust type for the cone.")
                return
            self.cone_dust.set("minSize", "{} micron".format(self.cone_dust_minsize))
            self.cone_dust.set("maxSize", "{} micron".format(self.cone_dust_maxsize))
            self.cone_dust.set("exponent", "{}".format(self.cone_dust_exponent))
        else:
            if self.cone_dust_type=='mrn77':
                self.cone_dust.set("minSize", "0.005 micron")
                self.cone_dust.set("maxSize", "0.25 micron")
                self.cone_dust.set("exponent", "3.5")
            else:
                print("Unrecognized dust type {}.".format(self.cone_dust_type))
                return

        #Set the cone dust optical depth. 
        self.cone.find("./normalization/OpticalDepthMaterialNormalization").set("opticalDepth","{}".format(2*self.cone_tau))

        #Set the base filename. 
        self.base_fname = "bHDPol"
        if self.cone_dust_type is None:
            self.base_fname += "_gs{}-{}_a{}".format(self.cone_dust_minsize, self.cone_dust_maxsize, self.cone_dust_exponent)
        else:
            self.base_fname += "_{}".format(self.cone_dust_type)

        return
    
    def write_files_parameter_grid(self, tor_oa_min, tor_oa_max, tor_doa, cone_oa_min, cone_oa_max, cone_doa, delta_eta=5, folder="scripts"):

        subprocess.call("mkdir {}".format(folder), shell=True)

        tor_oas  = np.arange(tor_oa_min , tor_oa_max +0.1*tor_doa , tor_doa )

        for tor_oa in tor_oas:

            #Set the torus opening angle. 
            self.tor.find("./geometry/TorusGeometry").set("openingAngle", "{} deg".format(90-tor_oa))

            #Set the inclinations. 
            eta_max = 90
            eta_min = tor_oa
            etas = np.arange(eta_min, eta_max+0.1*delta_eta, delta_eta)
            etas[0]+=0.2*delta_eta
            eta_temp = self.incs.find("FullInstrument")
            self.incs.remove(eta_temp)
            for eta in etas:
                eta_temp_use = copy.deepcopy(eta_temp)
                eta_temp_use.set("instrumentName","i{}".format(eta))
                eta_temp_use.set("inclination", "{} deg".format(eta))
                self.incs.append(eta_temp_use)

            cone_oa_max_use = np.min([cone_oa_max, tor_oa-5])
            cone_oas = np.arange(cone_oa_min, cone_oa_max_use+0.1*cone_doa, cone_doa)

            for cone_oa in cone_oas:

                #Set the cone geometry. 
                self.cone.find("./geometry/ConicalShellGeometry").set("maxAngle", "{} deg".format(90-cone_oa))
                self.cone.find("./geometry/ConicalShellGeometry").set("minAngle", "{} deg".format(90-cone_oa - self.cone_width))

                fname = folder+"/"+self.base_fname + "_tor_oa{}_con_oa{}-tauV{}.ski".format(tor_oa, cone_oa, self.cone_tau)

                self.tree.write(fname, xml_declaration=True, encoding='UTF-8')

        return