import xml.etree.ElementTree as ET
import numpy as np
import subprocess
import copy
import os

class SkiGen(object):

    def __init__(self, cone_tau, cone_width=10, cone_dust_type='mrn77', cone_dust_minsize=None, cone_dust_maxsize=None, cone_dust_exponent=None, cone_type="Full", SED_only=False):

        #Save the input parameters. 
        self.cone_tau = cone_tau
        self.cone_width = cone_width
        self.cone_dust_type = cone_dust_type
        self.cone_dust_minsize=cone_dust_minsize
        self.cone_dust_maxsize=cone_dust_maxsize
        self.cone_dust_exponent=cone_dust_exponent
        self.SED_only = SED_only
        self.cone_type = cone_type

        #The directory where this script lives. 
        skg_folder = os.path.dirname(os.path.realpath(__file__))

        #Set template and output instrument depending on whether we want SEDs only or the images as well. 
        if self.SED_only:
            self.output_instrument = "SEDInstrument"
            if cone_type=='Full':
                self.template   = "template_SEDInst.ski"
            elif cone_type=="Bottom":
                self.template   = "template_SEDInst_botCon.ski"
            elif cone_type=="Top":
                self.template   = "template_SEDInst_topCon.ski"
            else:
                print("Unrecognized cone type: ", cone_type)
                return
        else:
            self.output_instrument = "FullInstrument"
            if cone_type=='Full':
                self.template   = "template_FullInst.ski"
            else:
                print("For FullInstrument mode only a full cone_type can be used.")
                return

        #Start by reading the template.
        parser = ET.XMLParser(target=ET.TreeBuilder(insert_comments=True))
        self.tree = ET.parse(skg_folder+"/"+self.template, parser)
        self.root = self.tree.getroot()

        #Separate the main components. 
        self.tor, self.cone = self.root.findall("./MonteCarloSimulation/mediumSystem/MediumSystem/media/GeometricMedium")
        self.incs = self.root.findall("./MonteCarloSimulation/instrumentSystem/InstrumentSystem/instruments")[0]

        #Set the template to use for an inclination angle and remove it from the root.
        self.eta_temp = self.incs.find(self.output_instrument)
        self.incs.remove(self.eta_temp)

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
        if cone_type!="Full":
            self.base_fname += "_{}ConeOnly".format(cone_type)
        return
    
    def write_files_parameter_grid(self, tor_oa_min, tor_oa_max, tor_doa, cone_oa_min, cone_oa_max, cone_doa, etas = None, eta_max=90, eta_min=None, delta_eta=5, delta_eta_grazing=1, folder="scripts"):

        subprocess.call("mkdir {}".format(folder), shell=True)

        tor_oas  = np.arange(tor_oa_min , tor_oa_max +0.1*tor_doa , tor_doa )
        for tor_oa in tor_oas:

            #Set the torus opening angle. 
            self.tor.find("./geometry/TorusGeometry").set("openingAngle", "{} deg".format(90-tor_oa))

            #Remove all in the inclination blocks. These are left over from previous iterations. 
            for inc in self.incs.findall(self.output_instrument):
                self.incs.remove(inc)

            #Set the inclinations.
            if etas is None: 
                #eta_max = 90
                if eta_min is None or eta_min<tor_oa:
                    eta_min_use = tor_oa
                etas_use = np.arange(eta_min_use, eta_max+0.1*delta_eta, delta_eta)
                #etas[0]+=0.2*delta_eta
                if len(etas_use)>1:
                    etas_grazing = np.arange(etas_use[0], etas_use[1]-0.1*delta_eta_grazing, delta_eta_grazing)
                    etas_use = np.concatenate([etas_grazing, etas_use[1:]])
            else:
                etas_use = etas
            for eta in etas_use:
                eta_temp_use = copy.deepcopy(self.eta_temp)
                eta_temp_use.set("instrumentName","i{}".format(eta))
                eta_temp_use.set("inclination", "{} deg".format(eta))
                self.incs.append(eta_temp_use)

            cone_oa_max_use = np.min([cone_oa_max, tor_oa-5])
            cone_oas = np.arange(cone_oa_min, cone_oa_max_use+0.1*cone_doa, cone_doa)

            for cone_oa in cone_oas:

                #Set the cone geometry.
                if self.cone_type=="Full":
                    cone_prop_xml = "./geometry/ConicalShellGeometry"
                else:
                    cone_prop_xml = "./geometry/BoxClipGeometryDecorator/geometry/ConicalShellGeometry"
                self.cone.find(cone_prop_xml).set("maxAngle", "{} deg".format(90-cone_oa))
                self.cone.find(cone_prop_xml).set("minAngle", "{} deg".format(90-cone_oa - self.cone_width))

                fname = folder+"/"+self.base_fname + "_tor_oa{}_con_oa{}-tauV{}.ski".format(tor_oa, cone_oa, self.cone_tau)

                self.tree.write(fname, xml_declaration=True, encoding='UTF-8')

        return