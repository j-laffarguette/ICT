import os, sys
import random
import warnings

warnings.simplefilter(action='ignore', category=FutureWarning)
import pandas as pd

try:
    pid_file_path = os.path.join(os.environ.get('userprofile'), 'AppData', 'Local', 'Temp', 'raystation.pid')
    with open(pid_file_path) as f:
        os.environ['RAYSTATION_PID'] = f.read()
    script_client_path = r'C:\Program Files\RaySearch Laboratories\RayStation 11B\ScriptClient'
    sys.path.append(script_client_path)
except:
    print("")

from connect import *
import numpy as np
from easygui import multenterbox, ynbox
from datetime import datetime


# -------------------------------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------------------
# FONCTIONS
# -------------------------------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------------------

def check_roi(case, roi_to_check):
    # this method checks if a toi exists
    roi_check = False
    rois = case.PatientModel.RegionsOfInterest
    for roi in rois:
        if roi.Name == roi_to_check:
            roi_check = True
    return roi_check


def has_contour(case, examination, roi_to_check):
    """ Check if a structure is empty or not"""
    return case.PatientModel.StructureSets[examination].RoiGeometries[roi_to_check].HasContours()


def set_obj_function(plan, FunctionType, RoiName, DoseLevel, Weight, IsRobust, IsConstraint, HighDoseLevel):
    """
    Cette fonction permet de créer un objectif de dose et de remplir tous les paramètres
    :param plan: plan = get_current_plan()
    :param FunctionType: "MaxEud" ou autres
    :param RoiName: "nom de la roi"
    :param DoseLevel: entre 0 et 1 -> sera multiplié par la dose totale en cGy
    :param Weight: poids
    :param IsRobust: True ou False
    :param IsConstraint : True ou False
    :param HighDoseLevel: Entre 0 et 1 -> Pour le dose fall off
    :return:
    """

    # Création de la fonction
    plan.PlanOptimizations[0].AddOptimizationFunction(FunctionType=FunctionType, RoiName=RoiName,
                                                      IsConstraint=IsConstraint,
                                                      RestrictAllBeamsIndividually=False,
                                                      RestrictToBeam=None, IsRobust=IsRobust,
                                                      RestrictToBeamSet=None, UseRbeDose=False)

    # Remplissage des valeurs
    # On regarde le nombre de fonctions déjà présentes. En considérant que la fonction
    # venant juste d'être rentrée est la dernière. Attention, on ne cherche pas au même endroit pour les objectifs
    # et les contraintes

    if IsConstraint:
        n_functions = len(plan.PlanOptimizations[0].Constraints)
        dose = plan.PlanOptimizations[0].Constraints[n_functions - 1]
    else:
        n_functions = len(plan.PlanOptimizations[0].Objective.ConstituentFunctions)
        dose = plan.PlanOptimizations[0].Objective.ConstituentFunctions[n_functions - 1]

    if RoiName == dose.OfDoseGridRoi.OfRoiGeometry.OfRoi.Name:  # Vérification que tout va bien
        print(RoiName)

        if FunctionType == 'DoseFallOff':
            dose.DoseFunctionParameters.HighDoseLevel = HighDoseLevel
            dose.DoseFunctionParameters.LowDoseDistance = 0.5
            dose.DoseFunctionParameters.AdaptToTargetDoseLevels = True

        else:
            # On attribue la dose correspondante
            dose.DoseFunctionParameters.DoseLevel = DoseLevel
        dose.DoseFunctionParameters.Weight = Weight  # On met un weight même aux contraintes

    else:
        raise NameError("Le nom de la function dans le CSV ne correspond pas à la fonction créée. C'est la tuile...")
    print('-> Objectif créé !')


def round_to_nearest_half_int(num):
    # ✅ Round number to nearest 0.5
    return round(num * 2) / 2


# -------------------------------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------------------
# CLASSE
# -------------------------------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------------------

class Patient:
    def __init__(self):

        # Administratif Raystation

        self.case = get_current("Case")
        self.patient = get_current("Patient")
        self.db = get_current("PatientDB")
        self.external_name = None
        self.abdomen = []
        self.jonction = []
        self.pedia = False  # True si dossier pédia : trouvé par méthode self.pediatrique()
        self.poi_abdo = None
        self.poi_DSP = None
        self.lim_inf = None
        self.lim_sup = None
        self.upper_pallet = None
        self.direction = 'HFS'
        self.directory = r"G:\Commun\PHYSICIENS\JL\Projets\06 - Radixact\ICT\csv_files"

        # Récupération des informations des CT (nom + scan position HFS ou FFS)
        self.examinations = self.get_ct_list()  # dictionnaire {"exam_name":"FFS/HFS", ...}

        # Par défaut, à la première création de l'objet patient, le scanner Head First est pris en primary
        self.set_primary(self.examinations['HFS'])
        self.pediatrique()

        self.verifications()

    def verifications(self):
        # ----------------------------------------------------------
        # VERIF de la table. On ne vérifie que Upper Pallet
        # Verification de la présence de la roi upper pallet
        if not check_roi(self.case, 'Upper pallet'):
            raise NameError(f'Structures de table absentes. Veuillez les créer !')

        # --------------------------------------------------------------------------------
        # Verif du point "jonction" ou du point abdo sur les scanners non pédiatriques
        if not self.pedia:
            if not self.jonction:
                raise NameError(f'Point "jonction" manquant ou mal nommé. Veuillez le créer ou le modifier !')

            if not self.abdomen:
                raise NameError(f'Point "abdomen" manquant ou mal nommé. Veuillez le créer ou le modifier !')

        for direction in self.examinations:
            exam_name = self.examinations[direction]
            exam = self.case.Examinations[exam_name]

            # Verification de la présence de la table sur les deux scanners
            if not has_contour(self.case, exam_name, 'Upper pallet'):
                raise NameError(f'Structures de table absentes sur le scanner {exam_name}. Veuillez les créer !')

            # -------------------------------------------------------------------------------
            # verification de la presence des fichiers

            if direction == 'HFS':
                if not self.pediatrique:
                    if total_dose > 800:
                        filename = "adultesHFS.csv"
                    else:
                        filename = "2GyHFS.csv"
                else:
                    filename = "pediatrique.csv"
            else:
                filename = "adultesFFS.csv"

            if not os.path.isfile(os.path.join(self.directory, filename)):
                raise NameError(f'Fichier {filename} absent !')

    def pediatrique(self):
        # S'il n'y a qu'un scanner HFS et que le contour externe mesure moins de 130cm, on peut faire un seul plan
        self.pedia = False
        limite_de_taille = 130
        if 'FFS' not in self.examinations:
            print("\n---------------\n Ce dossier ne contient pas de scanner 'Feet First' ... ")
            self.lim_inf, self.lim_sup = [lim.z for lim in
                                          self.case.PatientModel.StructureSets[self.exam_name].RoiGeometries[
                                              self.external_name].GetBoundingBox()]
            if abs(self.lim_sup - self.lim_inf) <= limite_de_taille:
                self.pedia = ynbox("L'algo détecte que le patient mesure moins de 130 cm. Est-ce bien le cas ?")
                print("-------> Il s'agit d'un dossier pédiatrique")
            else:
                print("le patient semble trop grand (taille supérieure à 130cm). Nécessite une acquisition FFS")

    def set_primary(self, exam_name, direction=None):
        """Méthode très importante. Permet de définir un exam en 'primary et d'en récupérer les propriétés.
         Input : nom de l'examen à mettre en primary"""

        if direction is None:
            direction = self.direction
        else:
            self.direction = direction

        self.case.Examinations[exam_name].SetPrimary()
        self.examination = get_current("Examination")
        self.exam_name = self.examination.Name

        # -------------------------------------------------------------------------------------
        # Récupération de la position de tous les POI sur le scanner HFS
        # todo: à simplifier avec plusieurs typo

        if direction == 'HFS':
            round_it = True
        else:
            round_it = False

        self.jonction = self.get_point_coords('jonction', round_it)

        if not self.jonction:
            print('pas de point jonction!!!!!!!!!')

        self.poi_abdo = 'abdomen'
        self.abdomen = self.get_point_coords(self.poi_abdo)
        if not self.abdomen:
            self.poi_abdo = 'abdo'
            self.abdomen = self.get_point_coords(self.poi_abdo)

        # --------------------------------------------------------------------------------------
        # Création du point zero (abdo) situé à 0,0,hauteur table (valeur exprimée en mm convertie en cm)
        zero_scan = \
            self.case.Examinations[self.exam_name].GetStoredDicomTagValueForVerification(Group=0x0018, Element=0x1130)[
                'Table Height']
        self.zero_scan = -float(zero_scan) / 10

        # --------------------------------------------------------------------------------------
        # Récupération de la coordonnée la plus haute de la table (upper pallet) dans le but de mettre le laser vert
        # à 21 cm de ce point
        pallet = 'Upper pallet'

        if check_roi(self.case, pallet):
            if has_contour(self.case, self.exam_name, pallet):
                self.upper_pallet = self.case.PatientModel.StructureSets[self.exam_name].RoiGeometries[
                    pallet].GetBoundingBox()[0].y
            else:
                raise NameError('La structure de table est vide! Veuillez la créer!"')
        else:
            raise NameError("La structure de table n'existe pas. Veuillez la créer!")

        # --------------------------------------------------------------------------------------
        # Récupération du nom de l'externe
        for exam in self.case.PatientModel.StructureSets[self.exam_name].RoiGeometries:
            if exam.OfRoi.Type == 'External':
                self.external_name = exam.OfRoi.Name

        # Si pas d'external, faire un message d'erreur
        if self.external_name is None:
            raise NameError("Contour externe absent. Veuillez le créer ! \n"
                            "-> Remarques importantes :\n "
                            "- Vérifier les cavités d'aires (absence de trous)\n"
                            "- Vérifier que les doigts n'ont pas été supprimés lors de la simplification\n"
                            "- Vérifier l'absence d'externe dans les zones denses de l'omniboard")

            # external_name = "External_Auto"
            # retval_0 = self.case.PatientModel.CreateRoi(Name=external_name, Color="Green", Type="External",
            #                                             TissueName="",
            #                                             RbeCellTypeName=None, RoiMaterial=None)
            # retval_0.CreateExternalGeometry(Examination=self.examination, ThresholdLevel=-250)
            # self.external_name = external_name

        # si le contour existe, mais qu'il est vide, il faut le créer
        if not has_contour(self.case, self.exam_name, self.external_name):
            self.case.PatientModel.RegionsOfInterest[self.external_name].CreateExternalGeometry(
                Examination=self.examination, ThresholdLevel=-250)

        self.lim_inf, self.lim_sup = [lim.z for lim in
                                      self.case.PatientModel.StructureSets[self.exam_name].RoiGeometries[
                                          self.external_name].GetBoundingBox()]
        print(f"External name : {self.external_name}")

    def create_dsp(self):
        """
        Création du point DSP si pas pédiatrique et si n'existe pas déjà
        :return:
        """
        if self.examinations['FFS'] == self.exam_name:
            self.poi_DSP = [self.jonction[0], self.jonction[1] + 3, self.jonction[2]]
            self.create_poi('dsp', self.poi_DSP, color='Black')

    def get_zero_scan(self, scan_direction):
        """Méthode permettant de récupérer le zéro du scanner. À noter que sur le scanner siemens somatom, le zéro n'est
        pas à zéro en antépost, il est situé à (0,0,hauteur table)
        Input : scan_direction : 'HFS' ou 'FFS'
        Output: position du zéro en antépost (en cm)
        """
        zero_scan = \
            self.case.Examinations[self.examinations[scan_direction]].GetStoredDicomTagValueForVerification(
                Group=0x0018,
                Element=0x1130)[
                'Table Height']

        return -float(zero_scan) / 10

    def create_cylinder_ptv(self, roi_name, y_cranial, y_caudal):
        """Méthode permettant de déterminer le centre du cylindre à partir des coordonnées du point le plus haut et
        du point le plus bas en longitudinal. Le script est codé de sorte que l'on dise : je veux que le PTV démarre
        de ce point (ex: bille genou) et qu'il aille jusqu'à ce point (ex: bille jonction)
        IMPORTANT : cette méthode appelle la méthode cylindre permettant de créer les PTVS.

        Input : * roi_name = Nom du PTV attendu. Attention, la Roi doit déjà exister dans le roi set
                * y_cranial : point le plus supérieur du PTV (en cm)
                * y_caudal : point le plus inférieur (en cm)
        """

        x0, z0 = 0, self.zero_scan
        long = y_cranial - y_caudal  # longueur du volume
        center = y_cranial - long / 2  # centre du volume
        self.cylinder(roi_name, (x0, z0, center), longueur=abs(long))

    def cylinder(self, roi_name, coords, longueur=2, retraction=True):
        """Création des PTV en utilisant la méthode des cylindres. Un grand cylindre est créé (40cm de rayon dans le
        plan axial et une longueur par défaut de 2cm dans la direction longitudinale). Le cylindre ainsi créé est ensuite
        rétracté à la peau de 3 mm afin de créer le PTV final.\n
        - Inputs : * roi_name = Nom du PTV attendu. Attention, la Roi doit déjà exister dans le roi set.\
            * coords = coordonnées du centre du cylindre
            * longueur = 2 : longueur du PTV dans la direction longitudinale (en cm)
            * rétraction = True : Rétraction ou pas
        """
        x, y, z = coords
        self.case.PatientModel.RegionsOfInterest[roi_name].CreateCylinderGeometry(Radius=40,
                                                                                  Axis={'x': 0, 'y': 0, 'z': 1},
                                                                                  Length=longueur,
                                                                                  Examination=self.examination,
                                                                                  Center={'x': x,
                                                                                          'y': y,
                                                                                          'z': z},
                                                                                  Representation="TriangleMesh",
                                                                                  VoxelSize=None)

        if retraction:
            self.algebra(out_roi=roi_name, in_roiA=[roi_name], in_roiB=[obj_patient.external_name], margeB=-0.3,
                         ResultOperation="Intersection", derive=False)

    def get_point_coords(self, point_name, round_it=True):
        """ Méthode permettant de récupérer les coordonnées d'un point si celui-ci existe. \n
        - Input: nom du point \n
        - Outupt: coordonnées du point (x,y,z) ou None"""
        try:
            point = self.case.PatientModel.StructureSets[self.exam_name].PoiGeometries[point_name]
            coords = point.Point
        except:
            coords = None
            print(f'Point {point_name} absent de la collection ...')

        # Un point peut exister mais ne pas avoir de coordonnées dans l'exam!
        if coords is not None:
            print(f'---> Point {point_name} trouvé! ...')

            # Pour arrondir les coordonnées et modifier le point!
            if round_it:
                x, y, z = coords.x, coords.y, coords.z
                x, y, z = x, y, round_to_nearest_half_int(z)
                self.case.PatientModel.StructureSets[self.exam_name].PoiGeometries[point_name].Point = {'x': x, 'y': y,
                                                                                                        'z': z}

            return coords.x, coords.y, coords.z
        else:
            return None

    def algebra(self, out_roi, in_roiA=[], in_roiB=[], margeA=0, margeB=0, color='Blue', ResultOperation='None',
                Type="Organ", derive=True):
        """Méthode permettant de réaliser des algebra and marging de manière quasi généralisée """

        if not check_roi(self.case, out_roi):
            self.case.PatientModel.CreateRoi(Name=out_roi, Color=color, Type=Type, TissueName=None,
                                             RbeCellTypeName=None, RoiMaterial=None)

        # Ici, l'utilisateur entre une valeur réelle. Si elle est négative, RS demande à ce qu'elle soit positive avec
        # comme type : 'contract'
        if margeA < 0:
            typeA = 'Contract'
            margeA = abs(margeA)
        else:
            typeA = 'Expand'

        if margeB < 0:
            typeB = 'Contract'
            margeB = abs(margeB)
        else:
            typeB = 'Expand'

        # Si le volume de sortie dérive des volumes en entrée
        if derive:
            retval_0 = self.case.PatientModel.RegionsOfInterest[out_roi].SetAlgebraExpression(
                ExpressionA={'Operation': "Union", 'SourceRoiNames': in_roiA,
                             'MarginSettings': {'Type': typeA, 'Superior': margeA, 'Inferior': margeA,
                                                'Anterior': margeA,
                                                'Posterior': margeA, 'Right': margeA, 'Left': margeA}},
                ExpressionB={'Operation': "Union", 'SourceRoiNames': in_roiB,
                             'MarginSettings': {'Type': typeB, 'Superior': margeB, 'Inferior': margeB,
                                                'Anterior': margeB,
                                                'Posterior': margeB, 'Right': margeB, 'Left': margeB}},
                ResultOperation=ResultOperation,
                ResultMarginSettings={'Type': "Expand", 'Superior': 0, 'Inferior': 0, 'Anterior': 0, 'Posterior': 0,
                                      'Right': 0, 'Left': 0})

            retval_0.UpdateDerivedGeometry(Examination=self.examination, Algorithm="Auto")

        # Si le volume de sortie ne dérive pas.
        else:
            self.case.PatientModel.RegionsOfInterest[out_roi].CreateAlgebraGeometry(
                Examination=self.examination, Algorithm="Auto",
                ExpressionA={'Operation': "Union", 'SourceRoiNames': in_roiA,
                             'MarginSettings': {'Type': typeA, 'Superior': margeA, 'Inferior': margeA,
                                                'Anterior': margeA,
                                                'Posterior': margeA, 'Right': margeA, 'Left': margeA}},
                ExpressionB={'Operation': "Union", 'SourceRoiNames': in_roiB,
                             'MarginSettings': {'Type': typeB, 'Superior': margeB, 'Inferior': margeB,
                                                'Anterior': margeB,
                                                'Posterior': margeB, 'Right': margeB, 'Left': margeB}},
                ResultOperation=ResultOperation,
                ResultMarginSettings={'Type': "Expand", 'Superior': 0, 'Inferior': 0, 'Anterior': 0, 'Posterior': 0,
                                      'Right': 0, 'Left': 0})

    def get_ct_list(self):
        """ Méthode permettant de récupérer les noms des scanners en fonction de leur position (HFS ou FFS).
         Output: self.examinations = dictionnaire {"nom_duCT_HFS":"HFS", "nom_duCT_FFS":"FFS"}"""
        examinations = {}
        for exam in self.case.Examinations:
            # On ne travaille que sur les images CT
            if exam.EquipmentInfo.Modality == 'CT':
                data = exam.GetAcquisitionDataFromDicom()
                description = data['SeriesModule']['SeriesDescription']
                print(exam.Name, description)
                # On ne prend que les images reconstruites en mdd 
                if "mdd" in description.lower() :
                    examinations[exam.PatientPosition] = exam.Name
        self.examinations = examinations

        if not examinations:
            raise NameError("Aucun scanner trouvé. Deux causes possibles:\n"
                            "1 - Le ou les scanners n'ont pas 'HFS' ou 'FFS' en orientation\n"
                            "2 - Aucun scanner 'mdd' n'a été trouvé (scanner Br38 au lieu de Sn40)")


        return self.examinations

    def create_roi(self, roi_name, color=None, roi_type='Ptv'):
        """Méthode utilisée pour créer des volumes de type PTV
        Si on ne donne pas de couleur, celle-ci est prise aléatoirement"""

        if color is None:
            color = ["#" + ''.join([random.choice('ABCDEF0123456789') for i in range(6)])][0]

        self.case.PatientModel.CreateRoi(Name=roi_name, Color=color, Type=roi_type)

    def create_poi(self, poi_name, coords, color="128, 128, 255"):
        """Méthode utilisée pour créer des points (exemple: laser rouge et laser vert). Si le point n'est pas présent
        dans la poi list, le point est créé, sinon ses coordonnées sont modifiées

        Input:  * poi_name : nom du point
                * coords : coordonnées du point
                * color : couleur données au point (soit en valeurs RGB soit 'Yellow')
                """
        x, y, z = coords
        poi_list = [poi.Name for poi in self.case.PatientModel.PointsOfInterest]
        if not poi_name in poi_list:
            retval_0 = self.case.PatientModel.CreatePoi(Examination=self.examination,
                                                        Point={'x': x, 'y': y, 'z': z},
                                                        Name=poi_name, Color=color,
                                                        VisualizationDiameter=1, Type="Undefined")
        else:
            self.case.PatientModel.StructureSets[self.exam_name].PoiGeometries[poi_name].Point = {
                'x': x, 'y': y, 'z': z}


# -------------------------------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------------------
# MAIN
# -------------------------------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------------------

if __name__ == '__main__':

    # Création de l'objet patient. Définit par défaut le scanner HFS en primary
    obj_patient = Patient()

    # Enfant ou adulte? True or False
    pediatrique = obj_patient.pedia

    print("\n-----------------------------\nCREATION DE LA DOSIMETRIE ICT\n----------------------------- \n")
    print("Les fichiers csv contenant les objectifs et contraintes dosimétriques seront lus dans le répertoire :\n"
          rf"->  {obj_patient.directory}", "\n")

    # définition de la prescription (2Gy ou 12 Gy?)
    verif = 0
    msg = "Entrer les informations de la prescription"
    title = "Prescription"
    fields = ("Nombre de fractions", "Dose totale (Gy)")

    easygui = True  # mettre False pour ne pas afficher la fenêtre en mode programmation !

    # Cette section permet d'afficher une fenêtre permettant d'entrer la prescription
    # todo: remplacer par la méthode automatique de requête Mosaiq.
    if easygui:
        while verif != 1:
            mes_choix = multenterbox(msg, title, fields)
            number_of_fractions = int(mes_choix[0])
            total_dose = float(mes_choix[1]) * 100
            fraction_dose = total_dose / number_of_fractions

            if fraction_dose != 200:
                msg = "Une erreur a été commise. \nVérifiez les données entrées et recommencez !"
            else:
                break
    else:
        number_of_fractions = 6
        total_dose = 1200
        fraction_dose = 200

    # ------------------------------------------------------------------------------------
    # On commencera d'abord sur le scanner Head First puis on travaille sur le Feet First. SI pédiatrique, on ne
    # travaille que sur le HFS. La liste to_do est utilisée plus loin dans le script.
    if pediatrique:
        to_do = ['HFS']
    else:
        to_do = ['HFS', 'FFS']

    # do_it est mis sur False pour sauter toute la partie Patient Modeling pour la programmation du script
    do_it = True

    if do_it:
        #########################################################################
        #########################################################################
        ##################### TRAVAIL DE PREPARATION ############################
        #########################################################################
        #########################################################################

        if not pediatrique:
            # Première partie, réalisation du recalage rigide
            # Suppression pure et simple de tous les FOR registration déjà existants

            for reg in obj_patient.case.Registrations:
                FFOR = reg.FromFrameOfReference
                RFOR = reg.ToFrameOfReference
                obj_patient.case.RemoveFrameOfReferenceRegistration(
                    FloatingFrameOfReference=FFOR,
                    ReferenceFrameOfReference=RFOR)
                print('-> Un recalage a été supprimé')

            # Création du recalage entre le scanner FFS et le scanner HFS
            # 1. création du for reg
            print(f'Recalage rigide entre {obj_patient.examinations["HFS"]} et {obj_patient.examinations["FFS"]}')
            obj_patient.case.CreateNamedIdentityFrameOfReferenceRegistration(
                FromExaminationName=obj_patient.examinations['HFS'], ToExaminationName=obj_patient.examinations['FFS'],
                RegistrationName="HFS to FFS", Description=None)

            # recalage automatique
            obj_patient.case.ComputeGrayLevelBasedRigidRegistration(
                FloatingExaminationName=obj_patient.examinations['HFS'],
                ReferenceExaminationName=obj_patient.examinations[
                    'FFS'],
                UseOnlyTranslations=True, HighWeightOnBones=True,
                InitializeImages=True, FocusRoisNames=[],
                RegistrationName=None)
            print("-> Ok!")

            # Deuxième partie : réalisation du recalage élastique "rigide". On utilise la méthode discard intensity
            # avec la vessie comme volume d'intérêt. Pour cela, la vessie doit être copiée d'un scanner à l'autre.
            print(f"Copie de la structure Vessie d'un scanner à l'autre ...")
            obj_patient.case.PatientModel.CopyRoiGeometries(SourceExamination=obj_patient.examination,
                                                            TargetExaminationNames=[obj_patient.examinations['FFS']],
                                                            RoiNames=["Vessie"], ImageRegistrationNames=[],
                                                            TargetExaminationNamesToSkipAddedReg=[
                                                                obj_patient.examinations['FFS']])

            # Création du recalage élastique rigide
            print(f"Réalisation d'un recalage élastique (rigide, basé sur la vessie comme structure de contrôle) ...")
            reg_name = "elastique rigide"
            obj_patient.case.PatientModel.CreateHybridDeformableRegistrationGroup(RegistrationGroupName=reg_name,
                                                                                  ReferenceExaminationName=
                                                                                  obj_patient.examinations['HFS'],
                                                                                  TargetExaminationNames=[
                                                                                      obj_patient.examinations["FFS"]],
                                                                                  ControllingRoiNames=["Vessie"],
                                                                                  ControllingPoiNames=[],
                                                                                  FocusRoiNames=[],
                                                                                  AlgorithmSettings={
                                                                                      'NumberOfResolutionLevels': 1,
                                                                                      'InitialResolution': {'x': 0.5,
                                                                                                            'y': 0.5,
                                                                                                            'z': 0.5},
                                                                                      'FinalResolution': {'x': 0.25,
                                                                                                          'y': 0.25,
                                                                                                          'z': 0.5},
                                                                                      'InitialGaussianSmoothingSigma': 2,
                                                                                      'FinalGaussianSmoothingSigma': 0.333333333333333,
                                                                                      'InitialGridRegularizationWeight': 400,
                                                                                      'FinalGridRegularizationWeight': 400,
                                                                                      'ControllingRoiWeight': 0.5,
                                                                                      'ControllingPoiWeight': 0.1,
                                                                                      'MaxNumberOfIterationsPerResolutionLevel': 1000,
                                                                                      'ImageSimilarityMeasure': "None",
                                                                                      'DeformationStrategy': "Default",
                                                                                      'ConvergenceTolerance': 1E-05})
            print("-> Ok!")

            # mapping des POI d'un scanner vers l'autre
            print("Mapping des points jonction et abdo...")
            pois = ["jonction", obj_patient.poi_abdo]

            obj_patient.case.MapPoiGeometriesDeformably(PoiGeometryNames=pois,
                                                        CreateNewPois=False,
                                                        StructureRegistrationGroupNames=[reg_name],
                                                        ReferenceExaminationNames=[obj_patient.examinations['HFS']],
                                                        TargetExaminationNames=[obj_patient.examinations['FFS']],
                                                        ReverseMapping=False, AbortWhenBadDisplacementField=False)

        # Travail sur les volumes
        print("\nTravail sur les ROI...")
        # Création des ROI
        if not pediatrique:  # Patients adultes

            if total_dose > 800:  # Prescription 12 Gy
                ROI_LIST = ['PTV FFS', 'PTV_4', 'PTV_3', 'PTV_2', 'PTV_1', "PTV HFS", 'PTV poumons', 'PTV reins',
                            'opt PTV HFS']
                colors = ['#FF80FF', "#FF8080", "#FFFF80", "#00FF80", "#00FFFF", "#0080C0", "#66FFFF", "#666633",
                          '#696100']

            else:  # Prescriptions 2 ou 8 Gy. On n'a plus besoin d'OAR car la dose est inférieure ou égale
                # aux contraintes
                ROI_LIST = ['PTV FFS', 'PTV_4', 'PTV_3', 'PTV_2', 'PTV_1', "PTV HFS"]
                colors = ['#FF80FF', "#FF8080", "#FFFF80", "#00FF80", "#00FFFF", "#0080C0"]

        else:  # Enfants
            ROI_LIST = ['PTV HFS', 'PTV poumons', 'PTV robustesse', 'PTV reins', 'opt PTV HFS']
            colors = ["#0080C0", "#66FFFF", None, "#666633", '#696100']

        # Vérification de l'existance des differentes roi dans le case
        resultat = [check_roi(obj_patient.case, roi) for roi in ROI_LIST]
        # Création des ROI dans le roi set si ROI non existantes
        [obj_patient.create_roi(ROI_LIST[index], color=colors[index]) for index, roi in enumerate(resultat) if not
        roi]

        # ------------------------------------------------------------------------------------
        # On commence d'abord sur le scanner Head First puis on travaille sur le Feet First. Si pédiatrique, on ne
        # travaille que sur le HFS

        for direction in to_do:
            # ct_name pour simplification dans la boucle for
            ct_name = obj_patient.examinations[direction]
            ct = obj_patient.case.Examinations[ct_name]

            # ct étudié mis en primary (très important, car re-définit self.exam_name et self.examination dans la classe
            obj_patient.set_primary(ct_name, direction)

            # Retrait des trous dans externe
            obj_patient.case.PatientModel.StructureSets[ct_name].SimplifyContours(RoiNames=[obj_patient.external_name],
                                                                                  RemoveHoles3D=True,
                                                                                  RemoveSmallContours=False,
                                                                                  AreaThreshold=None,
                                                                                  ReduceMaxNumberOfPointsInContours=False,
                                                                                  MaxNumberOfPoints=None,
                                                                                  CreateCopyOfRoi=False,
                                                                                  ResolveOverlappingContours=True)

            # créer poumon G+D. On réalise +0.5 / -0.5 pour supprimer les défauts du volume
            obj_patient.algebra(out_roi="Poumons", in_roiA=["Poumon D", "Poumon G"], margeA=0.5, color="Aqua",
                                derive=False)
            obj_patient.algebra(out_roi="Poumons", in_roiA=["Poumons"], margeA=-0.5, color="Aqua",
                                derive=False)
            # Remove holes
            obj_patient.case.PatientModel.StructureSets[obj_patient.exam_name].SimplifyContours(RoiNames=["Poumons"],
                                                                                                RemoveHoles3D=True,
                                                                                                RemoveSmallContours=False,
                                                                                                AreaThreshold=None,
                                                                                                ReduceMaxNumberOfPointsInContours=False,
                                                                                                MaxNumberOfPoints=None,
                                                                                                CreateCopyOfRoi=False,
                                                                                                ResolveOverlappingContours=False)

            # créer somme des reins
            obj_patient.algebra(out_roi="Reins", in_roiA=["Rein D", "Rein G"], color="Yellow")

            if not pediatrique:
                # ------------------------------------------------------------------------------------
                # Création des cylindres
                # récupération du point "jonction", positionné sur la bille au niveau de la jonction
                _, _, y0 = obj_patient.jonction  # y de ref est au niveau de la jonction
                x0, z0 = 0, obj_patient.zero_scan

                # Réalisation des PTV 2 à 5, entourant la bille jonction
                for index, roi in enumerate(['PTV_1', 'PTV_2', 'PTV_3', 'PTV_4']):
                    # y0 est au niveau de la premiere coupe de PTV_2
                    # la premiere coupe de PTV_3 est donc située à -2; celle de D1 à +2 etc
                    # soit -2 ; 0 ; 2 ; 4 -> 4-2*i
                    yhaut = y0 + 4 - 2 * index
                    y_bas = yhaut - 2
                    obj_patient.create_cylinder_ptv(roi, yhaut, y_bas)

                # Liste qui contiendra pour chaque volume, son point de départ et son point d'arrivée
                from_to = []

                # lim sup et inf correspondent au points les plus hauts et bas du contour externe
                # Réalisation du "PTV HFS", au-dessus du dernier PTV de la jonction, jusqu'au sommet du crane
                haut_PTV1 = y0 + 4
                from_to.append(["PTV HFS", obj_patient.lim_sup, haut_PTV1])

                # Réalisation du PTV FFS : en dessous du PTV 4
                from_to.append(['PTV FFS', y0 - 4, obj_patient.lim_inf])

                # Une fois que toutes les bornes sont déterminées, on crée les PTVS
                for roi, f, t in from_to:
                    print(f'Création du PTV suivant : {roi} entre y = {str(f)} et y = {str(t)} ')
                    obj_patient.create_cylinder_ptv(roi, f, t)
                    print('-> OK!')

                # Dernière étape: soustraction des volumes deux à deux pour éviter le chevauchement des PTVS
                ROI_LIST = ['PTV FFS', 'PTV_4', 'PTV_3', 'PTV_2', 'PTV_1', "PTV HFS"]
                for i in range(len(ROI_LIST) - 1):
                    obj_patient.algebra(out_roi=ROI_LIST[i], in_roiA=[ROI_LIST[i]], in_roiB=[ROI_LIST[i + 1]],
                                        ResultOperation="Subtraction", derive=False)

                # Création d'un volume d'optimisation en sortie de jonction
                # HFS -> PTV_D6 devient "opt_jonction"
                # FFS -> PTV D1 devient "opt_jonction"

                if direction == "HFS":
                    source = "PTV FFS"
                    OAR_name = "opt_jonction"
                else:
                    source = "PTV HFS"
                    OAR_name = "opt_jonction"

                # Création du volume
                obj_patient.algebra(out_roi=OAR_name, in_roiA=[source], color='Magenta', Type="Organ")
                obj_patient.case.PatientModel.RegionsOfInterest[OAR_name].DeleteExpression()

            # Si pédiatrique, seul le PTV HFS compte et on prend juste external - 3mm
            elif pediatrique:
                obj_patient.algebra(out_roi="PTV HFS", in_roiA=[obj_patient.external_name], margeA=-0.3, derive=True)

                # On réalise aussi un PTV jambe pour la robustesse. Le PTV va du bout des pieds jusqu'à 15cm au dessus
                obj_patient.create_cylinder_ptv('PTV robustesse', obj_patient.lim_inf + 15, obj_patient.lim_inf)

            # simplification des volumes pour éviter overlaps
            roi_list = [roi.Name for roi in obj_patient.case.PatientModel.RegionsOfInterest if roi.Type == "Support"]

            obj_patient.case.PatientModel.StructureSets[obj_patient.exam_name].SimplifyContours(
                RoiNames=roi_list,
                RemoveHoles3D=False, RemoveSmallContours=False, AreaThreshold=None,
                ReduceMaxNumberOfPointsInContours=False, MaxNumberOfPoints=None, CreateCopyOfRoi=False,
                ResolveOverlappingContours=True)

            if direction == 'HFS' and total_dose > 800:

                # Création PTV poumons = poumons - 1cm # seulement pour le scanner HFS!
                obj_patient.algebra(out_roi='PTV poumons', in_roiA=["Poumons"], margeA=-1, Type="Ptv", color="Yellow")

                # PTV paroi thoracique = paroi thoracique + expansion du poumon de 1.5 cm + coeur. Ce volume est
                # utilisé pour que l'optimisation sur les poumons ne diminue pas trop la couverture aux alentours

                # 1. Expansion des poumons
                obj_patient.algebra(out_roi='opt pourtour poumons', in_roiA=["Poumons"], margeA=1.5,
                                    in_roiB=["Poumons"],
                                    margeB=0, ResultOperation="Subtraction", derive=False)
                # 2. Union du PTV, de l'expansion et du coeur puis rétraction à la peau
                obj_patient.algebra(out_roi='PTV Paroi thoracique', in_roiA=["opt pourtour poumons", "Coeur"],
                                    in_roiB=[obj_patient.external_name], margeB=-0.3, ResultOperation="Intersection",
                                    Type="Ptv", derive=False)

                # Création PTV reins = reins - 1cm # seulement pour le scanner HFS!
                obj_patient.algebra(out_roi='PTV reins', in_roiA=["Reins"], margeA=-1, Type="Ptv", color="Yellow")

                # PTV HFS - PTV Poumons - PTV reins -  = PTV HFS
                obj_patient.algebra(out_roi='PTV HFS', in_roiA=["PTV HFS"], in_roiB=["PTV poumons", 'PTV reins'],
                                    ResultOperation="Subtraction", derive=False)

                # opt PTV HFS = PTV HFS - Poumons - Reins
                obj_patient.algebra(out_roi='opt PTV HFS', in_roiA=["PTV HFS"], in_roiB=["Poumons", 'Reins'],
                                    ResultOperation="Subtraction", derive=False)
            else:
                try:
                    # Supprimer la dérivation
                    obj_patient.case.PatientModel.RegionsOfInterest['PTV poumons'].DeleteExpression()
                    obj_patient.case.PatientModel.RegionsOfInterest['PTV reins'].DeleteExpression()
                    obj_patient.case.PatientModel.RegionsOfInterest['Reins'].DeleteExpression()
                except:
                    print("Impossible de supprimer la dérivation")

    #########################################################################
    #########################################################################
    ################ TRAVAIL SUR LE PLAN DE TRAITEMENT ######################
    #########################################################################
    #########################################################################

    # Création des deux plans
    date = datetime.today().strftime('%Y%m%d')
    date = date[2:]  # Pour enlever les deux premiers chiffres de l'année (2022 -> 22)
    PLAN_NAMES = [f"{date}_HFS", f"{date}_FFS"]
    BS_NAMES = ["r1", "r1"]
    PATIENT_POSITIONS = ["HeadFirstSupine", "FeetFirstSupine"]
    machine_name = "Radixact1"

    # pour rappel: to_do == HFS pour pédia et [HFS,FFS] pour adultes > 130 cm
    for ind, direction in enumerate(to_do):

        # CT étudié en primary
        obj_patient.set_primary(obj_patient.examinations[direction])

        # Attribution de l'IVDT qui va bien
        if not obj_patient.examination.EquipmentInfo.ImagingSystemReference:
            obj_patient.examination.EquipmentInfo.SetImagingSystemReference(ImagingSystemName="SIEMENS_RT_mDD")
            print('~~~~ Saving case ~~~~')
            obj_patient.patient.Save()
        elif obj_patient.examination.EquipmentInfo.ImagingSystemReference.get('ImagingSystemName') == 'SIEMENS_RT_mDD':
            obj_patient.examination.EquipmentInfo.SetImagingSystemReference(ImagingSystemName="SIEMENS_RT_mDD")
            print('~~~~ Saving case ~~~~')
            obj_patient.patient.Save()

        # Attribution des noms de plans et de beamsets
        plan_name = PLAN_NAMES[ind]
        bs_name = BS_NAMES[ind]

        # vérification de la non-existence du plan avant création
        list_of_plans = [plan.Name for plan in obj_patient.case.TreatmentPlans]
        if plan_name not in list_of_plans:
            # création du plan
            retval_0 = obj_patient.case.AddNewPlan(PlanName=plan_name, PlannedBy="", Comment="",
                                                   ExaminationName=obj_patient.exam_name,
                                                   IsMedicalOncologyPlan=False, AllowDuplicateNames=False)

        # vérification de la non-existence du bs avant création
        list_of_bs = [bs.DicomPlanLabel for bs in obj_patient.case.TreatmentPlans[plan_name].BeamSets]
        if bs_name not in list_of_bs:
            # création du bs
            if direction == "HFS":
                obj_patient.case.TreatmentPlans[plan_name].AddNewBeamSet(Name=bs_name,
                                                                         ExaminationName=obj_patient.exam_name,
                                                                         MachineName=machine_name,
                                                                         Modality="Photons",
                                                                         TreatmentTechnique="TomoHelical",
                                                                         PatientPosition=PATIENT_POSITIONS[ind],
                                                                         NumberOfFractions=number_of_fractions,
                                                                         CreateSetupBeams=False,
                                                                         UseLocalizationPointAsSetupIsocenter=False,
                                                                         UseUserSelectedIsocenterSetupIsocenter=False,
                                                                         Comment="", RbeModelName=None,
                                                                         EnableDynamicTrackingForVero=False,
                                                                         NewDoseSpecificationPointNames=[],
                                                                         NewDoseSpecificationPoints=[],
                                                                         MotionSynchronizationTechniqueSettings={
                                                                             'DisplayName': None,
                                                                             'MotionSynchronizationSettings': None,
                                                                             'RespiratoryIntervalTime': None,
                                                                             'RespiratoryPhaseGatingDutyCycleTimePercentage': None,
                                                                             'MotionSynchronizationTechniqueType': "Undefined"},
                                                                         Custom=None, ToleranceTableLabel=None)
            else:
                obj_patient.case.TreatmentPlans[plan_name].AddNewBeamSet(Name=bs_name,
                                                                         ExaminationName=obj_patient.exam_name,
                                                                         MachineName=machine_name,
                                                                         Modality="Photons",
                                                                         TreatmentTechnique="TomoDirect",
                                                                         PatientPosition="FeetFirstSupine",
                                                                         NumberOfFractions=number_of_fractions,
                                                                         CreateSetupBeams=False,
                                                                         UseLocalizationPointAsSetupIsocenter=False,
                                                                         UseUserSelectedIsocenterSetupIsocenter=False,
                                                                         Comment="",
                                                                         RbeModelName=None,
                                                                         EnableDynamicTrackingForVero=False,
                                                                         NewDoseSpecificationPointNames=[],
                                                                         NewDoseSpecificationPoints=[],
                                                                         MotionSynchronizationTechniqueSettings={
                                                                             'DisplayName': None,
                                                                             'MotionSynchronizationSettings': None,
                                                                             'RespiratoryIntervalTime': None,
                                                                             'RespiratoryPhaseGatingDutyCycleTimePercentage': None,
                                                                             'MotionSynchronizationTechniqueType': "Undefined"},
                                                                         Custom=None, ToleranceTableLabel=None)

        # Sauvegarde du plan
        print('~~~~ Saving case ~~~~')
        obj_patient.patient.Save()

        # Le plan créé est mis en primary pour simplifier
        obj_patient.case.TreatmentPlans[plan_name].SetCurrent()
        plan = get_current('Plan')

        # Le bs créé est mis en primary pour simplifier
        plan.BeamSets[bs_name].SetCurrent()
        beam_set = get_current("BeamSet")

        # Optimization settings
        plan.PlanOptimizations[0].OptimizationParameters.Algorithm.OptimalityTolerance = 1e-7
        plan.PlanOptimizations[0].OptimizationParameters.Algorithm.MaxNumberOfIterations = 30
        plan.PlanOptimizations[0].OptimizationParameters.DoseCalculation.ComputeFinalDose = True

        # if ind = 0 = HFS la prescription est attribuée au PTV HFS, sinon au PTV FFS pour le plan FFS
        if direction == 'HFS':
            prescription_roi = 'PTV HFS'
        else:
            prescription_roi = 'PTV FFS'

        # Si la prescription existe déjà, ne pas la recréer
        if not beam_set.Prescription.PrescriptionDoseReferences:
            beam_set.AddRoiPrescriptionDoseReference(RoiName=prescription_roi, DoseVolume=0,
                                                     PrescriptionType="MedianDose",
                                                     DoseValue=total_dose, RelativePrescriptionLevel=1)
        elif not beam_set.Prescription.PrescriptionDoseReferences[0].OnStructure.Name == prescription_roi:
            beam_set.AddRoiPrescriptionDoseReference(RoiName=prescription_roi, DoseVolume=0,
                                                     PrescriptionType="MedianDose",
                                                     DoseValue=total_dose, RelativePrescriptionLevel=1)

        beam_set.SetAutoScaleToPrimaryPrescription(AutoScale=False)

        ########################################################
        # définition des points de ref

        # ------------------------------------------------------------------------------------
        # Positionnement du localization point HFS -> LASER ROUGE
        # Le loc point du traitement HFS est en Inf/Sup sur les billes thorax, à zéro en G/D
        # Il est situé en ant/post sur les billes cranes = zero scan = table height (voir __init__ de Patient())

        # modification du type de tous les points en Undefined au cazou
        for poi in obj_patient.case.PatientModel.PointsOfInterest:
            poi.Type = 'Undefined'

        laser_rouge_HFS = "laser rouge"

        # On positionne les lasers rouges et vert sur le point abdomen mais le rouge est décalé sur la bille
        # jonction en GD
        coords_laser_rouges_HFS = (obj_patient.jonction[0], obj_patient.zero_scan, obj_patient.abdomen[2])

        if direction == 'HFS':  # Plan HFS
            # Le laser vert est mis en antépost au centre vertical du volume "Poumons"
            AP_iso_HFS = obj_patient.case.PatientModel.StructureSets[obj_patient.examinations['HFS']].RoiGeometries[
                'Poumons'].GetCenterOfRoi().y
            coords_laser_vert = (0, AP_iso_HFS, obj_patient.abdomen[2])

        elif direction == 'FFS':  # Plan FFS
            # Pour le plan FFS, le laser vert est positionné en AP au centre du volume PTV ffs afin que le patient
            # soit bien centré en hauteur et n'ait pas les genoux qui dépassent du cadre. Mais attention, ne doit pas
            # dépasser 21.5 cm par rapport à la table

            AP_iso_FFS = obj_patient.case.PatientModel.StructureSets[obj_patient.examinations['FFS']].RoiGeometries[
                'PTV FFS'].GetCenterOfRoi().y

            if abs(AP_iso_FFS) > abs(obj_patient.upper_pallet) + 21:
                AP_iso_FFS = obj_patient.upper_pallet - 21

            coords_laser_vert = (0, AP_iso_FFS, obj_patient.abdomen[2])

        # création du laser rouge
        obj_patient.create_poi(laser_rouge_HFS, coords_laser_rouges_HFS, color='Red')
        # laser rouge devient localization point
        obj_patient.case.PatientModel.PointsOfInterest[laser_rouge_HFS].Type = "LocalizationPoint"
        # création du laser vert
        obj_patient.create_poi('laser vert', coords_laser_vert, color='Green')

        ########################################################
        # Création du faisceau

        plan = obj_patient.case.TreatmentPlans[plan_name]
        x, y, z = coords_laser_vert

        if direction == 'HFS':
            try:
                retval_0 = beam_set.CreatePhotonBeam(BeamQualityId="6", CyberKnifeCollimationType="Undefined",
                                                     CyberKnifeNodeSetName=None, CyberKnifeRampVersion=None,
                                                     CyberKnifeAllowIncreasedPitchCorrection=None,
                                                     IsocenterData={'Position': {'x': x, 'y': y, 'z': z},
                                                                    'NameOfIsocenterToRef': "",
                                                                    'Name': 'iso_laser_vert',
                                                                    'Color': "98, 184, 234"}, Name=bs_name,
                                                     Description="",
                                                     GantryAngle=0, CouchRotationAngle=0, CouchPitchAngle=0,
                                                     CouchRollAngle=0,
                                                     CollimatorAngle=0)

                retval_0.SetBolus(BolusName="")
                beam_set.BeamMU = 0
            except:
                # le faisceau existe déjà.
                Exception("le faisceau existe déjà")

            # modification des paramètres du faisceau
            pitch = 0.380
            gantry_period = 20  # impossible de mettre moins
            try:
                plan.PlanOptimizations[0].OptimizationParameters.TreatmentSetupSettings[0].BeamSettings[
                    0].TomoPropertiesPerBeam.EditTomoBasedBeamOptimizationSettings(JawMode="Dynamic",
                                                                                   PitchTomoHelical=pitch,
                                                                                   PitchTomoDirect=None,
                                                                                   BackJawPosition=2.1,
                                                                                   FrontJawPosition=-2.1,
                                                                                   MaxDeliveryTime=None,
                                                                                   MaxGantryPeriod=gantry_period,
                                                                                   MaxDeliveryTimeFactor=None)
            except:
                print("No changes to save.")

        elif direction == 'FFS':
            try:
                beam_names = ['Ant', 'Post']
                angles = [0, 180]

                for beam_name, angle in zip(beam_names, angles):
                    # Création du faisceau antérieur
                    retval_0 = beam_set.CreatePhotonBeam(BeamQualityId="6", CyberKnifeCollimationType="Undefined",
                                                         CyberKnifeNodeSetName=None, CyberKnifeRampVersion=None,
                                                         CyberKnifeAllowIncreasedPitchCorrection=None,
                                                         IsocenterData={'Position': {'x': x, 'y': y, 'z': z},
                                                                        'NameOfIsocenterToRef': bs_name,
                                                                        'Name': bs_name,
                                                                        'Color': "98, 184, 234"}, Name=beam_name,
                                                         Description="", GantryAngle=angle, CouchRotationAngle=0,
                                                         CouchPitchAngle=0, CouchRollAngle=0, CollimatorAngle=0)
                    # lignes obligatoires
                    retval_0.SetBolus(BolusName="")
                    beam_set.Beams[beam_name].BeamMU = 0

                # Attribution du champ 5 et du pitch de 0.5
                pitch = 0.5
                max_delivery_factor = 1.3

                for num in range(len(beam_names)):
                    plan.PlanOptimizations[0].OptimizationParameters.TreatmentSetupSettings[0].BeamSettings[
                        num].TomoPropertiesPerBeam.EditTomoBasedBeamOptimizationSettings(JawMode="Dynamic",
                                                                                         PitchTomoHelical=None,
                                                                                         PitchTomoDirect=pitch,
                                                                                         BackJawPosition=2.1,
                                                                                         FrontJawPosition=-2.1,
                                                                                         MaxDeliveryTime=None,
                                                                                         MaxGantryPeriod=None,
                                                                                         MaxDeliveryTimeFactor=max_delivery_factor)

                # modification du point de specification de dose pour le plan FFS (sinon export impossible)
                obj_patient.create_dsp()
                xd, yd, zd = obj_patient.poi_DSP
                retval_0 = beam_set.CreateDoseSpecificationPoint(Name="DSP", Coordinates={'x': xd,
                                                                                          'y': yd,
                                                                                          'z': zd},
                                                                 VisualizationDiameter=1)
                for beam in beam_set.Beams:
                    beam.SetDoseSpecificationPoint(Name="DSP")

            except:
                print("No changes to save.")

        beam_set.SetDefaultDoseGrid(VoxelSize={'x': 0.5, 'y': 0.5, 'z': 0.5})

        ###############################################################################################################
        ###############################################################################################################
        ############################################   OPTIMISATION     ###############################################
        ###############################################################################################################
        ###############################################################################################################

        # Les paramètres d'optimisation sont inscrits dans un fichier csv dans le g commun. Le fichier est lu avec
        # pandas

        if direction == 'HFS':
            if not pediatrique:
                if total_dose > 800:
                    filename = "adultesHFS.csv"
                else:
                    filename = "2GyHFS.csv"
            else:
                filename = "pediatrique.csv"
        else:
            filename = "adultesFFS.csv"

        csv_path = os.path.join(obj_patient.directory, filename)
        df = pd.read_csv(csv_path, sep=';')

        for row in df.iloc:
            FunctionType = row.FunctionType
            RoiName = row.RoiName
            DoseLevel = row.DoseLevel * total_dose
            Weight = row.Weight
            IsRobust = row.IsRobust
            IsConstraint = row.IsConstraint
            HighDoseLevel = row.HighDoseLevel * total_dose

            # Utilisation de la fonction set_obj_function (définie hors classe)
            set_obj_function(plan, FunctionType, RoiName, DoseLevel, Weight, IsRobust, IsConstraint, HighDoseLevel)

        ###############################################################################################################
        ###############################################################################################################
        ############################################   Robustess        ###############################################
        ###############################################################################################################
        ###############################################################################################################

        decalage = 1.5  # cm
        plan.PlanOptimizations[0].OptimizationParameters.SaveRobustnessParameters(PositionUncertaintyAnterior=0,
                                                                                  PositionUncertaintyPosterior=0,
                                                                                  PositionUncertaintySuperior=0,
                                                                                  PositionUncertaintyInferior=0,
                                                                                  PositionUncertaintyLeft=decalage,
                                                                                  PositionUncertaintyRight=decalage,
                                                                                  DensityUncertainty=0,
                                                                                  PositionUncertaintySetting="Universal",
                                                                                  IndependentLeftRight=True,
                                                                                  IndependentAnteriorPosterior=True,
                                                                                  IndependentSuperiorInferior=True,
                                                                                  ComputeExactScenarioDoses=False,
                                                                                  NamesOfNonPlanningExaminations=[],
                                                                                  PatientGeometryUncertaintyType="PerTreatmentCourse",
                                                                                  PositionUncertaintyType="PerTreatmentCourse",
                                                                                  TreatmentCourseScenariosFactor=1000)

    print('~~~~ Saving case ~~~~')
    obj_patient.patient.Save()
