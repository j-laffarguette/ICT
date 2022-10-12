import os, sys
import random

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
from easygui import multenterbox
from datetime import datetime


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


class Patient:
    def __init__(self):

        # Administratif Raystation
        self.case = get_current("Case")
        self.patient = get_current("Patient")
        self.db = get_current("PatientDB")

        try:
            # juste pour différencier les scanners de toulouse des nôtres (pour le zéro scan à zéro et non pas à
            # hauteur table)
            self.doctor = self.case.Physician.Name
        except:
            print('Pas de docteur')

        # Récupération des informations des CT (nom + scan position HFS ou FFS)
        self.examinations = self.get_ct_list()  # dictionnaire {"exam_name":"FFS/HFS", ...}

        # Par défaut, à la première création de l'objet patient, le scanner Head First est pris en primary
        self.set_primary(self.examinations['HFS'])
        print('pour débug')

    def set_primary(self, exam_name):
        """Méthode très importante. Permet de définir un exam en 'primary et d'en récupérer les propriétés.
         Input : nom de l'examen à mettre en primary"""

        self.case.Examinations[exam_name].SetPrimary()
        self.examination = get_current("Examination")
        self.exam_name = self.examination.Name

        # Récupération de la position de tous les POI sur le scanner HFS
        self.jonction = self.get_point_coords('jonction')
        self.genoux = self.get_point_coords('genoux')
        self.abdomen = self.get_point_coords('abdomen')

        # Création du point zero (crane) situé à 0,0,hauteur table (valeur exprimée en mm convertie en cm)
        zero_scan = \
            self.case.Examinations[self.exam_name].GetStoredDicomTagValueForVerification(Group=0x0018, Element=0x1130)[
                'Table Height']

        if self.doctor == 'IZAR^FRANCOISE':  # si c'est le cas de toulouse, on prend zéro
            print("Il s'agit du cas test de toulouse, on utilise donc z = 0")
            self.zero_scan = 0
        else:
            print("Il ne s'agit pas du cas test de toulouse, on utilise donc z = HT")
            self.zero_scan = -float(zero_scan) / 10

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
        de ce point (ex: bille genoux) et qu'il aille jusqu'à ce point (ex: bille jonction)
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
            self.algebra(out_roi=roi_name, in_roiA=[roi_name], in_roiB=["External"], margeB=-0.3,
                         ResultOperation="Intersection", derive=False)

    def get_point_coords(self, point_name):
        """ Méthode permettant de récupérer les coordonnées d'un point si celui-ci existe. \n
        - Input: nom du point \n
        - Outupt: coordonnées du point (x,y,z) ou None"""
        point = self.case.PatientModel.StructureSets[self.exam_name].PoiGeometries[point_name]
        coords = point.Point
        if coords is not None:
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
            if exam.EquipmentInfo.Modality == 'CT':
                examinations[exam.PatientPosition] = exam.Name
        self.examinations = examinations
        return self.examinations

    def create_ROI(self, roi_name, roi_type='Ptv'):
        """Méthode utilisée pour créer des volumes de type PTV"""
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


if __name__ == '__main__':

    # todo: contour automatique du contour externe. Impossible pour l'heure à cause de l'absence de densités

    # Création de l'objet patient. Définit par défaut le scanner HFS en primary
    obj_patient = Patient()

    # do_it est mis sur False pour passer toute la partie Patient Modeling pour la programmation du script
    do_it = True
    if do_it:
        #########################################################################
        #########################################################################
        ##################### TRAVAIL DE PREPARATION ############################
        #########################################################################
        #########################################################################

        if False:
            # Première partie, réalisation du recalage rigide
            # Suppression pure et simple de tous les FOR registration déjà existants
            for reg in obj_patient.case.Registrations:
                FFOR = reg.FromFrameOfReference
                RFOR = reg.ToFrameOfReference
                obj_patient.case.RemoveFrameOfReferenceRegistration(
                    FloatingFrameOfReference=FFOR,
                    ReferenceFrameOfReference=RFOR)

            # Création du recalage entre le scanner FFS et le scanner HFS
            # 1. création du for reg
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

            # Deuxième partie : réalisation du recalage élastique "rigide". On utilise la méthode discard intensity avec la
            # vessie comme volume d'intérêt. Pour cela, la vessie doit être copiée d'un scanner à l'autre.

            obj_patient.case.PatientModel.CopyRoiGeometries(SourceExamination=obj_patient.examination,
                                                            TargetExaminationNames=[obj_patient.examinations['FFS']],
                                                            RoiNames=["Vessie"], ImageRegistrationNames=[],
                                                            TargetExaminationNamesToSkipAddedReg=[
                                                                obj_patient.examinations['FFS']])

            # Création du recalage élastique rigide
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

            # mapping des POI d'un scanner vers l'autre
            pois = ["jonction", "genoux", "abdomen"]

            obj_patient.case.MapPoiGeometriesDeformably(PoiGeometryNames=pois,
                                                        CreateNewPois=False,
                                                        StructureRegistrationGroupNames=[reg_name],
                                                        ReferenceExaminationNames=[obj_patient.examinations['HFS']],
                                                        TargetExaminationNames=[obj_patient.examinations['FFS']],
                                                        ReverseMapping=False, AbortWhenBadDisplacementField=False)

        # Travail sur les volumes

        # copie des structures de table du scanner HFS vers le FFS
        # todo: résoudre le problème qui fait que les structures de table sont coupées lors de la copy du HFS vers le FFS
        obj_patient.case.PatientModel.CopyRoiGeometries(SourceExamination=obj_patient.examination,
                                                        TargetExaminationNames=[obj_patient.examinations["FFS"]],
                                                        RoiNames=["Mousse", "Renfort interne POM 1.4", "Nylon 1.15",
                                                                  "Carbon Fiber",
                                                                  "Renfort 1.2", "Upper pallet",
                                                                  "Lower pallet Radixact"],
                                                        ImageRegistrationNames=[],
                                                        TargetExaminationNamesToSkipAddedReg=[
                                                            obj_patient.examinations["FFS"]])

        # Création des ROI
        ROI_LIST = ['PTV FFS', 'PTV_4', 'PTV_3', 'PTV_2', 'PTV_1', "PTV HFS"]
        # Vérification de l'existance des differentes roi dans le case
        resultat = [check_roi(obj_patient.case, roi) for roi in ROI_LIST]
        # Création des ROI dans le roi set si ROI non existantes
        [obj_patient.create_ROI(ROI_LIST[index]) for index, roi in enumerate(resultat) if roi == False]

        # ------------------------------------------------------------------------------------
        # On commence d'abord sur le scanner Head First puis on travaille sur le Feet First

        for direction in ['HFS', 'FFS']:
            # ct_name pour simplification dans la boucle for
            ct_name = obj_patient.examinations[direction]
            ct = obj_patient.case.Examinations[ct_name]

            # ct étudié mis en primary (très important, car re-définit self.exam_name et self.examination dans la classe
            obj_patient.set_primary(ct_name)

            # Retrait des trous dans externe
            obj_patient.case.PatientModel.StructureSets[ct_name].SimplifyContours(RoiNames=["External"],
                                                                                  RemoveHoles3D=True,
                                                                                  RemoveSmallContours=False,
                                                                                  AreaThreshold=None,
                                                                                  ReduceMaxNumberOfPointsInContours=False,
                                                                                  MaxNumberOfPoints=None,
                                                                                  CreateCopyOfRoi=False,
                                                                                  ResolveOverlappingContours=True)

            # créer poumon G+D
            obj_patient.algebra(out_roi="Poumons", in_roiA=["Poumon D", "Poumon G"], color="Aqua")

            # créer somme des reins
            obj_patient.algebra(out_roi="Reins", in_roiA=["Rein D", "Rein G"], color="Yellow")

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
            lim_inf,lim_sup = [lim.z for lim in obj_patient.case.PatientModel.StructureSets[ct_name].RoiGeometries["External"].GetBoundingBox()]

            # Réalisation du "PTV HFS", au-dessus du dernier PTV de la jonction, jusqu'au sommet du crane
            haut_PTV1 = y0 + 4
            from_to.append(["PTV HFS", lim_sup, haut_PTV1])

            if direction == 'HFS':
                # Création PTVpoumons = poumons - 1cm # seulement pour le scanner HFS!
                obj_patient.algebra(out_roi='PTVpoumons', in_roiA=["Poumons"], margeA=-1, Type="Ptv",color="Yellow")
            else :
                # Supprimer la dérivation
                obj_patient.case.PatientModel.RegionsOfInterest['PTVpoumons'].DeleteExpression()

            # Réalisation du PTV FFS : en dessous du PTV E
            from_to.append(['PTV FFS', y0-4, lim_inf])

            # Une fois que toutes les bornes sont déterminées, on crée les PTVS
            for roi, f, t in from_to:
                print(f'Création du PTV suivant : {roi} entre y = {str(f)} et y = {str(t)} ')
                obj_patient.create_cylinder_ptv(roi, f, t)
                print('-> OK!')

            # PTV HFS - PTVPoumons = PTVB
            obj_patient.algebra(out_roi='PTV HFS', in_roiA=["PTV HFS"], in_roiB=["PTVpoumons"],
                                ResultOperation="Subtraction", derive=False)

            # Dernière étape : soustraction des volumes deux à deux pour éviter le chevauchement des PTVS
            ROI_LIST = ['PTV FFS', 'PTV_4', 'PTV_3', 'PTV_2', 'PTV_1', "PTV HFS"]
            for i in range(len(ROI_LIST) - 1):
                obj_patient.algebra(out_roi=ROI_LIST[i], in_roiA=[ROI_LIST[i]], in_roiB=[ROI_LIST[i + 1]],
                                    ResultOperation="Subtraction", derive=False)

            # simplification des volumes pour éviter overlaps
            roi_list = [roi.Name for roi in obj_patient.case.PatientModel.RegionsOfInterest if roi.Type == "Support"]

            obj_patient.case.PatientModel.StructureSets[obj_patient.exam_name].SimplifyContours(
                RoiNames=roi_list,
                RemoveHoles3D=False, RemoveSmallContours=False, AreaThreshold=None,
                ReduceMaxNumberOfPointsInContours=False, MaxNumberOfPoints=None, CreateCopyOfRoi=False,
                ResolveOverlappingContours=True)

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


    #########################################################################
    #########################################################################
    ################ TRAVAIL SUR LE PLAN DE TRAITEMENT ######################
    #########################################################################
    #########################################################################

    # Création des deux plans
    date = datetime.today().strftime('%Y%m%d')
    PLAN_NAMES = [f"HFS_{date}", f"FFS_{date}"]
    BS_NAMES = ["bs_HFS", "bs_FFS"]
    PATIENT_POSITIONS = ["HeadFirstSupine", "FeetFirstSupine"]
    machine_name = "Radixact1"

    ########################################################
    # définition de la prescription
    verif = 0

    msg = "Entrer les informations de la prescription"
    title = "Prescription"
    fields = ("Nombre de fractions", "Dose totale (Gy)")

    easygui = True

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

    for ind, direction in enumerate(['HFS', 'FFS']):

        # CT étudié en primary
        obj_patient.set_primary(obj_patient.examinations[direction])

        # Attribution de l'IVDT qui va bien
        if not obj_patient.examination.EquipmentInfo.ImagingSystemReference:
            obj_patient.examination.EquipmentInfo.SetImagingSystemReference(ImagingSystemName="SIEMENS_RT_mDD")
            obj_patient.patient.Save()
        elif obj_patient.examination.EquipmentInfo.ImagingSystemReference.get('ImagingSystemName') == 'SIEMENS_RT_mDD':
            obj_patient.examination.EquipmentInfo.SetImagingSystemReference(ImagingSystemName="SIEMENS_RT_mDD")
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

        obj_patient.patient.Save()
        obj_patient.case.TreatmentPlans[plan_name].SetCurrent()

        # Le plan créé est mis en primary pour simplifier
        plan = get_current('Plan')

        # Sauvegarde du plan (obligatoire)
        obj_patient.patient.Save()
        # Le bs créé est mis en primary pour simplifier
        plan.BeamSets[bs_name].SetCurrent()

        beam_set = get_current("BeamSet")

        # Optimization settings
        plan.PlanOptimizations[0].OptimizationParameters.Algorithm.OptimalityTolerance = 1e-7
        plan.PlanOptimizations[0].OptimizationParameters.Algorithm.MaxNumberOfIterations = 40
        plan.PlanOptimizations[0].OptimizationParameters.DoseCalculation.ComputeFinalDose = True

        # if ind = 0 = HFS la prescription est attribuée au PTV haut, sinon au PTV_E pour le plan FFS
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
        print('')

        ########################################################
        # définition des points de ref

        # ------------------------------------------------------------------------------------
        # Positionnement du localization point HFS -> LASER ROUGE
        # Le loc point du traitement HFS est en Inf/Sup sur les billes thorax, à zéro en G/D
        # Il est situé en ant/post sur les billes cranes = zero scan = table height (voir __init__ de Patient())

        # modification du type de tous les points en Undefined au cazou
        for poi in obj_patient.case.PatientModel.PointsOfInterest:
            poi.Type = 'Undefined'

        laser_rouge_HFS = "laser rouge HFS"

        # On positionne les lasers rouges et vert sur le point abdomen mais le rouge est décalé sur la bille genoux en GD
        coords_laser_rouges_HFS = (obj_patient.genoux[0], obj_patient.zero_scan, obj_patient.abdomen[2])

        if direction == 'HFS':  # Plan HFS
            coords_laser_vert = (0, obj_patient.zero_scan, obj_patient.abdomen[2])
        elif direction == 'FFS':  # Plan FFS
            # Pour le plan FFS, le laser vert est positionné en AP au centre du volume PTV_E afin que le patient
            # soit bien centré en hauteur et n'ait pas les genoux qui dépassent du cadre
            AP_iso_FFS = obj_patient.case.PatientModel.StructureSets[obj_patient.examinations['FFS']].RoiGeometries[
                'PTV FFS'].GetCenterOfRoi().y
            coords_laser_vert = (0, AP_iso_FFS, obj_patient.abdomen[2])

        # création du laser rouge
        obj_patient.create_poi(laser_rouge_HFS, coords_laser_rouges_HFS, color='Red')
        # laser rouge devient localization point
        obj_patient.case.PatientModel.PointsOfInterest[laser_rouge_HFS].Type = "LocalizationPoint"
        # création du laser vert
        obj_patient.create_poi('laser vert HFS', coords_laser_vert, color='Green')

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
                beam_names = ['Ant', 'Droit', 'Post', 'Gauche']
                angles = [0, 90, 180, 270]

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
                    # lignes semble-t-il obligatoires
                    retval_0.SetBolus(BolusName="")
                    beam_set.Beams[beam_name].BeamMU = 0

                # Attribution du champ 5 et du pitch de 0.5
                pitch = 0.5
                for num in range(4):
                    plan.PlanOptimizations[0].OptimizationParameters.TreatmentSetupSettings[0].BeamSettings[
                        num].TomoPropertiesPerBeam.EditTomoBasedBeamOptimizationSettings(JawMode="Dynamic",
                                                                                         PitchTomoHelical=None,
                                                                                         PitchTomoDirect=pitch,
                                                                                         BackJawPosition=2.1,
                                                                                         FrontJawPosition=-2.1,
                                                                                         MaxDeliveryTime=None,
                                                                                         MaxGantryPeriod=None,
                                                                                         MaxDeliveryTimeFactor=1.15)

            except:
                print("No changes to save.")

        beam_set.SetDefaultDoseGrid(VoxelSize={'x': 0.5, 'y': 0.5, 'z': 0.5})

        ################################################################################################################
        ################################################################################################################
        ############################################   OPTIMISATION     ################################################
        ################################################################################################################
        ################################################################################################################

        # ////////////////////////////////////////////////////////////////////////////////////////////////////////////
        # Plan HFS
        if ind == 0:  # Plan HFS

            # Ici on associe à chaque PTV un certain niveau par rapport à la dose de prescription. Travail avec numpy
            # pour pouvoir utiliser la méthode np.where, associant la dose à la roi
            ptv_list = np.array(['PTV HFS', 'PTV_1', 'PTV_2', 'PTV_3', 'PTV_4', 'PTVpoumons'])

            # Définition de la contrainte aux poumons. Si la dose prescrite est supérieure à 8Gy, on applique 7.5 Gy en
            # contrainte aux poumons (vu avec SP, ça marche mieux que de mettre 8)
            if total_dose > 800:
                dose_poumons = 750
            # Sinon on garde la même dose aux poumons qu'au reste.
            else:
                dose_poumons = total_dose

            objectives_list = [total_dose, 0.9 * total_dose, 0.63 * total_dose, 0.37 * total_dose, 0.16 * total_dose,
                               dose_poumons]

            weight_list = [10, 1, 1, 1, 1, 1]

            # ajout d'un max EUD au poumon avec fort poids
            contrainte_poumon = 800
            weight_poumon = 100

            # ajout d'un max EUD au reins avec fort poids
            contrainte_reins = 1000
            weight_reins = 100

            # ajout d'un "dose fall-off" sur les derniers volumes d'opt
            obj_DFO = 0.16 * total_dose
            weight_DFO = 50

            robustesse = False

        # ////////////////////////////////////////////////////////////////////////////////////////////////////////////
        # Plan FFS
        elif ind == 1:  # Plan FFS
            ptv_list = np.array(
                ['PTV FFS', 'PTV_4', 'PTV_3', 'PTV_2', 'PTV_1'])

            objectives_list = [total_dose, 0.9 * total_dose, 0.63 * total_dose, 0.37 * total_dose,
                               0.16 * total_dose]

            weight_list = [1, 1, 1, 1, 1]
            robustesse = True

        # Ici, on crée des fonctions d'optimisation pour chacun des PTVS. Initialement, ces fonctions sont vides,
        # elles sont remplies ensuite
        for roi in ptv_list:

            liste = ["PTV HFS", "PTVpoumons", "PTV FFS"]
            if any(roi in s for s in liste):  # On ne met l'uniform dose que sur ces volumes là
                plan.PlanOptimizations[0].AddOptimizationFunction(FunctionType="UniformDose", RoiName=str(roi),
                                                                  IsConstraint=False,
                                                                  RestrictAllBeamsIndividually=False,
                                                                  RestrictToBeam=None, IsRobust=False,
                                                                  RestrictToBeamSet=None, UseRbeDose=False)
            plan.PlanOptimizations[0].AddOptimizationFunction(FunctionType="MinDose", RoiName=str(roi),
                                                              IsConstraint=False, RestrictAllBeamsIndividually=False,
                                                              RestrictToBeam=None, IsRobust=robustesse,
                                                              RestrictToBeamSet=None, UseRbeDose=False)
            plan.PlanOptimizations[0].AddOptimizationFunction(FunctionType="MaxDose", RoiName=str(roi),
                                                              IsConstraint=False, RestrictAllBeamsIndividually=False,
                                                              RestrictToBeam=None, IsRobust=False,
                                                              RestrictToBeamSet=None, UseRbeDose=False)

            # On attribue la dose totale à tous les objectifs
            for dose in plan.PlanOptimizations[0].Objective.ConstituentFunctions:
                # On regarde dans le xième objectif à quelle roi il est associé
                associated_roi = dose.OfDoseGridRoi.OfRoiGeometry.OfRoi.Name
                print(associated_roi)
                # On recherche dans les np.array créée ci -dessus à quel objectif cela correspond
                indice = np.where(ptv_list == associated_roi)[0][0]
                print(indice)
                # On attribue la dose correspondante
                dose.DoseFunctionParameters.DoseLevel = objectives_list[indice]
                dose.DoseFunctionParameters.Weight = weight_list[indice]

            obj_patient.patient.Save()

        # Ajout du dose fall-of sur le volume opt_jonction (OAR_name)
        OAR_name = "opt_jonction"
        dose = plan.PlanOptimizations[0].AddOptimizationFunction(FunctionType="DoseFallOff", RoiName=OAR_name,
                                                                 IsConstraint=False,
                                                                 RestrictAllBeamsIndividually=False,
                                                                 RestrictToBeam=None, IsRobust=False,
                                                                 RestrictToBeamSet=None, UseRbeDose=False)
        dose.DoseFunctionParameters.HighDoseLevel = obj_DFO
        dose.DoseFunctionParameters.LowDoseDistance = 0.5
        dose.DoseFunctionParameters.Weight = weight_DFO
        dose.DoseFunctionParameters.AdaptToTargetDoseLevels = True

        # Finalement à la toute fin, on ajoute le max EUD sur le poumon pour le plan HFS
        try:
            if ind == 0:  # Plan HFS
                if total_dose >= 800:
                    dose = plan.PlanOptimizations[0].AddOptimizationFunction(FunctionType="MaxEud", RoiName="Poumons",
                                                                             IsConstraint=False,
                                                                             RestrictAllBeamsIndividually=False,
                                                                             RestrictToBeam=None, IsRobust=False,
                                                                             RestrictToBeamSet=None, UseRbeDose=False)
                    dose.DoseFunctionParameters.DoseLevel = contrainte_poumon
                    dose.DoseFunctionParameters.Weight = weight_poumon

                if total_dose >= 1000:
                    for rein in ['Rein D', 'Rein G']:
                        dose = plan.PlanOptimizations[0].AddOptimizationFunction(FunctionType="MaxEud", RoiName=rein,
                                                                                 IsConstraint=False,
                                                                                 RestrictAllBeamsIndividually=False,
                                                                                 RestrictToBeam=None, IsRobust=False,
                                                                                 RestrictToBeamSet=None,
                                                                                 UseRbeDose=False)
                        dose.DoseFunctionParameters.DoseLevel = contrainte_reins
                        dose.DoseFunctionParameters.Weight = weight_reins
        except:
            print("ça n'a pas marché le truc du poumon, là")

        # Robustesse
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
