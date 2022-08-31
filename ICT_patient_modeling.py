import os, sys
import random

pid_file_path = os.path.join(os.environ.get('userprofile'), 'AppData', 'Local', 'Temp', 'raystation.pid')
with open(pid_file_path) as f:
    os.environ['RAYSTATION_PID'] = f.read()
script_client_path = r'C:\Program Files\RaySearch Laboratories\RayStation 11B\ScriptClient'
sys.path.append(script_client_path)

from connect import *
import numpy as np


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
        self.examination = get_current("Examination")
        self.exam_name = self.examination.Name
        self.patient = get_current("Patient")
        self.db = get_current("PatientDB")

        # Récupération des informations des CT (nom + scan position HFS ou FFS)
        self.examinations = []  # dictionnaire {"exam_name":"FFS/HFS", ...}
        self.get_ct_list()  # méthode utilisée pour récupérer les données
        self.set_primary(self.examinations['HFS'])



        self.roi_list = []
        self.patient_id = self.patient.PatientID




        # Récupération de la position de tous les POI sur le scanner HFS

        self.jonction = self.get_point_coords('jonction')
        self.pubis = self.get_point_coords('pubis')
        self.genoux = self.get_point_coords('genoux')
        self.cou = self.get_point_coords('cou')
        self.thorax = self.get_point_coords('thorax')
        self.abdomen = self.get_point_coords('abdomen')  # todo: ajouter à la procédure
        # Création du point zero (crane) situé à 0,0,hauteur table (valeur exprimée en mm convertie en cm)
        zero_scan = \
        self.case.Examinations[self.exam_name].GetStoredDicomTagValueForVerification(Group=0x0018, Element=0x1130)[
            'Table Height']
        self.zero_scan = -float(zero_scan) / 10

        # todo création des points à la toute fin pour simplification
        # self.create_poi('crane', self.examinations['HFS'], (0, self.zero_scanHFS, 0))
        # self.create_poi('jonctionFFS', self.examinations['FFS'], (0, self.zero_scanFFS, 0))
        # print('')

    def get_zero_scan(self, scan_direction):
        zero_scan = \
        self.case.Examinations[self.examinations[scan_direction]].GetStoredDicomTagValueForVerification(Group=0x0018,
                                                                                                        Element=0x1130)[
            'Table Height']
        return -float(zero_scan) / 10

    def create_cylinder_ptv(self, roi_name, y_cranial, y_caudal):
        x0, z0 = 0, self.zero_scan
        long = y_cranial - y_caudal  # longueur du volume
        center = y_cranial - long / 2  # centre du volume
        obj_patient.cylinder(roi_name, (x0, z0, center), longueur=abs(long))


    def cylinder(self, roi_name, coords, longueur=2, retraction=True):
        x, y,z = coords
        self.case.PatientModel.RegionsOfInterest[roi_name].CreateCylinderGeometry(Radius=30,
                                                                                  Axis={'x': 0, 'y': 0, 'z': 1},
                                                                                  Length=longueur,
                                                                                  Examination=self.examination,
                                                                                  Center={'x': x,
                                                                                          'y': y,
                                                                                          'z': z},
                                                                                  Representation="TriangleMesh",
                                                                                  VoxelSize=None)

        if retraction:
            self.retraction(roi_name)


    def get_point_coords(self, point_name):
        point = self.case.PatientModel.StructureSets[self.exam_name].PoiGeometries[point_name]
        coords = point.Point
        if coords is not None:
            return coords.x, coords.y, coords.z
        else:
            return None

    def algebra_soustraction(self, roi_a, roi_b):
        # Simple soustraction entre deux volumes. Ecrase le volume de départ. Utilisé pour PTV_B-poumons = PTV_B
        self.case.PatientModel.RegionsOfInterest[roi_a].CreateAlgebraGeometry(
            Examination=self.examination, Algorithm="Auto",
            ExpressionA={'Operation': "Union",
                         'SourceRoiNames': [roi_a],
                         'MarginSettings': {
                             'Type': "Expand",
                             'Superior': 0,
                             'Inferior': 0,
                             'Anterior': 0,
                             'Posterior': 0, 'Right': 0,
                             'Left': 0}},
            ExpressionB={'Operation': "Union",
                         'SourceRoiNames': [
                             roi_b],
                         'MarginSettings': {
                             'Type': "Expand",
                             'Superior': 0,
                             'Inferior': 0,
                             'Anterior': 0,
                             'Posterior': 0, 'Right': 0,
                             'Left': 0}},
            ResultOperation="Subtraction",
            ResultMarginSettings={'Type': "Expand",
                                  'Superior': 0,
                                  'Inferior': 0,
                                  'Anterior': 0,
                                  'Posterior': 0,
                                  'Right': 0, 'Left': 0})

    def set_primary(self, exam_name):
        self.case.Examinations[exam_name].SetPrimary()
        self.examination = get_current("Examination")
        self.exam_name = self.examination.Name
        # Récupération de la position de tous les POI sur le scanner HFS

        self.jonction = self.get_point_coords('jonction')
        self.pubis = self.get_point_coords('pubis')
        self.genoux = self.get_point_coords('genoux')
        self.cou = self.get_point_coords('cou')
        self.thorax = self.get_point_coords('thorax')
        self.abdomen = self.get_point_coords('abdomen')  # todo: ajouter à la procédure
        # Création du point zero (crane) situé à 0,0,hauteur table (valeur exprimée en mm convertie en cm)
        zero_scan = \
        self.case.Examinations[self.exam_name].GetStoredDicomTagValueForVerification(Group=0x0018, Element=0x1130)[
            'Table Height']
        self.zero_scan = -float(zero_scan) / 10

    def get_ct_list(self):
        examinations = {}
        for exam in self.case.Examinations:
            if exam.EquipmentInfo.Modality == 'CT':
                examinations[exam.PatientPosition] = exam.Name
        self.examinations = examinations
        return self.examinations

    def get_irm_list(self):
        name, modality = [], []
        for exam in self.case.Examinations:
            if exam.EquipmentInfo.Modality == 'MR':
                name.append(exam.Name)
                modality.append(exam.GetProtocolName())
        self.irm_names = list(zip(name, modality))
        return self.irm_names

    def get_dixon_name(self):
        irm_list = self.get_irm_list()
        for irm in irm_list:
            irm_name = irm[0]
            irm_modality = irm[1]
            if "DIXON" in irm_modality and "RESPI LIBRE GADO TARDIF" in irm_modality:
                self.dixon_name.append(irm_name)
                print(irm_name)
        return self.dixon_name

    def get_roi_list(self):
        rois = self.case.PatientModel.RegionsOfInterest
        self.roi_list = [roi.Name for roi in rois]
        return self.roi_list

    def get_description(self, imageName):
        description = self.case.Examinations[imageName].GetStoredDicomTagValueForVerification(Group=0x0008,
                                                                                              Element=0x103E)
        description = description.__getitem__("SeriesDescription")
        return description

    def create_ROI(self, roi_name, roi_type='Ptv'):
        color = ["#" + ''.join([random.choice('ABCDEF0123456789') for i in range(6)])][0]
        self.case.PatientModel.CreateRoi(Name=roi_name, Color=color, Type=roi_type)

    def cylinder(self, roi_name, coords, longueur=2, retraction=True):
        x, y,z = coords
        self.case.PatientModel.RegionsOfInterest[roi_name].CreateCylinderGeometry(Radius=30,
                                                                                  Axis={'x': 0, 'y': 0, 'z': 1},
                                                                                  Length=longueur,
                                                                                  Examination=self.examination,
                                                                                  Center={'x': x,
                                                                                          'y': y,
                                                                                          'z': z},
                                                                                  Representation="TriangleMesh",
                                                                                  VoxelSize=None)

        if retraction:
            self.retraction(roi_name)

    def retraction(self, roi_name):
        self.case.PatientModel.RegionsOfInterest[roi_name].CreateAlgebraGeometry(Examination=self.examination,
                                                                            Algorithm="Auto",
                                                                            ExpressionA={'Operation': "Union",
                                                                                         'SourceRoiNames': [roi_name],
                                                                                         'MarginSettings': {
                                                                                             'Type': "Expand",
                                                                                             'Superior': 0,
                                                                                             'Inferior': 0,
                                                                                             'Anterior': 0,
                                                                                             'Posterior': 0, 'Right': 0,
                                                                                             'Left': 0}},
                                                                            ExpressionB={'Operation': "Union",
                                                                                         'SourceRoiNames': ["External"],
                                                                                         'MarginSettings': {
                                                                                             'Type': "Contract",
                                                                                             'Superior': 0.3,
                                                                                             'Inferior': 0.3,
                                                                                             'Anterior': 0.3,
                                                                                             'Posterior': 0.3,
                                                                                             'Right': 0.3,
                                                                                             'Left': 0.3}},
                                                                            ResultOperation="Intersection",
                                                                            ResultMarginSettings={'Type': "Expand",
                                                                                                  'Superior': 0,
                                                                                                  'Inferior': 0,
                                                                                                  'Anterior': 0,
                                                                                                  'Posterior': 0,
                                                                                                  'Right': 0,
                                                                                                  'Left': 0})

    def generate_poumons(self):
        if not check_roi(self.case, 'Poumons'):
            self.case.PatientModel.CreateRoi(Name="Poumons", Color="Aqua", Type="Organ", TissueName=None,
                                             RbeCellTypeName=None, RoiMaterial=None)

        retval_0 = self.case.PatientModel.RegionsOfInterest['Poumons'].SetAlgebraExpression(
            ExpressionA={'Operation': "Union", 'SourceRoiNames': ["Poumon D", "Poumon G"],
                         'MarginSettings': {'Type': "Expand", 'Superior': 0, 'Inferior': 0, 'Anterior': 0,
                                            'Posterior': 0, 'Right': 0, 'Left': 0}},
            ExpressionB={'Operation': "Union", 'SourceRoiNames': [],
                         'MarginSettings': {'Type': "Expand", 'Superior': 0, 'Inferior': 0, 'Anterior': 0,
                                            'Posterior': 0, 'Right': 0, 'Left': 0}}, ResultOperation="None",
            ResultMarginSettings={'Type': "Expand", 'Superior': 0, 'Inferior': 0, 'Anterior': 0, 'Posterior': 0,
                                  'Right': 0, 'Left': 0})
        retval_0.UpdateDerivedGeometry(Examination=self.examination, Algorithm="Auto")

    def generate_ptv_poumons(self):
        if not check_roi(self.case, 'PTVpoumons'):
            self.create_ROI('PTVpoumons')
        retval_1 = self.case.PatientModel.RegionsOfInterest['PTVpoumons'].SetAlgebraExpression(
            ExpressionA={'Operation': "Union", 'SourceRoiNames': ["Poumons"],
                         'MarginSettings': {'Type': "Contract", 'Superior': 1, 'Inferior': 1, 'Anterior': 1,
                                            'Posterior': 1,
                                            'Right': 1, 'Left': 1}},
            ExpressionB={'Operation': "Union", 'SourceRoiNames': [],
                         'MarginSettings': {'Type': "Expand", 'Superior': 0, 'Inferior': 0, 'Anterior': 0,
                                            'Posterior': 0,
                                            'Right': 0, 'Left': 0}}, ResultOperation="None",
            ResultMarginSettings={'Type': "Expand", 'Superior': 0, 'Inferior': 0, 'Anterior': 0, 'Posterior': 0,
                                  'Right': 0,
                                  'Left': 0})

        retval_1.UpdateDerivedGeometry(Examination=self.examination, Algorithm="Auto")

    def create_poi(self, poi_name, examination, coords, color="128, 128, 255"):
        x, y, z = coords
        poi_list = [poi.Name for poi in self.case.PatientModel.PointsOfInterest]
        if not poi_name in poi_list:
            retval_0 = self.case.PatientModel.CreatePoi(Examination=examination,
                                                        Point={'x': x, 'y': y, 'z': z},
                                                        Name=poi_name, Color=color,
                                                        VisualizationDiameter=1, Type="Undefined")
        else:
            self.case.PatientModel.StructureSets[examination].PoiGeometries[poi_name].Point = {
                'x': x, 'y': y, 'z': z}


if __name__ == '__main__':

    # todo: contour automatique du contour externe. Impossible pour l'heure à cause de l'absence de densités
    # todo: faire soustraction PTVD6 - PTVD5

    obj_patient = Patient()

    #########################################################################
    #########################################################################
    ##################### TRAVAIL DE PREPARATION ############################
    #########################################################################
    #########################################################################

    # # Suppression pure et simple de tous les FOR registration déjà existants
    # for reg in obj_patient.case.Registrations:
    #     FFOR = reg.FromFrameOfReference
    #     RFOR = reg.ToFrameOfReference
    #     obj_patient.case.RemoveFrameOfReferenceRegistration(
    #         FloatingFrameOfReference=FFOR,
    #         ReferenceFrameOfReference=RFOR)
    #
    # # Création du recalage entre le scanner FFS et le scanner HFS
    # obj_patient.case.CreateNamedIdentityFrameOfReferenceRegistration(
    #     FromExaminationName=obj_patient.examinations['HFS'], ToExaminationName=obj_patient.examinations['FFS'],
    #     RegistrationName="HFS to FFS", Description=None)
    #
    # obj_patient.case.ComputeGrayLevelBasedRigidRegistration(FloatingExaminationName=obj_patient.examinations['HFS'],
    #                                                         ReferenceExaminationName=obj_patient.examinations['FFS'],
    #                                                         UseOnlyTranslations=True, HighWeightOnBones=True,
    #                                                         InitializeImages=True, FocusRoisNames=[],
    #                                                         RegistrationName=None)

    # Travail sur les volumes

    # Création des ROI
    ROI_LIST = ['PTV_D1', 'PTV_D2', 'PTV_D3', 'PTV_D4', 'PTV_D5', 'PTV_D6', 'PTV_BOLUS', 'PTV_A', 'PTV_B', 'PTV_C',
                'PTV_E']
    # Vérification de l'existance des differentes roi dans le case
    resultat = [check_roi(obj_patient.case, roi) for roi in ROI_LIST]
    # Création des ROI dans le roi set si ROI non existante
    [obj_patient.create_ROI(ROI_LIST[index]) for index, roi in enumerate(resultat) if roi == False]

    # ------------------------------------------------------------------------------------
    # On commence d'abord sur le scanner Head First puis on travaille sur le Feet First


    for direction in ['FFS']:
        # ct_name pour simplification dans la boucle for
        ct_name = obj_patient.examinations[direction]
        ct = obj_patient.case.Examinations[ct_name]

        # ct étudié mis en primary (très important, car redéfinit self.exam_name et self.examination dans la classe
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
        obj_patient.generate_poumons()

        # Création PTVpoumons = poumons - 1cm # seulement pour le scanner HFS!
        if direction == 'HFS':
            obj_patient.generate_ptv_poumons()

        # ------------------------------------------------------------------------------------
        # Création des cylindres

        # récupération du point "jonction", positionné sur la bille au niveau de la jonction
        _, _, y0 = obj_patient.jonction  # y de ref est au niveau de la jonction
        x0, z0 = 0, obj_patient.zero_scan

        # Réalisation des PTV 2 à 5, entourant la bille jonction
        for index, roi in enumerate(['PTV_D2', 'PTV_D3', 'PTV_D4', 'PTV_D5']):
            # Réalisation de plusieurs cylindres espacés de 2 cm en commençant de sorte que la bille soit au niveau
            # de la premiere coupe du PTV D3
            y = y0 + 3 - 2 * index
            # Création du cylindre + retraction à la peau
            obj_patient.cylinder(roi, (x0, z0, y))

        from_to = []

        # Réalisation du PTV D6 : de la bille genoux jusqu'au bas du PTV D5
        bas_PTV5 = y0 - 4  # 4 cm en dessous de la bille jonction
        _, _, y_genoux = obj_patient.genoux

        # Réalisation du PTV D1 : de la bille pubis jusqu'au haut du PTV D2
        haut_PTVD2 = y0 + 4  # 4 cm au dessus de la bille jonction
        _, _, y_pubis = obj_patient.pubis

        from_to.append(['PTV_D6', y_genoux, bas_PTV5])
        from_to.append(['PTV_D1', haut_PTVD2, y_pubis])

        if direction == 'HFS':
            # Réalisation du PTV C : au dessus du PTV D1 jusqu'à 3 cm en dessous des poumons
            sous_poumons = \
                obj_patient.case.PatientModel.StructureSets[ct_name].RoiGeometries["Poumons"].GetBoundingBox()[
                    0].z  # relou, raystation intervertit y et z
            sous_poumons -= 3  # 3cm en dessous du plus bas du poumon

            # Réalisation du PTV B : au dessus du PTV C jusqu'à la bille cou
            # Edit 18/08/2022 -> soustraction des poumons
            _, _, y_cou = obj_patient.cou

            # TODO VOIR SI Nécessaire de mettre ça a la fin
            obj_patient.algebra_soustraction('PTV_B', "PTVpoumons")

            # Réalisation du PTV A : au dessus du PTV B
            tout_en_haut = y_cou + 50


            from_to.append(['PTV_C', sous_poumons, y_pubis])
            from_to.append(['PTV_B', y_cou, sous_poumons])
            from_to.append(['PTV_A',tout_en_haut, y_cou])


        elif direction == 'FFS':
            # Réalisation du PTV E : en dessous du PTV E
            tout_en_bas = y_genoux - 100

            from_to.append(['PTV_E', y_genoux, tout_en_bas])


        print(from_to)
        for roi, f, t in from_to:
            obj_patient.create_cylinder_ptv(roi, f, t)

        print('')
    # #########################################################################
    # #########################################################################
    # ################ TRAVAIL SUR LE PLAN DE TRAITEMENT ######################
    # #########################################################################
    # #########################################################################
    #
    # if not obj_patient.examination.EquipmentInfo.ImagingSystemReference:
    #     obj_patient.examination.EquipmentInfo.SetImagingSystemReference(ImagingSystemName="SIEMENS_RT_mDD")
    #     obj_patient.patient.Save()
    # elif obj_patient.examination.EquipmentInfo.ImagingSystemReference.get('ImagingSystemName') == 'SIEMENS_RT_mDD':
    #     obj_patient.examination.EquipmentInfo.SetImagingSystemReference(ImagingSystemName="SIEMENS_RT_mDD")
    #     obj_patient.patient.Save()
    #
    # # Création des deux plans
    # PLAN_NAMES = ["Plan ICT Head First", "Plan ICT Feet First"]
    # BS_NAMES = ["bs_HFS", "bs_FFS"]
    # PATIENT_POSITIONS = ["HeadFirstSupine", "FeetFirstSupine"]
    #
    # for ind in range(2):
    #
    #     # Attribution des noms de plans et de beamsets
    #     plan_name = PLAN_NAMES[ind]
    #     bs_name = BS_NAMES[ind]
    #
    #     # vérification de la non-existence du plan avant création
    #     list_of_plans = [plan.Name for plan in obj_patient.case.TreatmentPlans]
    #     if plan_name not in list_of_plans:
    #         # création du plan
    #         retval_0 = obj_patient.case.AddNewPlan(PlanName=plan_name, PlannedBy="", Comment="",
    #                                                ExaminationName=obj_patient.exam_name,
    #                                                IsMedicalOncologyPlan=False, AllowDuplicateNames=False)
    #
    #     # vérification de la non-existence du bs avant création
    #     list_of_bs = [bs.DicomPlanLabel for bs in obj_patient.case.TreatmentPlans[plan_name].BeamSets]
    #     if bs_name not in list_of_bs:
    #         # création du bs
    #         obj_patient.case.TreatmentPlans[plan_name].AddNewBeamSet(Name=bs_name,
    #                                                                  ExaminationName=obj_patient.exam_name,
    #                                                                  MachineName="Radixact1",
    #                                                                  Modality="Photons",
    #                                                                  TreatmentTechnique="TomoHelical",
    #                                                                  PatientPosition=PATIENT_POSITIONS[ind],
    #                                                                  NumberOfFractions=6, CreateSetupBeams=False,
    #                                                                  UseLocalizationPointAsSetupIsocenter=False,
    #                                                                  UseUserSelectedIsocenterSetupIsocenter=False,
    #                                                                  Comment="", RbeModelName=None,
    #                                                                  EnableDynamicTrackingForVero=False,
    #                                                                  NewDoseSpecificationPointNames=[],
    #                                                                  NewDoseSpecificationPoints=[],
    #                                                                  MotionSynchronizationTechniqueSettings={
    #                                                                      'DisplayName': None,
    #                                                                      'MotionSynchronizationSettings': None,
    #                                                                      'RespiratoryIntervalTime': None,
    #                                                                      'RespiratoryPhaseGatingDutyCycleTimePercentage': None,
    #                                                                      'MotionSynchronizationTechniqueType': "Undefined"},
    #                                                                  Custom=None, ToleranceTableLabel=None)
    #     obj_patient.patient.Save()
    #     obj_patient.case.TreatmentPlans[plan_name].SetCurrent()
    #     plan = get_current('Plan')
    #
    #     ########################################################
    #     # définition de la prescription
    #     fraction_dose = 200
    #     number_of_fractions = 6
    #     total_dose = fraction_dose * number_of_fractions
    #
    #     obj_patient.patient.Save()
    #     plan.BeamSets[bs_name].SetCurrent()
    #
    #     beam_set = get_current("BeamSet")
    #
    #     # if ind = 1 = HFS -> PTV A, sinon PTV_E
    #     if ind == 0:
    #         prescription_roi = 'PTV_A'
    #     else:
    #         prescription_roi = 'PTV_E'
    #
    #     # Si la prescription existe déjà, ne pas la recréer
    #     if not beam_set.Prescription.PrescriptionDoseReferences:
    #         beam_set.AddRoiPrescriptionDoseReference(RoiName=prescription_roi, DoseVolume=0,
    #                                                  PrescriptionType="MedianDose",
    #                                                  DoseValue=total_dose, RelativePrescriptionLevel=1)
    #     elif not beam_set.Prescription.PrescriptionDoseReferences[0].OnStructure.Name == prescription_roi:
    #         beam_set.AddRoiPrescriptionDoseReference(RoiName=prescription_roi, DoseVolume=0,
    #                                                  PrescriptionType="MedianDose",
    #                                                  DoseValue=total_dose, RelativePrescriptionLevel=1)
    #
    #     beam_set.SetAutoScaleToPrimaryPrescription(AutoScale=False)
    #     print('')
    #
    #     ########################################################
    #     # définition des points de ref
    #
    #     # ------------------------------------------------------------------------------------
    #     # Positionnement du localization point HFS -> LASER ROUGE
    #     # Le loc point du traitement HFS est en Inf/Sup sur les billes thorax, à zéro en G/D
    #     # Il est situé en ant/post sur les billes cranes = zero scan = table height (voir __init__ de Patient())
    #
    #     # modification du type de tous les points en Undefined au cazou
    #     for poi in obj_patient.case.PatientModel.PointsOfInterest:
    #         poi.Type = 'Undefined'
    #
    #     laser_rouge_HFS = "laser rouge HFS"
    #
    #     choix1 = "Toulouse"  # on positionne le laser rouge à deux endroits différents = relou
    #     choix2 = "Lille"  # On positionne les lasers rouges et vert sur le point abdomen mais le rouge est décalé sur la bille genoux en GD
    #
    #     choix = choix2
    #
    #     if choix == choix1:
    #         coords_laser_rouges_HFS = (0, obj_patient.zero_scan, obj_patient.thorax[2])
    #         # laser vert positionné 13 cm en dessous du laser rouge
    #         shift = -13
    #         coords_laser_vert = [0, coords_laser_rouges_HFS[1], obj_patient.thorax[2] + shift]
    #     elif choix == choix2:
    #         coords_laser_rouges_HFS = (obj_patient.genoux[0], obj_patient.zero_scan, obj_patient.abdomen[2])
    #         coords_laser_vert = (0, obj_patient.zero_scan, obj_patient.abdomen[2])
    #
    #     # création du laser rouge
    #     obj_patient.create_poi(laser_rouge_HFS, coords_laser_rouges_HFS, color='Red')
    #     # laser rouge devient localization point
    #     obj_patient.case.PatientModel.PointsOfInterest[laser_rouge_HFS].Type = "LocalizationPoint"
    #     # création du laser vert
    #     obj_patient.create_poi('laser vert HFS', coords_laser_vert, color='Green')
    #
    #     ########################################################
    #     # Création du faisceau
    #
    #     plan = obj_patient.case.TreatmentPlans[plan_name]
    #     x, y, z = coords_laser_vert
    #
    #     try:
    #         retval_0 = beam_set.CreatePhotonBeam(BeamQualityId="6", CyberKnifeCollimationType="Undefined",
    #                                              CyberKnifeNodeSetName=None, CyberKnifeRampVersion=None,
    #                                              CyberKnifeAllowIncreasedPitchCorrection=None,
    #                                              IsocenterData={'Position': {'x': x, 'y': y, 'z': z},
    #                                                             'NameOfIsocenterToRef': "", 'Name': 'iso_laser_vert',
    #                                                             'Color': "98, 184, 234"}, Name=bs_name, Description="",
    #                                              GantryAngle=0, CouchRotationAngle=0, CouchPitchAngle=0,
    #                                              CouchRollAngle=0,
    #                                              CollimatorAngle=0)
    #
    #         retval_0.SetBolus(BolusName="")
    #         beam_set.BeamMU = 0
    #     except:
    #         # le faisceau existe déjà.
    #         Exception("le faisceau existe déjà")
    #
    #     # modification des paramètres du faisceau
    #     pitch = 0.38
    #     gantry_period = 22
    #     try:
    #         plan.PlanOptimizations[0].OptimizationParameters.TreatmentSetupSettings[0].BeamSettings[
    #             0].TomoPropertiesPerBeam.EditTomoBasedBeamOptimizationSettings(JawMode="Dynamic",
    #                                                                            PitchTomoHelical=pitch,
    #                                                                            PitchTomoDirect=None,
    #                                                                            BackJawPosition=2.1,
    #                                                                            FrontJawPosition=-2.1,
    #                                                                            MaxDeliveryTime=None,
    #                                                                            MaxGantryPeriod=gantry_period,
    #                                                                            MaxDeliveryTimeFactor=None)
    #     except:
    #         print("No changes to save.")
    #
    #     #################################################################################################################
    #     #################################################################################################################
    #     ############################################   OPTIMISATION     #################################################
    #     #################################################################################################################
    #     #################################################################################################################
    #
    #     if ind == 0:  # Plan HFS
    #
    #         # Ici on associe à chaque PTV un certain niveau par rapport à la dose de prescription. Travail avec numpy pour
    #         # pouvoir utiliser la méthode np.where, associant la dose à la roi
    #         ptv_list = np.array(
    #             ['PTV_A', 'PTV_B', 'PTV_C', 'PTV_D1', 'PTV_D2', 'PTV_D3', 'PTV_D4', 'PTV_D5', 'PTV_D6', 'PTVpoumons'])
    #
    #         # Définition de la contrainte aux poumons. Si la dose prescrite est supérieure à 8Gy, on applique 7.5 Gy en
    #         # contrainte aux poumons (vu avec SP, ça marche mieux que de mettre 8)
    #         if total_dose > 800:
    #             dose_poumons = 750
    #         # Sinon on garde la même dose aux poumons qu'au reste.
    #         else:
    #             dose_poumons = total_dose
    #
    #         objectives_list = [total_dose, total_dose, total_dose, total_dose, 0.9 * total_dose, 0.63 * total_dose,
    #                            0.37 * total_dose, 0.16 * total_dose, 0.05 * total_dose, dose_poumons]
    #     elif ind == 1:  # Plan FFS
    #         ptv_list = np.array(
    #             ['PTV_E', 'PTV_D6', 'PTV_D5', 'PTV_D4', 'PTV_D3', 'PTV_D2', 'PTV_D1'])
    #
    #         objectives_list = [total_dose, total_dose, 0.9 * total_dose, 0.63 * total_dose, 0.37 * total_dose,
    #                            0.16 * total_dose, 0.05 * total_dose, ]
    #
    #     for roi in ptv_list:
    #         robustesse = False
    #         plan.PlanOptimizations[0].AddOptimizationFunction(FunctionType="UniformDose", RoiName=str(roi),
    #                                                           IsConstraint=False, RestrictAllBeamsIndividually=False,
    #                                                           RestrictToBeam=None, IsRobust=robustesse,
    #                                                           RestrictToBeamSet=None, UseRbeDose=False)
    #         plan.PlanOptimizations[0].AddOptimizationFunction(FunctionType="MinDose", RoiName=str(roi),
    #                                                           IsConstraint=False, RestrictAllBeamsIndividually=False,
    #                                                           RestrictToBeam=None, IsRobust=robustesse,
    #                                                           RestrictToBeamSet=None, UseRbeDose=False)
    #         plan.PlanOptimizations[0].AddOptimizationFunction(FunctionType="MaxDose", RoiName=str(roi),
    #                                                           IsConstraint=False, RestrictAllBeamsIndividually=False,
    #                                                           RestrictToBeam=None, IsRobust=robustesse,
    #                                                           RestrictToBeamSet=None, UseRbeDose=False)
    #
    #         # On attribue la dose totale à tous les objectifs
    #         for dose in plan.PlanOptimizations[0].Objective.ConstituentFunctions:
    #             # On regarde dans le xième objectif à quelle roi il est associé
    #             associated_roi = dose.OfDoseGridRoi.OfRoiGeometry.OfRoi.Name
    #             print(associated_roi)
    #             # On recherche dans les np.array créée ci -dessus à quel objectif cela correspond
    #             indice = np.where(ptv_list == associated_roi)[0][0]
    #             print(indice)
    #             # On attribue la dose correspondante
    #             dose.DoseFunctionParameters.DoseLevel = objectives_list[indice]
    #
    #         obj_patient.patient.Save()
    #
    #         # # Robustesse
    #         #
    #         # decalage = 0  # cm
    #         # plan.PlanOptimizations[0].OptimizationParameters.SaveRobustnessParameters(PositionUncertaintyAnterior=0,
    #         #                                                                           PositionUncertaintyPosterior=0,
    #         #                                                                           PositionUncertaintySuperior=0,
    #         #                                                                           PositionUncertaintyInferior=0,
    #         #                                                                           PositionUncertaintyLeft=decalage,
    #         #                                                                           PositionUncertaintyRight=decalage,
    #         #                                                                           DensityUncertainty=0,
    #         #                                                                           PositionUncertaintySetting="Universal",
    #         #                                                                           IndependentLeftRight=True,
    #         #                                                                           IndependentAnteriorPosterior=True,
    #         #                                                                           IndependentSuperiorInferior=True,
    #         #                                                                           ComputeExactScenarioDoses=False,
    #         #                                                                           NamesOfNonPlanningExaminations=[],
    #         #                                                                           PatientGeometryUncertaintyType="PerTreatmentCourse",
    #         #                                                                           PositionUncertaintyType="PerTreatmentCourse",
    #         #                                                                           TreatmentCourseScenariosFactor=1000)
