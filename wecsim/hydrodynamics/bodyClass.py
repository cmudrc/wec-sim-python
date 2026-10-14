# -*- coding: utf-8 -*-
"""
Created on Fri Apr  3 07:52:07 2020

@author: SungJun Won

This code is written based on WEC-sim.
wonsungjun0000@gmail.com

Note: All the struct() has been changed to dictionary

Note: trimesh is used to obtain bodyGeometry properties

"""
import h5py
import numpy as np
import warnings
import trimesh

from numpy.linalg import inv
from scipy import interpolate
from scipy.signal import fftconvolve


def _directional_interp(frequencies, directions, values, query_frequencies, query_directions):
    """Interpolate a direction-by-frequency hydrodynamic table linearly."""
    frequencies = np.asarray(frequencies).ravel()
    directions = np.asarray(directions).ravel()
    query_frequencies = np.atleast_1d(query_frequencies)
    query_directions = np.atleast_1d(query_directions)
    frequency_grid, direction_grid = np.meshgrid(query_frequencies, query_directions)
    points = np.column_stack((direction_grid.ravel(), frequency_grid.ravel()))
    interpolator = interpolate.RegularGridInterpolator(
        (directions, frequencies), np.asarray(values), bounds_error=False,
        fill_value=None,
    )
    return interpolator(points).reshape(direction_grid.shape)


class BodyClass:
    def hdf5FileProperties(self):#hdf5 file properties
        """
        Hydrodynamic data from BEM or user defined.
        
        """
        self.hydroData={'simulation_parameters':{'scaled':[],
                                                 'wave_dir':[],
                                                 'water_depth':[],
                                                 'w':[],
                                                 'T':[]},
                        'properties':{'name':[],
                                      'body_number':[],
                                      'cg':[],
                                      'cb':[],
                                      'disp_vol':[],
                                      'dof_start':[],
                                      'dof_end':[]},
                        'hydro_coeffs':{'linear_restoring_stiffness':[],
                                        'excitation':{'re':[],
                                                      'im':[],
                                                      'impulse_response_fun':{'f':[],
                                                                              't':[]}},
                                        'added_mass':{'all':[],
                                                      'inf_freq':[]},
                                        'radiation_damping':{'all':[],
                                                             'impulse_response_fun':{'K':[],
                                                                                     't':[]},
                                                             'state_space':{'it':[],
                                                                            'A':{'all':[]},
                                                                            'B':{'all':[]},
                                                                            'C':{'all':[]},
                                                                            'D':{'all':[]}}},
                                        'mean_drift':[]},
                        'gbm':{'mass':[],
                               'stiffness':[],
                               'damping':[]}}                                                  

    def inputFileProperties(self):#input file properties
        self.name              = []                                            # Body name. For WEC bodies this is given in the h5 file.
        self.mass              = []                                            # Mass in kg or specify 'equilibrium' to have mass= dis vol * density
        self.momOfInertia      = []                                            # Moment of inertia [Ixx, Iyy, Izz] in kg*m^2
        self.cg                = []                                            # Center of gravity [x, y, z] in meters. For WEC bodies this is given in the h5 file.
        self.cb                = []                                            # Center of buoyancy [x, y, z] in meters. For WEC bodies this is given in the h5 file.
        self.dispVol           = []                                            # Displaced volume at equilibrium position in meters cubed. For WEC bodies this is given in the h5 file.
        self.dof               = []                                            # Number of DOFs. For WEC bodies this is given in the h5 file. IF not, default is 6
        self.dof_gbm           = []                                            # Number of DOFs for GBM.
        self.dof_start         = []                                            # Index of DOF starts. For WEC bodies this is given in the h5 file. IF not, default is (bodyNumber-1)*6
        self.dof_end           = []                                            # Index of DOF ends. For WEC bodies this is given in the h5 file. IF not, default is (bodyNumber-1)*6+5
        self.geometryFile      = 'NONE'                                        # Location of geomtry stl files
        self.viscDrag          = {                                             # Structure defining the viscous (quadratic) drag
                                  'Drag':                 [0,0,0,0,0,0],       # Viscous (quadratic) drag, matrix 6x6
                                  'cd':                   [0,0,0,0,0,0],       # Viscous (quadratic) drag cd coefficient, vector length 6
                                  'characteristicArea':   [0,0,0,0,0,0]}       # Characteristic area for viscous drag, vector length 6
        self.initDisp          = {                                             # Structure defining the initial displacement
                                  'initLinDisp':          [0,0,0],             # Initial displacement of center fo gravity - used for decay tests (format: [displacment in m], default = [0, 0, 0])
                                  'initAngularDispAxis':  [0,1,0],             # Initial displacement of cog - axis of rotation - used for decay tests (format: [x y z], default = [1, 0, 0])
                                  'initAngularDispAngle': 0}                   # Initial displacement of cog - Angle of rotation - used for decay tests (format: [radians], default = 0)
        self.hydroStiffness    = [0,0,0,0,0,0]                                 # Hydrostatic stiffness matrix overrides BEMIO definition, matrix 6x6
        self.linearDamping     = [0,0,0,0,0,0]                                 # Linear damping coefficient, matrix size of 6x6
        self.userDefinedExcIRF = []                                            # Excitation IRF from BEMIO used for User-Defined Time-Series
        self.viz               = {                                             # Structure defining visualization properties
                                  'color': [1,1,0],                            # Visualization color for either SimMechanics Explorer or Paraview.
                                  'opacity': 1}                                # Visualization opacity for either SimMechanics Explorer or Paraview.
        self.morisonElement   = {                                              # Structure defining the Morrison Elements
                                 'cd':                 [0,0,0],                # Viscous (quadratic) drag cd, vector length 3
                                 'ca':                 [0,0,0],                # Added mass coefficent for Morrison Element (format [Ca_x, Ca_y, Ca_z], default = [0 0 0])
                                 'characteristicArea': [0,0,0],                # Characteristic area for Morrison Elements calculations (format [Area_x, Area_y, Area_z], default = [0 0 0])
                                 'VME':                 0     ,                # Characteristic volume for Morrison Element (default = 0)
                                 'rgME':               [0,0,0]}                # Vector from center of gravity to point of application for Morrison Element (format [X, Y, Z], default = [0 0 0]).
        self.nhBody            = 0                                             # Flag for non-hydro body.
        self.flexHydroBody     = 0                                             # Flag for flexible body. 
        self.meanDriftForce    = 0                                             # Flag for mean drift force. 0: No 1: from control surface 2: from momentum conservation.
        
    def bodyGeometryFileProperties(self):                                      # body geometry stl file properties
        self.bodyGeometry = {                                                  # Structure defining body's mesh
                             'numFace': [],                                    # Number of faces
                             'numVertex': [],                                  # Number of vertices
                             'vertex': [],                                     # List of vertices
                             'face': [],                                       # List of faces
                             'norm': [],                                       # List of normal vectors
                             'area': [],                                       # List of cell areas
                             'center': []}                                     # List of cell centers
        self.meshFile    = []                                                  # body geometry stl file from Trimesh api
    
    def internalProperties(self,filename):                                     # internal properties 
        self.hydroForce        = {'linearHydroRestCoef':[],                    # Hydrodynamic forces and coefficients used during simulation.
                                  'visDrag':[],
                                  'linearDamping':[],
                                  'userDefinedFe':[],
                                  'fExt':{'re':[],
                                          'im':[],
                                          'md':[]},
                                  'fAddedMass':[],
                                  'fDamping':[],
                                  'irkb':[],
                                  'ssRadf':{'A':[],
                                            'B':[],
                                            'C':[],
                                            'D':[]},
                                  'gbm':{'mass_ff':[],
                                         'mass_ff_inv':[],
                                         'stiffness':[],
                                         'damping':[],
                                         'state_space':{}},
                                  'storage':{'mass':[],
                                             'momOfInertia':[],
                                             'fAddedMass':[],
                                             'output_forceAddedMass':[],
                                             'output_forceTotal':[]}}           
        self.tmp               ={'fadm':[],                                    # temporary file
                                 'adjmass':[],
                                 'mass':[],
                                 'momOfInertia':[],
                                 'hydroForce_fAddedMass':[]}                    
        self.h5File            = filename                                      # hdf5 file containing the hydrodynamic data
        self.hydroDataBodyNum  = []                                            # Body number within the hdf5 file.
        self.massCalcMethod    = []                                            # Method used to obtain mass: 'user', 'fixed', 'equilibrium'
        self.bodyNumber        = []                                            # bodyNumber in WEC-Sim as defined in the input file. Can be different from the BEM body number.
        self.bodyTotal         = 0                                             # Total number of WEC-Sim bodies (body block iterations)
        self.lenJ              = []                                            # Matrices length. 6 for no body-to-body interactions. 6*numBodies if body-to-body interactions.

    def __init__(self,filename):
        """
        Initialize Body Class
        Takes string parameter called filename
        filename: name of h5 file

        """
        self.bodyGeometryFileProperties()
        self.internalProperties(filename)
        self.inputFileProperties()
        self.hdf5FileProperties()
        self.meanDriftForce = 0

    def readH5file(self):
        """
        Read and record properties of h5 file in self.hydroData

        """
        f = h5py.File(self.h5File, 'r')
        name = '/body' + str(self.bodyNumber)
        self.cg = np.transpose(np.array(f.get(name + '/properties/cg')))
        self.cb = np.transpose(np.array(f.get(name + '/properties/cb')))
        self.dispVol = np.array(f.get(name + '/properties/disp_vol'))
        self.name = f[name + '/properties/name'].asstr()[()]
        self.hydroData['simulation_parameters']['scaled'] = np.array(f.get('/simulation_parameters/scaled'))
        self.hydroData['simulation_parameters']['wave_dir'] = np.transpose(np.array(f.get('/simulation_parameters/wave_dir')))
        if np.array(f.get('/simulation_parameters/water_depth')).dtype.kind in 'fi':
            self.hydroData['simulation_parameters']['water_depth'] = np.array(f.get('/simulation_parameters/water_depth'))
        else:            
            self.hydroData['simulation_parameters']['water_depth'] = f['/simulation_parameters/water_depth'].asstr()[()]
        self.hydroData['simulation_parameters']['w'] = np.transpose(np.array(f.get('/simulation_parameters/w')))
        self.hydroData['simulation_parameters']['T'] = np.transpose(np.array(f.get('/simulation_parameters/T')))
        self.hydroData['properties']['name'] = self.name
        self.hydroData['properties']['body_number'] = np.array(f.get(name + '/properties/body_number'))
        self.hydroData['properties']['cg'] = np.transpose(np.array(f.get(name + '/properties/cg')))
        self.hydroData['properties']['cb'] = np.transpose(np.array(f.get(name + '/properties/cb')))
        self.hydroData['properties']['disp_vol'] = np.array(f.get(name + '/properties/disp_vol'))
        if f.get(name + '/properties/dof') is not None:
            self.hydroData['properties']['dof'] = np.array(f.get(name +'/properties/dof')) 
        else:
            self.hydroData['properties']['dof'] = np.array(6)
        if f.get(name + '/properties/dof_start') is not None:
            self.hydroData['properties']['dof_start'] = np.array(f.get(name + '/properties/dof_start'))
        else:
            self.hydroData['properties']['dof_start'] = np.array((self.bodyNumber-1)*6+1)
        if f.get(name + '/properties/dof_end') is not None:
            self.hydroData['properties']['dof_end'] = np.array(f.get(name + '/properties/dof_end'))
        else:
            self.hydroData['properties']['dof_end'] = np.array((self.bodyNumber-1)*6+6)
        self.dof       = np.asarray(self.hydroData['properties']['dof']).reshape(-1)
        self.dof_start = np.asarray(self.hydroData['properties']['dof_start']).reshape(-1)
        self.dof_end   = np.asarray(self.hydroData['properties']['dof_end']).reshape(-1)
        self.dof_gbm   = self.dof-6
        restoring = np.array(f.get(name + '/hydro_coeffs/linear_restoring_stiffness'))
        # The 2D HDF5 matrix is already in the orientation used by MATLAB's
        # h5load. Transposing it moves asymmetric OWC yaw couplings.
        self.hydroData['hydro_coeffs']['linear_restoring_stiffness'] = (
            restoring[:, :, 0] if restoring.ndim == 3 and restoring.shape[2] == 1
            else restoring
        )
        self.hydroData['hydro_coeffs']['excitation']['re'] = np.array(f.get(name +  '/hydro_coeffs/excitation/re'))
        self.hydroData['hydro_coeffs']['excitation']['im'] = np.array(f.get(name + '/hydro_coeffs/excitation/im'))
        if f.get(name + '/hydro_coeffs/excitation/impulse_response_fun/f') is not None:
            self.hydroData['hydro_coeffs']['excitation']['impulse_response_fun']['f'] = np.array(f.get(name + '/hydro_coeffs/excitation/impulse_response_fun/f'))
        if f.get(name + '/hydro_coeffs/excitation/impulse_response_fun/t') is not None:
            self.hydroData['hydro_coeffs']['excitation']['impulse_response_fun']['t'] = np.array(f.get(name + '/hydro_coeffs/excitation/impulse_response_fun/t'))
        self.hydroData['hydro_coeffs']['added_mass']['all'] = np.array(f.get(name + '/hydro_coeffs/added_mass/all'))
        self.hydroData['hydro_coeffs']['added_mass']['inf_freq'] = np.array(f.get(name + '/hydro_coeffs/added_mass/inf_freq'))
        self.hydroData['hydro_coeffs']['radiation_damping']['all'] = np.array(f.get(name + '/hydro_coeffs/radiation_damping/all'))
        if f.get(name + '/hydro_coeffs/radiation_damping/impulse_response_fun/K') is not None:
            self.hydroData['hydro_coeffs']['radiation_damping']['impulse_response_fun']['K'] = np.array(f.get(name + '/hydro_coeffs/radiation_damping/impulse_response_fun/K'))
        if f.get(name + '/hydro_coeffs/radiation_damping/impulse_response_fun/t') is not None:
            self.hydroData['hydro_coeffs']['radiation_damping']['impulse_response_fun']['t'] = np.transpose(np.array(f.get(name +'/hydro_coeffs/radiation_damping/impulse_response_fun/t')))
        if f.get(name + '/hydro_coeffs/radiation_damping/state_space/it') is not None:
            self.hydroData['hydro_coeffs']['radiation_damping']['state_space']['it'] = np.array(f.get(name + '/hydro_coeffs/radiation_damping/state_space/it'))
        if f.get(name + '/hydro_coeffs/radiation_damping/state_space/A/all') is not None:
            self.hydroData['hydro_coeffs']['radiation_damping']['state_space']['A']['all'] = np.array(f.get(name + '/hydro_coeffs/radiation_damping/state_space/A/all'))
        if f.get(name + '/hydro_coeffs/radiation_damping/state_space/B/all') is not None:
            self.hydroData['hydro_coeffs']['radiation_damping']['state_space']['B']['all'] = np.array(f.get(name + '/hydro_coeffs/radiation_damping/state_space/B/all'))
        if f.get(name + '/hydro_coeffs/radiation_damping/state_space/C/all') is not None:
            self.hydroData['hydro_coeffs']['radiation_damping']['state_space']['C']['all'] = np.array(f.get(name + '/hydro_coeffs/radiation_damping/state_space/C/all'))
        if f.get(name + '/hydro_coeffs/radiation_damping/state_space/D/all') is not None:
            self.hydroData['hydro_coeffs']['radiation_damping']['state_space']['D']['all'] = np.array(f.get(name +'/hydro_coeffs/radiation_damping/state_space/D/all'))
        gbm_start = int(np.asarray(self.dof_start).item()) + 5
        gbm_end = int(np.asarray(self.dof_end).item())
        for property_name in ('mass', 'stiffness', 'damping'):
            dataset = f.get(name + '/properties/' + property_name)
            if dataset is not None:
                # HDF5 reverses MATLAB's matrix dimension order.
                matrix = np.asarray(dataset)
                matrix = (matrix[:, :, 0] if matrix.ndim == 3 and matrix.shape[2] == 1
                          else matrix.T)
                self.hydroData['gbm'][property_name] = matrix[
                    gbm_start:gbm_end, gbm_start:gbm_end
                ]
        if self.meanDriftForce == 0:
            self.hydroData['hydro_coeffs']['mean_drift'] = 0.*self.hydroData['hydro_coeffs']['excitation']['re']
        elif self.meanDriftForce == 1:
            self.hydroData['hydro_coeffs']['mean_drift'] = np.array(f.get(name + '/hydro_coeffs/mean_drift/control_surface/val'))
        elif self.meanDriftForce == 2:
            self.hydroData['hydro_coeffs']['mean_drift'] = np.array(f.get(name + '/hydro_coeffs/mean_drift/momentum_conservation/val'))
        else:
            warnings.warn("Wrong flag for mean drift force.",DeprecationWarning)
        f.close()
        
    def loadHydroData(self, hydroData):
        """
        Loads user defiend hydroData structure as alternative
        to reading the h5 file. Used in wecSimMCR

        """
        self.hydroData = hydroData
        self.cg        = hydroData['properties']['cg']
        self.cb        = hydroData['properties']['cb']
        self.dispVol   = hydroData['properties']['disp_vol']
        self.name      = hydroData['properties']['name']
        self.dof       = np.asarray(self.hydroData['properties']['dof']).reshape(-1)
        self.dof_start = np.asarray(self.hydroData['properties']['dof_start']).reshape(-1)
        self.dof_end   = np.asarray(self.hydroData['properties']['dof_end']).reshape(-1)
        self.dof_gbm   = self.dof-6
    
    def hydroForcePre(self,w,waveDir,CIkt,CTTime,numFreq,dt,rho,g,waveType,waveAmpTime,iBod,numBod,ssCalc,nlHydro,B2B):
        """
        This is the most important method of this class
        HydroForce Pre-processing calculations
        1. Set the linear hydrodynamic restoring coefficient, viscous
            drag, and linear damping matrices
        2. Set the wave excitation force    

        """
        self.setMassMatrix(rho,nlHydro)
        if self.dof_gbm > 0:
            dof = int(np.asarray(self.dof).item())

            def extend_matrix(value):
                matrix = np.asarray(value, dtype=float)
                if matrix.ndim == 1:
                    matrix = np.diag(matrix)
                if matrix.ndim != 2 or max(matrix.shape) > dof:
                    raise ValueError("GBM damping matrix exceeds body DOFs")
                extended = np.zeros((dof, dof))
                extended[:matrix.shape[0], :matrix.shape[1]] = matrix
                return extended

            self.linearDamping = extend_matrix(self.linearDamping)
            self.viscDrag['Drag'] = extend_matrix(self.viscDrag['Drag'])
            for key in ('cd', 'characteristicArea'):
                values = np.asarray(self.viscDrag[key], dtype=float).ravel()
                if len(values) > dof:
                    raise ValueError("GBM drag vector exceeds body DOFs")
                self.viscDrag[key] = np.pad(values, (0, dof - len(values)))

        if np.any(self.hydroStiffness):  # check for a user-defined stiffness
            self.hydroForce['linearHydroRestCoef'] = self.hydroStiffness
        else:
            k = self.hydroData['hydro_coeffs']['linear_restoring_stiffness']#(:,self.dof_start:self.dof_end)
            self.hydroForce['linearHydroRestCoef'] = k*rho*g

        if np.any(self.viscDrag['Drag']):  # check for a user-defined drag matrix
            self.hydroForce['visDrag'] = self.viscDrag['Drag']
        else:
            self.hydroForce['visDrag'] = np.diag(
                0.5 * rho * np.asarray(self.viscDrag['cd'])
                * np.asarray(self.viscDrag['characteristicArea'])
            )

        self.hydroForce['linearDamping'] = self.linearDamping
        self.hydroForce['userDefinedFe'] = np.zeros(
            (len(waveAmpTime[1]), int(np.asarray(self.dof).item()))
        )  # initialize userDefinedFe for non-imported wave cases
        if waveType == 'noWave':
            self.noExcitation()
            self.constAddedMassAndDamping(w,CIkt,rho,B2B)
        elif waveType == 'noWaveCIC':
            self.noExcitation()
            self.irfInfAddedMassAndDamping(CIkt,CTTime,ssCalc,rho,B2B)
        elif waveType == 'regular':
            self.regExcitation(w,waveDir,rho,g)
            self.constAddedMassAndDamping(w,CIkt,rho,B2B)
        elif waveType == 'regularCIC':
            self.regExcitation(w,waveDir,rho,g)
            self.irfInfAddedMassAndDamping(CIkt,CTTime,ssCalc,rho,B2B)
        elif waveType == 'irregular' or waveType == 'spectrumImport':
            self.irrExcitation(w,numFreq,waveDir,rho,g)
            self.irfInfAddedMassAndDamping(CIkt,CTTime,ssCalc,rho,B2B)
        elif waveType in ('etaImport', 'elevationImport'):
            self.userDefinedExcitation(waveAmpTime,dt,waveDir,rho,g)
            self.irfInfAddedMassAndDamping(CIkt,CTTime,ssCalc,rho,B2B)

        gbmDOF = int(np.asarray(self.dof_gbm).item())
        if gbmDOF > 0:
            flexible = slice(6, 6 + gbmDOF)
            gbm = self.hydroForce['gbm']
            gbm['stiffness'] = np.asarray(self.hydroData['gbm']['stiffness'])
            gbm['damping'] = np.asarray(self.hydroData['gbm']['damping'])
            gbm['mass_ff'] = (self.hydroForce['fAddedMass'][flexible, flexible]
                              + np.asarray(self.hydroData['gbm']['mass']))
            # MATLAB's flexible state-space block carries this mass, so its
            # applied added-mass matrix no longer contains the flex block.
            self.hydroForce['fAddedMass'][flexible, flexible] = 0
            gbm['mass_ff_inv'] = inv(gbm['mass_ff'])
            gbm['state_space']['A'] = np.block([
                [np.zeros((gbmDOF, gbmDOF)), np.eye(gbmDOF)],
                [-gbm['mass_ff_inv'] @ gbm['stiffness'],
                 -gbm['mass_ff_inv'] @ gbm['damping']],
            ])
            gbm['state_space']['B'] = np.eye(2 * gbmDOF)
            gbm['state_space']['C'] = np.eye(2 * gbmDOF)
            gbm['state_space']['D'] = np.zeros((2 * gbmDOF, 2 * gbmDOF))
            self.flexHydroBody = 1
            self.nhBody=0

    def listInfo(self):
        """
        List body info{:f}\n'.format(self.T)
        """
        print('\n\t***** Body Number ', self.hydroData['properties']['body_number'][0][0] ,', Name: ' , self.hydroData['properties']['name'] , ' *****\n')
        print('\tBody CG                          (m) = [',self.hydroData['properties']['cg'][0],']\n')
        print('\tBody Mass                       (kg) = ',self.mass[0][0],' \n')
        print('\tBody Diagonal MOI              (kgm2)= [',self.momOfInertia,']\n')

    def bodyGeo(self,fname):
        """
        Reads mesh file and calculates areas and centroids
        Use trimesh to get geometry values besides using the method from WEC-Sim
        Takes in string parameter called fname. fname is name of stl file

        """
        your_mesh = trimesh.load_mesh(fname)
        self.meshFile = your_mesh
        v = your_mesh.triangles.reshape(int(np.size(your_mesh.triangles)/3),3)
        [vertex,faces] = np.unique(v,return_inverse=True,axis=0)
        face = faces.reshape(int(np.size(faces)/3),3)+1
        self.bodyGeometry['vertex'] = vertex
        self.bodyGeometry['numVertex'] = np.size(vertex,0)
        self.bodyGeometry['face'] = face
        self.bodyGeometry['numFace'] = np.size(face,0)
        self.bodyGeometry['norm'] = your_mesh.face_normals
        self.bodyGeometry['center'] = your_mesh.triangles_center
        self.bodyGeometry['area'] = np.transpose([your_mesh.area_faces])        
        
    def plotStl(self):
        """
        Plots the body's mesh

        """
        # Plots the body's mesh
        your_mesh = self.meshFile
        your_mesh.show()
                    
    def checkinputs(self):
        """
        Checks the user inputs

        """
        # hydro data file
        if self.h5File == None and self.nhBody == 0:
            warnings.warn("The hdf5 file does not exist")
        # geometry file
        if self.geometryFile == None:
            warnings.warn("Could not locate and open geometry file")
        

    def noExcitation(self):
        """
        Set exciation force for no excitation case
        Gets used in hydroForcePre for noWave condition
        """
        nDOF = int(self.dof[0])
        self.hydroForce['fExt']['re'] = np.zeros(nDOF)
        self.hydroForce['fExt']['im'] = np.zeros(nDOF)
    
    
    def regExcitation(self,w,waveDir,rho,g):
        """
        Regular wave excitation force
        Used by hydroForcePre       

        """
        nDOF = int(self.dof[0])
        re = self.hydroData['hydro_coeffs']['excitation']['re']*rho*g
        im = self.hydroData['hydro_coeffs']['excitation']['im']*rho*g
        md = self.hydroData['hydro_coeffs']['mean_drift']*rho*g
        self.hydroForce['fExt']['re'] = np.zeros(nDOF)
        self.hydroForce['fExt']['im'] = np.zeros(nDOF)
        self.hydroForce['fExt']['md'] = np.zeros(nDOF)
        for ii in range(nDOF):
            if np.size(self.hydroData['simulation_parameters']['wave_dir']) > 1:
                x = self.hydroData['simulation_parameters']['w'][0]
                y = self.hydroData['simulation_parameters']['wave_dir'][0]
                self.hydroForce['fExt']['re'][ii] = _directional_interp(x, y, np.squeeze(re[ii]), w, waveDir).item()
                self.hydroForce['fExt']['im'][ii] = _directional_interp(x, y, np.squeeze(im[ii]), w, waveDir).item()
                self.hydroForce['fExt']['md'][ii] = _directional_interp(x, y, np.squeeze(md[ii]), w, waveDir).item()
            elif self.hydroData['simulation_parameters']['wave_dir'] == waveDir:
                x = self.hydroData['simulation_parameters']['w'][0]
                frequency = np.asarray(w).item()
                s1 = interpolate.CubicSpline(x, np.squeeze(re[ii][0]))# interpolate using CubicSline to get interpolation of spline 3d space
                s2 = interpolate.CubicSpline(x, np.squeeze(im[ii][0]))
                s3 = interpolate.CubicSpline(x, np.squeeze(md[ii][0]))
                self.hydroForce['fExt']['re'][ii] = s1(frequency)
                self.hydroForce['fExt']['im'][ii] = s2(frequency)
                self.hydroForce['fExt']['md'][ii] = s3(frequency)

    def irrExcitation(self,wv,numFreq,waveDir,rho,g):
        """
        Irregular wave excitation force
        Used by hydroForcePre

        """
        nDOF = int(self.dof[0])
        re = self.hydroData['hydro_coeffs']['excitation']['re']*rho*g
        im = self.hydroData['hydro_coeffs']['excitation']['im']*rho*g
        md = self.hydroData['hydro_coeffs']['mean_drift']*rho*g
        self.hydroForce['fExt']['re'] = np.zeros((np.size(waveDir),numFreq,nDOF))
        self.hydroForce['fExt']['im'] = np.zeros((np.size(waveDir),numFreq,nDOF))
        self.hydroForce['fExt']['md'] = np.zeros((np.size(waveDir),numFreq,nDOF))
        for ii in range(nDOF):
            if np.size(self.hydroData['simulation_parameters']['wave_dir']) > 1:
                x = self.hydroData['simulation_parameters']['w'][0]
                y = self.hydroData['simulation_parameters']['wave_dir'][0]
                self.hydroForce['fExt']['re'][:,:,ii] = _directional_interp(x, y, np.squeeze(re[ii]), wv, waveDir)
                self.hydroForce['fExt']['im'][:,:,ii] = _directional_interp(x, y, np.squeeze(im[ii]), wv, waveDir)
                self.hydroForce['fExt']['md'][:,:,ii] = _directional_interp(x, y, np.squeeze(md[ii]), wv, waveDir)
            elif self.hydroData['simulation_parameters']['wave_dir'] == waveDir:
                x = self.hydroData['simulation_parameters']['w'][0]
                s1 = interpolate.CubicSpline(x, np.squeeze(re[ii][0]))# interpolate using CubicSline to get interpolation of spline 3d space
                s2 = interpolate.CubicSpline(x, np.squeeze(im[ii][0]))
                s3 = interpolate.CubicSpline(x, np.squeeze(md[ii][0]))
                self.hydroForce['fExt']['re'][:,:,ii] = s1(wv)
                self.hydroForce['fExt']['im'][:,:,ii] = s2(wv)
                self.hydroForce['fExt']['md'][:,:,ii] = s3(wv)
    
    def userDefinedExcitation(self,waveAmpTime,dt,waveDir,rho,g):
        """
        Calculated User-Defined wave excitation force with non-causal convolution
        Used by hydroForcePre
        
        """
        elevation = np.asarray(waveAmpTime, dtype=float)
        if (elevation.ndim != 2 or elevation.shape[0] != 2
                or not np.isfinite(elevation).all() or dt <= 0):
            raise ValueError("imported wave elevation must be finite 2-by-time data")
        nDOF = int(np.asarray(self.dof).item())
        irf = self.hydroData['hydro_coeffs']['excitation']['impulse_response_fun']
        kernel = np.asarray(irf['f'], dtype=float) * rho * g
        source_lag = np.asarray(irf['t'], dtype=float).ravel()
        count = round((source_lag[-1] - source_lag[0]) / dt)
        lag = source_lag[0] + np.arange(count + 1) * dt
        directions = np.asarray(
            self.hydroData['simulation_parameters']['wave_dir'], dtype=float,
        ).ravel()
        direction = float(np.asarray(waveDir).item())
        if len(directions) == 1 and not np.isclose(directions[0], direction):
            raise ValueError("imported wave direction differs from the hydro database")
        for ii in range(nDOF):
            if len(directions) > 1:
                sampled_kernel = _directional_interp(
                    source_lag, directions, kernel[ii], lag, direction,
                ).ravel()
            else:
                sampled_kernel = np.interp(lag, source_lag, kernel[ii, 0])
            self.userDefinedExcIRF = sampled_kernel
            self.hydroForce['userDefinedFe'][:, ii] = (
                fftconvolve(elevation[1], sampled_kernel, mode='same') * dt
            )
        
        self.hydroForce['fExt']['re'] = np.zeros(nDOF)
        self.hydroForce['fExt']['im'] = np.zeros(nDOF)
        self.hydroForce['fExt']['md'] = np.zeros(nDOF)
    
    
    def constAddedMassAndDamping(self,w,CIkt,rho,B2B):
        """
        Set added mass and damping for a specific frequency
        Used by hydroForcePre
        
        """
        frequency = np.asarray(w).item()
        am = self.hydroData['hydro_coeffs']['added_mass']['all']*rho
        rd = self.hydroData['hydro_coeffs']['radiation_damping']['all']*rho
        for i in range(len(self.hydroData['simulation_parameters']['w'][0])):
            rd[:,:,i] *= self.hydroData['simulation_parameters']['w'][0][i]
        # Change matrix size: B2B [6x6n], noB2B [6x6]
        if B2B == 1:
            lenJ = 6*int(np.asarray(self.bodyTotal).item())
            self.hydroForce['fAddedMass'] = np.zeros((6,lenJ))
            self.hydroForce['fDamping'] = np.zeros((6,lenJ))
            self.hydroForce['totDOF']  = np.zeros((6,lenJ))
            for ii in range(6):
                for jj in range(lenJ):
                    s1 = interpolate.CubicSpline(self.hydroData['simulation_parameters']['w'][0], np.squeeze(am[ii,jj,:]))
                    self.hydroForce['fAddedMass'][ii,jj] = s1(frequency)
                    s2 = interpolate.CubicSpline(self.hydroData['simulation_parameters']['w'][0], np.squeeze(rd[ii,jj,:]))
                    self.hydroForce['fDamping'][ii,jj] = s2(frequency)
        else: #B2B =2
            nDOF = int(self.dof[0])
            self.hydroForce['fAddedMass'] = np.zeros((nDOF,nDOF))
            self.hydroForce['fDamping'] = np.zeros((nDOF,nDOF))
            self.hydroForce['totDOF']  = np.zeros((nDOF,nDOF))
            for ii in range(nDOF):
                for jj in range(nDOF):
                    jjj = int(self.dof_start[0])-1+jj
                    s1 = interpolate.CubicSpline(self.hydroData['simulation_parameters']['w'][0], np.squeeze(am[ii,jjj,:]))
                    self.hydroForce['fAddedMass'][ii,jj] = s1(frequency)
                    s2 = interpolate.CubicSpline(self.hydroData['simulation_parameters']['w'][0], np.squeeze(rd[ii,jjj,:]))
                    self.hydroForce['fDamping'][ii,jj] = s2(frequency)
    
    
    def irfInfAddedMassAndDamping(self,CIkt,CTTime,ssCalc,rho,B2B):
        """
        Set radiation force properties using impulse response function\
        Added mass at infinite frequency
        Convolution integral raditation dampingiBod
        State space formulation
        Used by hydroForcePre
        
        """
        nDOF = int(self.dof[0])
        if B2B == 1:
            LDOF = int(np.asarray(self.bodyTotal).item())*6
        else:
            LDOF = int(self.dof[0])
        
        # Convolution integral formulation
        if B2B == 1:
            self.hydroForce['fAddedMass'] = self.hydroData['hydro_coeffs']['added_mass']['inf_freq']*rho
        else:
            self.hydroForce['fAddedMass'] = self.hydroData['hydro_coeffs']['added_mass']['inf_freq'][:,int(self.dof_start[0])-1:int(self.dof_end[0])]*rho
        
        # Radition IRF
        self.hydroForce['fDamping'] = np.zeros((nDOF,LDOF))
        irfk = self.hydroData['hydro_coeffs']['radiation_damping']['impulse_response_fun']['K']*rho
        irft = self.hydroData['hydro_coeffs']['radiation_damping']['impulse_response_fun']['t'][0]
        self.hydroForce['irkb'] = np.zeros((len(CTTime),nDOF,LDOF))
        if B2B == 1:
            for ii in range(nDOF):
                for jj in range(LDOF):
                    s1 = interpolate.CubicSpline(irft, np.squeeze(irfk[ii,jj,:]))
                    self.hydroForce['irkb'][:,ii,jj] = s1(CTTime)
        else:
            for ii in range(nDOF):
                for jj in range(LDOF):
                    jjj = int(self.dof_start[0])-1+jj
                    s1 = interpolate.CubicSpline(irft, np.squeeze(irfk[ii,jjj,:]))
                    self.hydroForce['irkb'][:,ii,jj] = s1(CTTime)
                
        # State Space Formulation
        if ssCalc == 1:
            if B2B == 1:
                for ii in range(nDOF):
                    for jj in range(LDOF):
                        arraySize = int(self.hydroData['hydro_coeffs']['radiation_damping']['state_space']['it'][ii,jj])
                        if ii == 0 and jj == 0: # Begin construction of combined state, input, and output matrices
                            Af = np.zeros((arraySize,arraySize))
                            Bf = np.zeros((arraySize,LDOF))
                            Cf = np.zeros((nDOF,arraySize))
                            Af[:arraySize,:arraySize] = self.hydroData['hydro_coeffs']['radiation_damping']['state_space']['A']['all'][ii,jj,:arraySize,:arraySize]
                            Bf[:arraySize,jj]         = self.hydroData['hydro_coeffs']['radiation_damping']['state_space']['B']['all'][ii,jj,:arraySize,0]
                            Cf[ii,:arraySize]         = self.hydroData['hydro_coeffs']['radiation_damping']['state_space']['C']['all'][ii,jj,0,:arraySize]
                        else:
                            Ay = np.size(Af,0)
                            Ax = np.size(Af,1)
                            Af = np.pad(Af, ((0,arraySize),(0,arraySize)), mode='constant', constant_values=0)
                            Af[Ay:Ay+arraySize,Ax:Ax+arraySize] = self.hydroData['hydro_coeffs']['radiation_damping']['state_space']['A']['all'][ii,jj,:arraySize,:arraySize]
                            By = np.size(Bf,0)
                            Bf = np.pad(Bf, ((0,arraySize),(0,0)), mode='constant', constant_values=0)
                            Bf[By:By+arraySize,jj] = self.hydroData['hydro_coeffs']['radiation_damping']['state_space']['B']['all'][ii,jj,:arraySize,0]
                            Cx = np.size(Cf,1)
                            Cf = np.pad(Cf, ((0,0),(0,arraySize)), mode='constant', constant_values=0)
                            Cf[ii,Cx:Cx+arraySize]  = self.hydroData['hydro_coeffs']['radiation_damping']['state_space']['C']['all'][ii,jj,0,:arraySize]
                self.hydroForce['ssRadf']['D'] = np.zeros((nDOF,LDOF))
            else:
                for ii in range(nDOF):
                    for jj in range(int(self.dof_start[0])-1,int(self.dof_end[0])):
                        jInd = jj-int(self.dof_start[0])+1
                        arraySize = int(self.hydroData['hydro_coeffs']['radiation_damping']['state_space']['it'][ii,jj])
                        if ii == 0 and jInd == 0: # Begin construction of combined state, input, and output matrices
                            Af = np.zeros((arraySize,arraySize))
                            Bf = np.zeros((arraySize,LDOF))
                            Cf = np.zeros((nDOF,arraySize))
                            Af[:arraySize,:arraySize] = self.hydroData['hydro_coeffs']['radiation_damping']['state_space']['A']['all'][ii,jj,:arraySize,:arraySize]
                            Bf[:arraySize,jInd]       = self.hydroData['hydro_coeffs']['radiation_damping']['state_space']['B']['all'][ii,jj,:arraySize,0]
                            Cf[ii,:arraySize]         = self.hydroData['hydro_coeffs']['radiation_damping']['state_space']['C']['all'][ii,jj,0,:arraySize]
                        else:
                            Ay = np.size(Af,0)
                            Ax = np.size(Af,1)
                            Af = np.pad(Af, ((0,arraySize),(0,arraySize)), mode='constant', constant_values=0)
                            Af[Ay:Ay+arraySize,Ax:Ax+arraySize] = self.hydroData['hydro_coeffs']['radiation_damping']['state_space']['A']['all'][ii,jj,:arraySize,:arraySize]
                            By = np.size(Bf,0)
                            Bf = np.pad(Bf, ((0,arraySize),(0,0)), mode='constant', constant_values=0)
                            Bf[By:By+arraySize,jInd] = self.hydroData['hydro_coeffs']['radiation_damping']['state_space']['B']['all'][ii,jj,:arraySize,0]
                            Cx = np.size(Cf,1)
                            Cf = np.pad(Cf, ((0,0),(0,arraySize)), mode='constant', constant_values=0)
                            Cf[ii,Cx:Cx+arraySize]  = self.hydroData['hydro_coeffs']['radiation_damping']['state_space']['C']['all'][ii,jj,0,:arraySize]
                self.hydroForce['ssRadf']['D'] = np.zeros((nDOF,nDOF))
            self.hydroForce['ssRadf']['A'] = Af
            self.hydroForce['ssRadf']['B'] = Bf
            self.hydroForce['ssRadf']['C'] = Cf*rho
    
    def setMassMatrix(self, rho, nlHydro):
        """
        Sets mass for the special cases of body at equilibrium or fixed
        Used by hydroForcePre
        
        """
        if self.mass == 'equilibrium':
            self.massCalcMethod = self.mass
            if nlHydro == 0:
                self.mass = self.hydroData['properties']['disp_vol'] * rho
            else:
                cg_tmp = self.hydroData['properties']['cg'][0]
                z = np.conj(np.transpose(np.array(self.bodyGeometry['center'])))[2] + cg_tmp[2]
                zr = [0 if x > 0 else x for x in z]
                area = np.conj(np.transpose(self.bodyGeometry['area']))[0]
                av = np.array([area,area,area])*-1*np.conj(np.transpose(self.bodyGeometry['norm']))
                tmp = rho*np.array([zr, zr, zr])*-1*av
                self.mass = sum(tmp[2])
        elif self.mass == 'fixed':
            self.massCalcMethod = self.mass
            self.mass = 999
            self.momOfInertia = [999, 999, 999]
        else:
            self.massCalcMethod = 'user'
