'''
Integrator for several rigid bodies.
'''
import numpy as np
import math as m
import scipy.sparse.linalg as spla
from functools import partial
import copy

from quaternion import Quaternion
from stochastic_forcing import stochastic_forcing as stochastic
from mobility import mobility as mob
import utils

import scipy

class QuaternionIntegrator(object):
  '''
  Integrator that timesteps using deterministic forwars Euler scheme.
  '''  
  def __init__(self, bodies, Nblobs, scheme, tolerance = None, domain = 'single_wall'): 
    ''' 
    Init object 
    '''
    self.bodies = bodies
    self.Nblobs = Nblobs
    self.scheme = scheme
    self.mobility_bodies = np.empty((len(bodies), 3, 3))

    # Other variables
    self.get_blobs_r_vectors = None
    self.mobility_blobs = None
    self.force_torque_calculator = None
    self.calc_K_matrix_bodies = None
    self.linear_operator = None
    self.eta = None
    self.a = None
    self.velocities = None
    self.velocities_previous_step = None
    self.first_step = True
    self.kT = 0.0
    self.tolerance = 1e-08
    self.rf_delta = 1e-05
    self.invalid_configuration_count = 0
    self.det_iterations_count = 0
    self.stoch_iterations_count = 0
    self.domain = domain

    # Optional variables
    self.build_stochastic_block_diagonal_preconditioner = None
    self.periodic_length = None
    self.calc_slip = None
    self.calc_force_torque = None
    self.mobility_inv_blobs = None
    self.first_guess = None
    self.preconditioner = None
    self.mobility_vector_prod = None    
    if tolerance is not None:
      self.tolerance = tolerance
    return 

  def advance_time_step(self, dt, *args, **kwargs):
    '''
    Advance time step with integrator self.scheme
    '''
    return getattr(self, self.scheme)(dt, *args, **kwargs)
    

  def deterministic_forward_euler(self, dt, *args, **kwargs): 
    ''' 
    Take a time step of length dt using the deterministic forward Euler scheme. 
    The function uses gmres to solve the rigid body equations.

    The linear and angular velocities are sorted like
    velocities = (v_1, w_1, v_2, w_2, ...)
    where v_i and w_i are the linear and angular velocities of body i.
    ''' 
    while True: 
      # Call preprocess
      preprocess_result = self.preprocess(self.bodies)
      print('das me innat')
      
      rand_rhs = np.random.normal(0.0, 1.0, len(self.bodies) * 6) 

      # Solve mobility problem
      sol_precond = self.solve_mobility_problem(x0 = self.first_guess,noise = rand_rhs, save_first_guess = True, update_PC = self.update_PC, step = kwargs.get('step'))
      
      # Extract velocities
      velocities = np.reshape(sol_precond[3*self.Nblobs: 3*self.Nblobs + 3*len(self.bodies)], (len(self.bodies) * 3))

      # Update location orientation 
      for k, b in enumerate(self.bodies):
        b.location_new = b.location + velocities[3*k:3*k+3] * dt

      # Call postprocess
      postprocess_result = self.postprocess(self.bodies)

      # Check positions, if valid, return 
      if self.check_positions(new = 'new', old = 'current', update_in_success = True, domain = self.domain) is True:
        return

    return
      


  def solve_mobility_problem(self, RHS = None, noise = None, noise_FT = None, AB = None, x0 = None, save_first_guess = False, PC_partial = None, *args, **kwargs): 
    ''' 
    Solve the mobility problem using preconditioned GMRES. Compute 
    velocities on the bodies subject to active slip and enternal 
    forces-torques.

    The linear and angular velocities are sorted like
    velocities = (v_1, w_1, v_2, w_2, ...)
    where v_i and w_i are the linear and angular velocities of body i.
    ''' 
    while True: 
      System_size = self.Nblobs * 3 + len(self.bodies) * 3

      # Get blobs coordinates
      r_vectors_blobs = self.get_blobs_r_vectors(self.bodies, self.Nblobs)

      # If RHS = None set RHS = [slip, -force_torque]
      if RHS is None:
        # Calculate slip on blobs
        if self.calc_slip is not None:
          slip = self.calc_slip(self.bodies, self.Nblobs)
        else:
          slip = np.zeros((self.Nblobs, 3))
        # Calculate force-torque on bodies
        force_torque = np.zeros((self.Nblobs, 3)) ######self.force_torque_calculator(self.bodies, r_vectors_blobs)
        # Add noise to the force/torque
        if noise_FT is not None:
          force_torque += noise_FT
        # Set right hand side
        RHS = np.reshape(np.concatenate([slip, -force_torque]), (System_size))

      # Add noise to the slip
      if noise is not None:
        RHS[0:System_size] -= noise

      # Calculate K matrix
      K = self.calc_K_matrix_bodies(self.bodies, self.Nblobs)

      # Set linear operators 
      linear_operator_partial = partial(self.linear_operator, 
                                        bodies = self.bodies, 
                                        r_vectors = r_vectors_blobs, 
                                        eta = self.eta, 
                                        a = self.a, 
                                        K_bodies = K,
                                        periodic_length=self.periodic_length)
      A = spla.LinearOperator((System_size, System_size), matvec = linear_operator_partial, dtype='float64')
      
#      e_i = np.identity(System_size)
#      for k in range(System_size):
#		print A.matvec(e_i[:,k])

      # Set preconditioner 
      if PC_partial is None:
        PC_partial = self.build_block_diagonal_preconditioner(self.bodies, r_vectors_blobs, self.Nblobs, self.eta, self.a, *args, **kwargs)
      PC = spla.LinearOperator((System_size, System_size), matvec = PC_partial, dtype='float64')

      # Scale RHS to norm 1
      RHS_norm = np.linalg.norm(RHS)
      if RHS_norm > 0:
        RHS = RHS / RHS_norm

      # Solve preconditioned linear system
      counter = gmres_counter(print_residual = self.print_residual)
      (sol_precond, info_precond) = utils.gmres(A, RHS, x0=x0, tol=self.tolerance, M=PC, maxiter=1000, restart=60, callback=counter) 
      self.det_iterations_count += counter.niter

      if save_first_guess:
        self.first_guess = sol_precond  

      # Scale solution with RHS norm
      if RHS_norm > 0:
        sol_precond = sol_precond * RHS_norm
      else:
        sol_precond[:] = 0.0
      
      # Return solution
      return sol_precond


  def solve_mobility_problem_dense_algebra(self, *args, **kwargs): 
    ''' 
    Solve the mobility problem using dense algebra methods. Compute 
    velocities on the bodies subject to active slip and enternal 
    forces-torques.
    
    The linear and angular velocities are sorted like
    velocities = (v_1, w_1, v_2, w_2, ...)
    where v_i and w_i are the linear and angular velocities of body i.
    ''' 
    while True: 
      # Calculate slip on blobs
      if self.calc_slip is not None:
        slip = self.calc_slip(self.bodies, self.Nblobs)
      else:
        slip = np.zeros((self.Nblobs, 3))

      # Get blobs coordinates
      r_vectors_blobs = self.get_blobs_r_vectors(self.bodies, self.Nblobs)

      # Calculate mobility (M) at the blob level
      mobility_blobs = self.mobility_blobs(r_vectors_blobs, self.eta, self.a)

      # Calculate resistance at the blob level (use np.linalg.inv or np.linalg.pinv)
      resistance_blobs = np.linalg.inv(mobility_blobs)

      # Calculate force-torque on bodies
      force_torque = self.force_torque_calculator(self.bodies, r_vectors_blobs)

      # Calculate block-diagonal matrix K
      K = self.calc_K_matrix(self.bodies, self.Nblobs)

      # Add slip force = K^T * M^{-1} * slip
      force_torque -= np.reshape(np.dot(K.T,np.dot(resistance_blobs, np.reshape(slip, (3*self.Nblobs,1)))), force_torque.shape)
   
      # Calculate mobility (N) at the body level. Use np.linalg.inv or np.linalg.pinv
      mobility_bodies = np.linalg.pinv(np.dot(K.T, np.dot(resistance_blobs, K)), rcond=1e-14)

      # Compute velocities, return velocities and bodies' mobility
      return (np.dot(mobility_bodies, np.reshape(force_torque, 6*len(self.bodies))), mobility_bodies)


  def solve_mobility_problem_DLA(self, *args, **kwargs): 
    ''' 
    Solve the mobility problem using dense algebra methods. Compute 
    velocities on the bodies subject to active slip and enternal 
    forces-torques.
    
    The linear and angular velocities are sorted like
    velocities = (v_1, w_1, v_2, w_2, ...)
    where v_i and w_i are the linear and angular velocities of body i.
    ''' 
    while True: 
      # Calculate slip on blobs
      if self.calc_slip is not None:
        slip = self.calc_slip(self.bodies, self.Nblobs)
      else:
        slip = np.zeros((self.Nblobs, 3))

      # Get blobs coordinates
      r_vectors_blobs = self.get_blobs_r_vectors(self.bodies, self.Nblobs)

      # Calculate mobility (M) at the blob level
      mobility_blobs = self.mobility_blobs(r_vectors_blobs, self.eta, self.a)

      # Calculate resistance at the blob level (use np.linalg.inv or np.linalg.pinv)
      resistance_blobs = np.linalg.inv(mobility_blobs)

      # Calculate block-diagonal matrix K
      K = self.calc_K_matrix(self.bodies, self.Nblobs)
     
      # Calculate constraint force due to slip l = M^{-1}*slip
      force_slip = np.dot(K.T,np.dot(resistance_blobs, np.reshape(slip, (3*self.Nblobs,1))))
      
      # Calculate force-torque on bodies
      force_torque = self.force_torque_calculator(self.bodies, r_vectors_blobs)
      
      # Calculate RHS
      FT = np.reshape(force_torque, 6*len(self.bodies))
      FTS = FT + np.reshape(force_slip, 6*len(self.bodies))

      # Calculate mobility (N) at the body level. Use np.linalg.inv or np.linalg.pinv
      mobility_bodies = np.linalg.pinv(np.dot(K.T, np.dot(resistance_blobs, K)), rcond=1e-14)

      # Compute velocities
      return (np.dot(mobility_bodies, FTS), mobility_bodies, mobility_blobs, resistance_blobs, K, r_vectors_blobs)



  def check_positions(self, new = None, old = None, update_in_success = None, update_in_failure = None, domain = 'single_wall'):
    '''
    This function checks if the configuration is valid calling
    body.check_function. If necessary it updates the configuration
    of body.location and body.orientation.
    '''
    # Check positions, if valid return 
    valid_configuration = True
    if domain == 'single_wall':
      if new == 'current':
        for b in self.bodies:
          valid_configuration = b.check_function(b.location, b.orientation)
          if valid_configuration is False:
            self.invalid_configuration_count += 1
            print 'Invalid configuration number ', self.invalid_configuration_count
            break
      elif new == 'new':
        for b in self.bodies:
          valid_configuration = b.check_function(b.location_new, b.orientation_new)
          if valid_configuration is False:
            self.invalid_configuration_count += 1
            print 'Invalid configuration number ', self.invalid_configuration_count
            break

    # Update position if necessary
    if (valid_configuration is False) and (update_in_failure is True):
      if old == 'old':
        for b in self.bodies:
          np.copyto(b.location, b.location_old)
          b.orientation = copy.copy(b.orientation_old)
      if old == 'current':
        for b in self.bodies:
          b.location = b.location
          b.orientation = b.orientation          
    elif (valid_configuration is True) and (update_in_success is True):
      if new == 'current':
        for b in self.bodies:
          b.location = b.location
          b.orientation = b.orientation
      if new == 'new':
        for b in self.bodies:
          np.copyto(b.location, b.location_new)
          b.orientation = copy.copy(b.orientation_new)

    # Return true or false
    return valid_configuration


class gmres_counter(object):
  '''
  Callback generator to count iterations. 
  '''
  def __init__(self, print_residual = False):
    self.print_residual = print_residual
    self.niter = 0
  def __call__(self, rk=None):
    self.niter += 1
    if self.print_residual is True:
      if self.niter == 1:
        print 'gmres =  0 1'
      print 'gmres = ', self.niter, rk

      
