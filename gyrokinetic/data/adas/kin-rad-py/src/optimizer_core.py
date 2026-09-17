# created and designed by jaime caballero for my phd at differ and the gkeyll team
# to mimic original https://github.com/jRoeltgen/radiation_operator/tree/main by J. Roeltgen
# code written with some help of ai 
#
# optimizer_core.py
# this script contains the core logic for the non-linear optimization of the fit parameters
# script with functions to compute 
# safe_integrand: the integrand in eq. 12 in Roeltgen's paper
#### useful stuff #####
# c_const = 8*electron_charge/sqrt(pi)*(2*electron_charge/mass)^(xn(5)/2);
#######################
import math
import numpy as np
from scipy.integrate import quad
from scipy.optimize import minimize, Bounds, LinearConstraint
from scipy.optimize._numdiff import approx_derivative

# just in case matlab optimizer is being used (no longer needed)
# 1. STATE CONTAINER FOR MATLAB BRIDGE
class MatlabBridgeState:
    Te = None
    Target = None
    Weight = None

# 2. TOP-LEVEL MATLAB CALLBACK
def matlab_objective(x_matlab):
    params = [float(val) for val in list(x_matlab)]
    cost = objective_function(params, MatlabBridgeState.Te, MatlabBridgeState.Target, MatlabBridgeState.Weight)
    return float(cost), []

# 3. INITIALIZATION HOOK
def init_matlab_globals(te, target, w):
    MatlabBridgeState.Te = np.array(te).flatten()
    MatlabBridgeState.Target = np.array(target).flatten()
    MatlabBridgeState.Weight = float(w)

def safe_integrand(v_bar, Te, A_bar, alpha, beta, V0, gamma):
    # function to eq. 13 
    # at this stage, the optimizer is not in its final implementation 
    # and some asymptotic limits may need to be handled 

    # 1. Calculate the dimensionless velocity ratio
    x = v_bar / V0

    #alpha = alpha * 1e3 # Remember to scale alpha back up for the integral
    
    # 2. Handle absolute zero to prevent 0.0 ** negative_power
    if x == 0.0:
        return 0.0
        
    # 3. The Numerically stable fraction
    if x < 1.0:
        # Multiply top and bottom by x^alpha
        num = (alpha + beta) * (x ** alpha)
        den = beta + alpha * (x ** (alpha + beta))
        fraction = num / den
    else:
        # Multiply top and bottom by x^-beta
        num = (alpha + beta) * (x ** -beta)
        den = beta * (x ** -(alpha + beta)) + alpha
        fraction = num / den

    # fraction = (alpha + beta) / (beta*(x**(-alpha)) + alpha*(x**(beta)))
        
    # 4. Construct the full integrand from eq. 13
    term1 = (v_bar ** 4) / (Te ** 1.5)
    term2 = A_bar * fraction * (v_bar ** gamma)
    term3 = np.exp(-(v_bar ** 2) / Te)
    
    return term1 * term2 * term3

def objective_function(params, Te_data, Li_target_data, weight_w):
    """
    The equivalent of MATLAB's objective.m
    Calculates the sum of the squared weighted residuals.
    """
    # Unpack the parameters exactly as MATLAB does
    A_scaled, alpha, beta, V0, gamma = params
    #alpha = alpha_guess * 1e3
    
    total_squared_error = 0.0
    
    for j in range(len(Te_data)):

        Te_j = Te_data[j]
        y_j = Li_target_data[j] # This is already L_z * 1e30

        v_max = np.sqrt(30.0 * Te_j)
    
        # integral split in two parts, having V0 as the split point
        # to handle numerical instabilities around the peak of the integrand

        if V0 < v_max:

            val_1, _ = quad(safe_integrand, 0, V0, 
                            args=(Te_j, A_scaled, alpha, beta, V0, gamma),
                            epsabs=1e-8, epsrel=1e-8)
                            
            val_2, _ = quad(safe_integrand, V0, v_max, 
                            args=(Te_j, A_scaled, alpha, beta, V0, gamma),
                            epsabs=1e-8, epsrel=1e-8)
                            
            calcY = val_1 + val_2

        else:
            # If V0 is larger than the expected peak, we can just do a single integral to save time
            calcY, _ = quad(safe_integrand, 0, v_max, 
                            args=(Te_j, A_scaled, alpha, beta, V0, gamma),
                            epsabs=1e-8, epsrel=1e-8)

        # single integral 
        # calcY, _ = quad(safe_integrand, 0, np.inf, 
        #                 args=(Te_j, A_scaled, alpha, beta, V0, gamma),
        #                 epsabs=1e-8, epsrel=1e-8)
        # calcY = integral_val
        
        # MATLAB residual logic: ((calcY - y) * y^(weight_power - 1))^2
        residual = (calcY - y_j) * (y_j ** (weight_w - 1))

        #print(residual)
        
        # Add the squared residual to the total error
        total_squared_error += residual ** 2
        
    return total_squared_error

def run_single_optimization(initial_guess, Te_data, Li_target_data, weight_w, optimizer_choice='ipopt', eng=None):
    """
    Sets up the constraints, bounds, and runs the optimizer.
    Equivalent to the fmincon setup in radiation_operator.m.
    """
    # 1. Scale the target data internally
    # This keeps the scaling trick hidden from the user/main script
    # right now is done in main script for development purposes
    Bs = 1e30
    target_data_scaled = Li_target_data # * Bs
    
    # 2. Define constraints (Equivalent to MATLAB's xl and xu)
    # Order: [A_scaled, alpha, beta, V0, gamma]
    # Using the safe physics bounds from radiation_wrapper.m
    lower_bounds = [1e-12, 0.01, 0.08, 0.1, -20.0]
    upper_bounds = [np.inf, np.inf, 70.0,  80.0, 20.0]
    bounds = Bounds(lower_bounds, upper_bounds)
    
    # 3. Define Linear Constraints (Equivalent to MATLAB's A and b matrices)
    # The constraint matrix multiplies our parameter array:
    # [A_scaled, alpha, beta, V0, gamma]
    #
    # Row 1: 1*alpha + 1*gamma >= 0  (alpha + gamma > 0)
    # Row 2: -1*beta + 1*gamma <= 2  (gamma - beta < 2)
    constraint_matrix = [
        [0,  1,  0, 0, 1], 
        [0,  0, -1, 0, 1]  
    ]
    # Lower and upper limits for the two rows
    constraint_lb = [0.0, -np.inf]
    constraint_ub = [np.inf, 2.0]
    linear_constraint = LinearConstraint(constraint_matrix, constraint_lb, constraint_ub)
    
    # 4. Create a wrapper for the objective function
    # minimize() only wants to pass 'params', so we wrap the extra arguments
    def obj_wrapper(params):
        return objective_function(params, Te_data, target_data_scaled, weight_w)
    
    # Calculate MATLAB's sqrt(eps) for the FiniteDifferenceStepSize
    fd_step = np.sqrt(np.finfo(float).eps)
    fd_step = 1e-5
    
    def jac_wrapper(params):
        # Perfectly replicates 'FiniteDifferenceType', 'central' AND 'FiniteDifferenceStepSize', sqrt(eps)
        return approx_derivative(obj_wrapper, params, method='3-point', abs_step=fd_step)
    
    # 5. Run the Optimizer

    # ---------------------------------------------------------
    # 1: IPOPT (The open-source interior-point replica)
    # ---------------------------------------------------------
    if optimizer_choice == 'ipopt':

        lower_bounds = [1e-12, 0.01, 0.08, 0.1, -20.0]
        upper_bounds = [np.inf, np.inf, 70.0,  80.0, 20.0]
        bounds = Bounds(lower_bounds, upper_bounds)
        
        constraint_matrix = [
            [0,  1,  0, 0, 1], # alpha + gamma >= 0
            [0,  0, -1, 0, 1]  # gamma - beta <= 2
        ]
        linear_constraint = LinearConstraint(constraint_matrix, [0.0, -np.inf], [np.inf, 2.0])
        
        from cyipopt import minimize_ipopt

        return minimize_ipopt(
            obj_wrapper, x0=initial_guess, jac=jac_wrapper,
            bounds=bounds, constraints=[linear_constraint],
            options={'disp': 5, 'max_iter': 8000, 'tol': 1e-6, 'acceptable_tol': 1e-5}
        )
    
    elif optimizer_choice == 'trust-constr':
        matlab_sqrt_eps = float(np.sqrt(np.finfo(float).eps))
        
        result = minimize(
            obj_wrapper,
            x0=initial_guess,
            method='trust-constr',
            bounds=bounds,
            constraints=[linear_constraint],
            options={
                'disp': False, 
                'maxiter': 8000,
                'xtol': 1e-12,                # Matched to MATLAB StepTolerance
                'gtol': 1e-12,                # Matched to MATLAB OptimalityTolerance
                'barrier_tol': 1e-12,         # Strict interior-point barrier
                'finite_diff_rel_step': matlab_sqrt_eps # Matched to MATLAB StepSize
            }
        )
        return result

    # ---------------------------------------------------------
    # 2: SLSQP (SciPy Native)
    # ---------------------------------------------------------
    elif optimizer_choice == 'slsqp':
        # SLSQP setup from previous steps...

        # Calculate standard machine epsilon square root (~1.49e-08)
        matlab_sqrt_eps = float(np.sqrt(np.finfo(float).eps))

        result = minimize(
            obj_wrapper,
            x0=initial_guess,
            method='SLSQP', #method='trust-constr',
            bounds=bounds,
            constraints=[linear_constraint],
            options={
                'disp': True,         # Set to True if you want to see iterations in the terminal
                'maxiter': 8000,       # Max iterations
                'ftol': 1e-12,
                'eps': matlab_sqrt_eps,            # Step size for numerical gradient approximation (if needed)
                # ,         # Function value tolerance
                # 'xtol': 1e-12,         # Step tolerance
                # 'gtol': 1e-12          # Optimality tolerance
            }
        )

        return result

    # ---------------------------------------------------------
    # 3: MATLAB FMINCON (via pyfmincon)
    # ---------------------------------------------------------

    if optimizer_choice == 'fmincon':
        
        # this is only if we want to use the MATLAB bridge for optimization, otherwise we can just use the scipy minimize function directly
        import matlab
        from opt import fmincon

        # Force MATLAB to reload this Python module, clearing any broken cache
        eng.eval("mod = py.importlib.import_module('optimizer_core'); py.importlib.reload(mod);", nargout=0)

        # 1. Load data into the MATLAB workspace
        eng.workspace['te_mat'] = matlab.double(Te_data.tolist())
        eng.workspace['target_mat'] = matlab.double(target_data_scaled.tolist())
        eng.workspace['w_mat'] = float(weight_w)
        
        # 2. Initialize the Class state in MATLAB's embedded Python
        eng.eval("py.optimizer_core.init_matlab_globals(te_mat, target_mat, w_mat)", nargout=0)

        # Format bounds and constraints
        lb = [1e-12, 0.01, 0.01, 0.1, -20.0]
        ub = [math.inf, math.inf, 70.0,  80.0, 20.0]
        
        A = [
            [0.0, -1.0,  0.0, 0.0, -1.0], 
            [0.0,  0.0, -1.0, 0.0,  1.0]
        ]
        b = [[0.0], [2.0]]
        
        # Calculate standard machine epsilon square root (~1.49e-08)
        matlab_sqrt_eps = float(np.sqrt(np.finfo(float).eps))

        # The Exact 1:1 options mapping from the paper's MATLAB script
        options = {
            'Algorithm': 'interior-point',
            'Display': 'iter', 
            'Diagnostics': 'off',
            'FunValCheck': 'on',
            'MaxIterations': 8000.0,
            'MaxFunctionEvaluations': 12000.0,
            'FiniteDifferenceType': 'central',
            'FiniteDifferenceStepSize': matlab_sqrt_eps,
            'DiffMaxChange': 0.1,
            'StepTolerance': 1e-12,
            'FunctionTolerance': 1e-12, 
            'ConstraintTolerance': 1e-12,
            'OptimalityTolerance': 1e-18
        }
        
        # 4. PASS THE FUNCTION AS A STRING
        xopt, fopt, exitflag, output = fmincon(
            'optimizer_core.matlab_objective', # <- Must be a string reference
            initial_guess, 
            lb, 
            ub, 
            options=options, 
            A=A, 
            b=b, 
            eng=eng
        )
        
        # Package into standard return object
        class DummyResult: pass
        result = DummyResult()
        
        result.x = np.array(xopt).flatten()
        result.success = (exitflag > 0)
        
        return result

    
