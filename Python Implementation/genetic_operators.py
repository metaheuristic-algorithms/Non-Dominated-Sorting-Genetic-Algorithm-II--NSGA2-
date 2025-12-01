import numpy as np
import random

def genetic_operators(p1_pos, p2_pos, var_min, var_max, pc, eta_c, pm, eta_m):
    """
    Performs SBX Crossover and Polynomial Mutation
    """
    n_var = len(p1_pos)
    c1 = p1_pos.copy()
    c2 = p2_pos.copy()
    
    # SBX Crossover
    # Full SBX implementation (commented out for clarity)
    """
    if random.random() <= pc:
        for j in range(n_var):
            if random.random() <= 0.5:
                if abs(p1_pos[j] - p2_pos[j]) > 1e-6: # Avoid div by zero
                    if p1_pos[j] < p2_pos[j]:
                        y1, y2 = p1_pos[j], p2_pos[j]
                    else:
                        y1, y2 = p2_pos[j], p1_pos[j]
                    
                    yl, yu = var_min[j], var_max[j]
                    rand_val = random.random()
                    
                    beta = 1.0 + (2.0 * (y1 - yl) / (y2 - y1))
                    alpha = 2.0 - beta**-(eta_c + 1.0)
                    
                    if rand_val <= (1.0 / alpha):
                        betaq = (rand_val * alpha)**(1.0 / (eta_c + 1.0))
                    else:
                        betaq = (1.0 / (2.0 - rand_val * alpha))**(1.0 / (eta_c + 1.0))
                        
                    c1[j] = 0.5 * ((y1 + y2) - betaq * (y2 - y1))
                    c2[j] = 0.5 * ((y1 + y2) + betaq * (y2 - y1))
                    
                else:
                   # Standard average if too close
                   beta = 0 # Dummy
            else:
               # Simple random exchange if not doing SBX on this var
               pass 
           """    
               
    # Simplified SBX (Matches Matlab code logic more closely)
    # The block above is full SBX, but here is the simpler version from your Matlab code:
    if random.random() <= pc:
        for j in range(n_var):
            u = random.random()
            if u <= 0.5:
                beta = (2 * u)**(1 / (eta_c + 1))
            else:
                beta = (1 / (2 * (1 - u)))**(1 / (eta_c + 1))
                
            val1 = 0.5 * ((1 + beta) * p1_pos[j] + (1 - beta) * p2_pos[j])
            val2 = 0.5 * ((1 - beta) * p1_pos[j] + (1 + beta) * p2_pos[j])
            c1[j] = val1
            c2[j] = val2

    # Polynomial Mutation
    for child in [c1, c2]:
        for j in range(n_var):
            if random.random() <= pm: 
                u = random.random()
                if u < 0.5:
                    delta = (2 * u)**(1 / (eta_m + 1)) - 1
                else:
                    delta = 1 - (2 * (1 - u))**(1 / (eta_m + 1))
                
                child[j] += delta * (var_max[j] - var_min[j])

    # Clipping to bounds
    c1 = np.clip(c1, var_min, var_max)
    c2 = np.clip(c2, var_min, var_max)
    
    return c1, c2