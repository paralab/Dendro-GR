D = ["alpha", "chi", "K", "Gt0", "Gt1", "Gt2", "beta0", "beta1", "beta2",
     "B0", "B1", "B2", "gt0", "gt1", "gt2", "gt3", "gt4", "gt5",
     "At0", "At1", "At2", "At3", "At4", "At5" ]

# second derivs required for RHS
DD = ["gt0", "gt1", "gt2", "gt3", "gt4", "gt5", "chi",
      "alpha", "beta0", "beta1", "beta2" ]

def arithmetic_intensity(d, m, k):
    """
    d stencil order (7 pt)
    m block width (1d) 
    k padding width 
    """
    overall_ai = 0

    for var in D:
        read  = m**3
        ops   = (2*d-1) * ((m-2*k)*(m**2) * 6)
        write = ((m-2*k)*(m**2)*6)
        if var in DD:
            ops   += (2*d-1)* ((m-2*k)*(m**2) * 6)
            read  += m**3 * 2
            write += ((m-2*k)*(m**2)*6)
        overall_ai += ops/(read + write)/8
        
        print("%s \t AI = %.2f"%(var, (ops/(read + write)/8)))

            
arithmetic_intensity(7, 13, 3)    