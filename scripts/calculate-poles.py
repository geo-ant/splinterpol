# algorithm 9 in the B&M paper
def compute_z_transform_coefficients(n: uint):

    if n < 2:
        raise Exception("spline degree must be >= 2")
    
    ntilde = n // 2
    # these are b^(n) (k) in the paper and
    # d^(n) (k)
    bn = [0] * (ntilde+2)
    dn = [0] * (ntilde+1)

    bn[0] = 1;
    dn[0] = 1/2;
        
    for m in range(1,n):
        mtilde = m // 2
        # these are the b^(m) (k) and d^(m) (k)
        # temporary arrays in the B&M paper
        bm = [0] * (ntilde+1)
        dm = [0] * (ntilde+1)
        for k in range(0,mtilde+1):
            # the last line is a hack to get d^(m)_{-1} = d^(m)_0, see
            # appendix A.1 in the paper
            bm[k] = ( (m+1)/2 + k ) * dn[k]+ ( (m+1)/2 - k)* dn[k-1 if k > 0 else 0]
            bm[k] /= m
            dm[k] = ( (m+2)/2 + k) * bn[k+1] + (m/2 - k) * bn[k]
            dm[k] /= m

        # reassign the temporaries to the original arrays
        for k in range(0,mtilde+1):
            bn[k] = bm[k]
            dn[k] = dm[k]

    for k in range(0,ntilde+1):
        bn[k] = ( (n+1)/2 + k) * dn[k] + ( (n+1)/2 - k)* dn[k-1 if k > 0 else 0]
        bn[k] /= n

    # now mirror and expand the coefficients such that
    # coeff[k] is the coefficient for power z^k in the
    # polynomial \tilde{B}^(n)(z) = coeff[k]*z^k
    # in the B&M article appendix A.1

    coeff = [0] * (n+1 if n % 2 == 0 else n)

    for k in range(0,ntilde+1):
        coeff[k] = bn[ntilde-k]
    for k in range(ntilde+1,2*ntilde+1):
        coeff[k] = coeff[2*ntilde-k]
        
    return coeff

coeff = compute_z_transform_coefficients(5)
print(coeff)
