# Nombre: Sergio Martinez Olivera
# Grupo: DG MAT-INF
# para esta actividad necesitamos la primera (cuerpos_finitos.py) completa
import cuerpos_finitos as cf
import sys
import time

# input: fpx -> anillo_fp_x
# input: f -> polinomio (objeto opaco creado por fpx)
# input: g -> polinomio (objeto opaco creado por fpx)
# output: f*g calculado usando el método de Karatsuba
def fp_x_mult_karatsuba(fpx, f, g):

    if (len(f) == 1 or len(g) == 1):
        return fpx.mult(f, g)
    else:
        m = max(len(f), len(g))
        
        f_izq = f[:m//2]
        f_der = f[m//2:]
        g_izq = g[:m//2]
        g_der = g[m//2:]
        
        if not f_der: f_der = fpx.cero()
        if not g_der: g_der = fpx.cero()

        c0 = fp_x_mult_karatsuba(fpx, f_izq, g_izq)
        c2 = fp_x_mult_karatsuba(fpx, f_der, g_der)
        
        c1_1 = fp_x_mult_karatsuba(fpx, fpx.suma(f_izq, f_der), fpx.suma(g_izq, g_der))
        c1_2 = fpx.suma(c0, c2)
        c1 = fpx.suma(c1_1, fpx.inv_adit(c1_2))
        
        return tuple(fpx.suma(fpx.suma(c0, (fpx.fp.cero(),)* (m//2) + c1), (fpx.fp.cero(),)*2*(m//2) + c2))

# añadimos esta función a la clase (sin sobreescribir la que ya teníamos)
cf.anillo_fp_x.mult_fast = fp_x_mult_karatsuba

# input: fp -> cuerpo_fp
# input: n >= 1 (int)
# input: a -> tupla de longitud n de elementos de fp (primera columna de una
#    matriz de Toeplitz inferior T de nxn)
# input: b -> tupla de longitud n de elementos de fp (vector)
# output: T*b -> tupla de longitud n de elementos de fp (vector)
# se debe utilizar fp_x_mult_karatsuba internamente
def fp_toep_inf_vec(fp, n, a, b):
    res = fp_x_mult_karatsuba(cf.anillo_fp_x(fp), a, b)[:n]
    expected_len = max(1, len(a) + len(b) - 1)
    if len(res) < expected_len:
        res += (fp.cero(),) * (expected_len - len(res))
    return res[:n]
    '''T = [[0 for i in range(n)] for j in range(n)]
    for i in range(n):
        T[i] = [a[i-j] if i>=j else 0 for j in range(n)]
    return tuple([sum(T[i][j]*b[j]) for j in range(n)] for i in range(n))'''
    

# input: fp -> cuerpo_fp
# input: n >= 1 (int)
# input: a -> tupla de longitud n de elementos de fp (primera fila de una
#    matriz de Toeplitz superior T de nxn)
# input: b -> tupla de longitud n de elementos de fp (vector)
# output: T*b -> tupla de longitud n de elementos de fp (vector)
# se debe utilizar fp_x_mult_karatsuba internamente
def fp_toep_sup_vec(fp, n, a, b):
    res = fp_x_mult_karatsuba(cf.anillo_fp_x(fp), a[::-1], b)
    expected_len = max(1, len(a) + len(b) - 1)
    if len(res) < expected_len:
        res += (fp.cero(),) * (expected_len - len(res))
    return res[n - 1:]
    ''' T = [[0 for i in range(n)] for j in range(n)]
    for i in range(n):
        T[i] = [a[n - 1 + i - j] if j>=i else 0 for j in range(n)]
    return tuple([sum(T[i][j]*b[j]) for j in range(n)] for i in range(n))'''
    
# input: fp -> cuerpo_fp
# input: n >= 1 (int)
# input: a -> tupla de longitud 2*n-1 de elementos de fp (primera fila de una
#    matriz de Toeplitz completa T de nxn seguida de la primera columna
#    excepto el elemento de la esquina)
# input: b -> tupla de longitud n de elementos de fp (vector)
# output: T*b -> tupla de longitud n de elementos de fp (vector)
# se debe utilizar fp_x_mult_karatsuba internamente
def fp_toep_vec(fp, n, a, b):
    aux1 = fp_toep_inf_vec(fp, n, (a[0],) + a[n:], b)
    aux2 = fp_toep_sup_vec(fp, n, (fp.cero(),) + a[1:n], b)
    return tuple(fp.suma(aux1[i], aux2[i]) for i in range(n))

# input: fp -> cuerpo_fp
# input: n >= 1 (int)
# input: a -> tupla de longitud n de elementos de fp (primera columna de una
#    matriz de Toeplitz inferior T de nxn)... suponemos a[0] != 0
# output: primera columna de T^(-1) -> tupla de longitud n de elementos de
#    fp (vector)
# utilizar un método recursivo que "divida el problema a la mitad"
# recordar que T^(-1) es también una matriz de Toeplitz inferior
def fp_toep_inf_inv(fp, n, a):
    if n == 1:
        return (fp.inv_mult(a[0]),)
    else:
        k = n // 2
        a1 = a[:k]
        
        T1_inv = fp_toep_inf_inv(fp, k, a1)
        
        a_row = a[k:0:-1] 
        a_col = a[k+1:n] 
        

        if len(a_row) < k:
            a_row += (fp.cero(),) * (k - len(a_row))
        if len(a_col) < k - 1:
            a_col += (fp.cero(),) * (k - 1 - len(a_col))

        a_2 = a_row + a_col
        T_aux = tuple(fp.inv_adit(elem) for elem in fp_toep_inf_vec(fp, k, T1_inv, fp_toep_vec(fp, k, a_2, T1_inv)))
        
        return T1_inv + T_aux

# input: fp -> cuerpo_fp
# input: n >= 1 (int)
# input: a -> tupla de longitud n de elementos de fp (primera fila de una
#    matriz de Toeplitz superior T de nxn)... suponemos a[0] != 0
# output: primera fila de T^(-1) -> tupla de longitud n de elementos de
#    fp (vector)
# utilizar un método recursivo que "divida el problema a la mitad"
# recordar que T^(-1) es también una matriz de Toeplitz superior
def fp_toep_sup_inv(fp, n, a):
    if n == 1:
        return (fq.inv_mult(a[0]),)
    else:
        k = n // 2
        a1 = a[:k]
        
        T1_inv = fp_toep_sup_inv(fp, k, a1)
        

        a_row = a[k+1:n] 
        a_col = a[k:0:-1] 
        
        if len(a_row) < k:
            a_row += (fp.cero(),) * (k - len(a_row))
        if len(a_col) < k - 1:
            a_col += (fp.cero(),) * (k - 1 - len(a_col))

        a_2 = a_row + a_col
        T_aux = tuple(fp.inv_adit(elem) for elem in fp_toep_sup_vec(fp, k, T1_inv, fp_toep_vec(fp, k, a_2, T1_inv)))
        
        return T1_inv + T_aux

# input: fpx -> anillo_fp_x
# input: f -> polinomio (objeto opaco creado por fpx)
# input: g -> polinomio no nulo (objeto opaco creado por fpx)
# output: q -> cociente
# output: r -> resto
# se cumple que f = g*q+r, r=0 o deg(r)<deg(g)
# reformular el problema en términos de matrices de Toeplitz y luego usar
# las funciones de arriba para obtener q y r
def fp_x_divmod(fpx, f, g):
    # Grados
    m = len(f) - 1
    n = len(g) - 1

    if m < n:
        # q = 0, r = f
        return (fpx.cero(), f ) 

    # Grado del cociente q
    k = m - n
    N = k + 1 
    
    # N debe ser una potencia de 2 para fp_toep_inf_inv
    N_pad = 1
    while N_pad < N:
        N_pad *= 2
    # --- FIN DE LA CORRECCIÓN ---

    g_rev = g[::-1] # (g_n, g_{n-1}, ..., g_0)
    a = g_rev[:N]
    # Rellenar con ceros si g_rev es más corto que N_pad
    if len(a) < N_pad:
        a = a + (fpx.fp.cero(),) * (N_pad - len(a))

    f_rev = f[::-1] # (f_m, f_{m-1}, ..., f_0)
    b = f_rev[:N]
    # Rellenar con ceros si b es más corto que N_pad
    if len(b) < N_pad:
        b = b + (fpx.fp.cero(),) * (N_pad - len(b))
    
    # 3. Resolver T * q_rev = b  =>  q_rev = T^{-1} * b
    
    # 3a. Obtener la 1ª columna de T^{-1} (usando el tamaño N_pad)
    t_inv_col = fp_toep_inf_inv(fpx.fp, N_pad, a)
    
    # 3b. Calcular q_rev (usando el tamaño N_pad)
    q_rev_pad = fp_toep_inf_vec(fpx.fp, N_pad, t_inv_col, b)
    
    # Recortar al tamaño original N
    q_rev = q_rev_pad[:N]
    
    # 4. Obtener q
    # Invertimos q_rev = (q_k, ..., q_0) para obtener q = (q_0, ..., q_k)
    q = q_rev[::-1]
    q = fpx.quitar_ceros(q) # Limpiar por si acaso

    # 5. Calcular el resto: r = f - g*q
    # Usamos la multiplicación rápida
    g_q = fpx.mult_fast(tuple(g), tuple(q))
    r = fpx.suma(tuple(f), fpx.inv_adit(tuple(g_q)))

    return (q, r)
    

# añadimos esta función a la clase (sin sobreescribir la que ya teníamos)
cf.anillo_fp_x.divmod_fast = fp_x_divmod

# input: fp -> cuerpo_fp
# input: g -> elemento del grupo multiplicativo fp* de orden n (objeto opaco)
# input: k >= 0 tal que n = 2**k divide a p-1
# input: a -> tupla de longitud n de elementos de fp
# output: DFT_{n,g}(a) -> tupla de longitud n de elementos de fp
# utilizar el algoritmo de Cooley-Tuckey
def fp_fft(fp, g, k, a):
    n = 2**k
    if k == 0:
        return (a[0],)
    
    pol_par = [a[2*i] for i in range(n//2)]
    pol_impar = [a[2*i+1] for i in range(n//2)]
    t_par = fp_fft(fp, fp.pot(g, 2), k-1, pol_par)
    t_impar = fp_fft(fp, fp.pot(g, 2), k-1, pol_impar)
    t = [fp.cero()]*n
    pot = fp.uno()
    for i in range(n//2):
        aux = fp.mult(t_impar[i], pot)
        t[i] = fp.suma(t_par[i], aux)
        t[i + n//2] = fp.suma(t_par[i], fp.inv_adit(aux))
        pot = fp.mult(pot, g)
    return tuple(t)


# input: fp -> cuerpo_fp
# input: g -> elemento del grupo multiplicativo fp* de orden n (objeto opaco)
# input: k >= 0 tal que n = 2**k divide a p-1
# input: a -> tupla de longitud n de elementos de fp
# output: IDFT_{n,g}(a) -> tupla de longitud n de elementos de fp
# recordar que IDFT_{n,g} = n^(-1) * DFT_{n,g^(-1)}
def fp_ifft(fp, g, k, a):
    n = 2**k
    return tuple(fp.mult(fp.inv_mult(n), x) for x in fp_fft(fp, fp.inv_mult(g), k, a))

# input: fqx -> anillo_fq_x
# input: f -> polinomio (objeto opaco creado por fqx)
# input: g -> polinomio (objeto opaco creado por fqx)
# output: f*g calculado usando el método de Karatsuba
def fq_x_mult_karatsuba(fqx, f, g):
    if (len(f) == 1 or len(g) == 1):
        return fqx.mult(f, g)
    else:
        m = max(len(f), len(g))
        
        f_izq = f[:m//2]
        f_der = f[m//2:]
        g_izq = g[:m//2]
        g_der = g[m//2:]
        
        if not f_der: f_der = fqx.cero()
        if not g_der: g_der = fqx.cero()

        c0 = fq_x_mult_karatsuba(fqx, f_izq, g_izq)
        c2 = fq_x_mult_karatsuba(fqx, f_der, g_der)
        
        c1_1 = fq_x_mult_karatsuba(fqx, fqx.suma(f_izq, f_der), fqx.suma(g_izq, g_der))
        c1_2 = fqx.suma(c0, c2)
        c1 = fqx.suma(c1_1, fqx.inv_adit(c1_2))
        
        return fqx.suma(fqx.suma(c0, fqx.cero()* (m//2) + c1), fqx.cero()*2*(m//2) + c2)

    

# añadimos esta función a la clase (sin sobreescribir la que ya teníamos)
cf.anillo_fq_x.mult_fast = fq_x_mult_karatsuba

# input: fq -> cuerpo_fq
# input: n >= 1 (int)
# input: a -> tupla de longitud n de elementos de fq (primera columna de una
#    matriz de Toeplitz inferior T de nxn)
# input: b -> tupla de longitud n de elementos de fq (vector)
# output: T*b -> tupla de longitud n de elementos de fq (vector)
# se debe utilizar fq_x_mult_karatsuba internamente
def fq_toep_inf_vec(fq, n, a, b):
    res = fq_x_mult_karatsuba(cf.anillo_fq_x(fq), a, b)
    expected_len = max(1, len(a) + len(b) - 1)
    if len(res) < expected_len:
        res += (fq.cero(),) * (expected_len - len(res))
    return res[:n]
# input: fq -> cuerpo_fq
# input: n >= 1 (int)
# input: a -> tupla de longitud n de elementos de fq (primera fila de una
#    matriz de Toeplitz superior T de nxn)
# input: b -> tupla de longitud n de elementos de fq (vector)
# output: T*b -> tupla de longitud n de elementos de fq (vector)
# se debe utilizar fq_x_mult_karatsuba internamente
def fq_toep_sup_vec(fq, n, a, b):
    res = fq_x_mult_karatsuba(cf.anillo_fq_x(fq), a[::-1], b)
    expected_len = max(1, len(a) + len(b) - 1)
    if len(res) < expected_len:
        res += (fq.cero(),) * (expected_len - len(res))
    return res[n - 1:]

# input: fq -> cuerpo_fq
# input: n >= 1 (int)
# input: a -> tupla de longitud 2*n-1 de elementos de fq (primera fila de una
#    matriz de Toeplitz completa T de nxn seguida de la primera columna
#    excepto el elemento de la esquina)
# input: b -> tupla de longitud n de elementos de fq (vector)
# output: T*b -> tupla de longitud n de elementos de fq (vector)
# se debe utilizar fq_x_mult_karatsuba internamente
def fq_toep_vec(fq, n, a, b):
    aux1 = fq_toep_inf_vec(fq, n, (a[0],) + a[n:], b)
    aux2 = fq_toep_sup_vec(fq, n, (fq.cero(),) + a[1:n], b)
    return tuple(fq.suma(aux1[i], aux2[i]) for i in range(n))


# input: fq -> cuerpo_fq
# input: n >= 1 (int)
# input: a -> tupla de longitud n de elementos de fq (primera columna de una
#    matriz de Toeplitz inferior T de nxn)... suponemos a[0] != 0
# output: primera columna de T^(-1) -> tupla de longitud n de elementos de
#    fq (vector)
# utilizar un método recursivo que "divida el problema a la mitad"
# recordar que T^(-1) es también una matriz de Toeplitz inferior
def fq_toep_inf_inv(fq, n, a):
    if n == 1:
        return (fq.inv_mult(a[0]),)
    else:
        k = n // 2
        a1 = a[:k]
        
        T1_inv = fq_toep_inf_inv(fq, k, a1)
        
        a_row = a[k:0:-1] 
        a_col = a[k+1:n] 
        
        if len(a_row) < k:
            a_row += (fq.cero(),) * (k - len(a_row))
        if len(a_col) < k - 1:
            a_col += (fq.cero(),) * (k - 1 - len(a_col))

        a_2 = a_row + a_col
        T_aux = tuple(fq.inv_adit(elem) for elem in fq_toep_inf_vec(fq, k, T1_inv, fq_toep_vec(fq, k, a_2, T1_inv)))
        
        return T1_inv + T_aux

# input: fq -> cuerpo_fq
# input: n >= 1 (int)
# input: a -> tupla de longitud n de elementos de fq (primera fila de una
#    matriz de Toeplitz superior T de nxn)... suponemos a[0] != 0
# output: primera fila de T^(-1) -> tupla de longitud n de elementos de
#    fq (vector)
# utilizar un método recursivo que "divida el problema a la mitad"
# recordar que T^(-1) es también una matriz de Toeplitz superior
def fq_toep_sup_inv(fq, n, a):
    if n == 1:
        return (fq.inv_mult(a[0]),)
    else:
        k = n // 2
        a1 = a[:k]
        
        T1_inv = fq_toep_sup_inv(fq, k, a1)
        
        a_row = a[k+1:n] 
        a_col = a[k:0:-1] 
        
        if len(a_row) < k:
            a_row += (fq.cero(),) * (k - len(a_row))
        if len(a_col) < k - 1:
            a_col += (fq.cero(),) * (k - 1 - len(a_col))

        a_2 = a_row + a_col
        T_aux = tuple(fq.inv_adit(elem) for elem in fq_toep_sup_vec(fq, k, T1_inv, fq_toep_vec(fq, k, a_2, T1_inv)))
        
        return T1_inv + T_aux

# input: fqx -> anillo_fq_x
# input: f -> polinomio (objeto opaco creado por fqx)
# input: g -> polinomio no nulo (objeto opaco creado por fqx)
# output: q -> cociente
# output: r -> resto
# se cumple que f = g*q+r, r=0 o deg(r)<deg(g)
# reformular el problema en términos de matrices de Toeplitz y luego usar
# las funciones de arriba para obtener q y r
def fq_x_divmod(fqx, f, g):
    # Grados
    m = len(f) - 1
    n = len(g) - 1

    if m < n:
        # q = 0, r = f
        return (fqx.cero(), f ) 

    # Grado del cociente q
    k = m - n
    N = k + 1 

    # --- INICIO DE LA CORRECCIÓN ---
    # N debe ser una potencia de 2 para fq_toep_inf_inv
    N_pad = 1
    while N_pad < N:
        N_pad *= 2

    g_rev = g[::-1] # (g_n, g_{n-1}, ..., g_0)
    a = g_rev[:N]
    # Rellenar con ceros si g_rev es más corto que N_pad
    if len(a) < N_pad:
        a = a + (fqx.fq.cero(),) * (N_pad - len(a))

    f_rev = f[::-1] # (f_m, f_{m-1}, ..., f_0)
    b = f_rev[:N]
    # Rellenar con ceros si b es más corto que N_pad
    if len(b) < N_pad:
        b = b + (fqx.fq.cero(),) * (N_pad - len(b))
    
    # 3. Resolver T * q_rev = b  =>  q_rev = T^{-1} * b
    
    # 3a. Obtener la 1ª columna de T^{-1}
    t_inv_col = fq_toep_inf_inv(fqx.fq, N_pad, a)
    
    # 3b. Calcular q_rev
    q_rev_pad = fq_toep_inf_vec(fqx.fq, N_pad, t_inv_col, b)
    
    # Recortar al tamaño original N
    q_rev = q_rev_pad[:N]
    
    # 4. Obtener q
    # Invertimos q_rev = (q_k, ..., q_0) para obtener q = (q_0, ..., q_k)
    q = q_rev[::-1]
    q = fqx.quitar_ceros(q) # Limpiar por si acaso

    # 5. Calcular el resto: r = f - g*q
    # Usamos la multiplicación rápida
    g_q = fqx.mult_fast(tuple(g), tuple(q))
    r = fqx.suma(tuple(f), fqx.inv_adit(tuple(g_q)))

    return (q, r)

# añadimos esta función a la clase (sin sobreescribir la que ya teníamos)
cf.anillo_fq_x.divmod_fast = fq_x_divmod

# input: fq -> cuerpo_fq
# input: g -> elemento del grupo multiplicativo fq* de orden n (objeto opaco)
# input: k >= 0 tal que n = 2**k divide a q-1
# input: a -> tupla de longitud n de elementos de fq
# output: DFT_{n,g}(a) -> tupla de longitud n de elementos de fq
# utilizar el algoritmo de Cooley-Tuckey
def fq_fft(fq, g, k, a):
    n = 2**k
    if k == 0:
        return (a[0],)
    
    pol_par = [a[2*i] for i in range(n//2)]
    pol_impar = [a[2*i+1] for i in range(n//2)]
    t_par = fq_fft(fq, fq.pot(g, 2), k-1, pol_par)
    t_impar = fq_fft(fq, fq.pot(g, 2), k-1, pol_impar)
    t = [fq.cero()]*n
    pot = fq.uno()
    for i in range(n//2):
        aux = fq.mult(t_impar[i], pot)
        t[i] = fq.suma(t_par[i], aux)
        t[i + n//2] = fq.suma(t_par[i], fq.inv_adit(aux))
        pot = fq.mult(pot, g)
    return tuple(t)


# input: fq -> cuerpo_fq
# input: g -> elemento del grupo multiplicativo fq* de orden n (objeto opaco)
# input: k >= 0 tal que n = 2**k divide a p-1
# input: a -> tupla de longitud n de elementos de fq
# output: IDFT_{n,g}(a) -> tupla de longitud n de elementos de fq
# recordar que IDFT_{n,g} = n^(-1) * DFT_{n,g^(-1)}
def fq_ifft(fq, g, k, a):
    n = 2**k
    return tuple(fq.mult(fq.inv_mult(fq.elem_de_tuple((n,))), x) for x in fq_fft(fq, fq.inv_mult(g), k, a))

