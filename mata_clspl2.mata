* version 1.1.0  20aug2025

version 15.1
clear all

loc RS        real scalar
loc RV        real colvector
loc RM        real matrix
loc SS        string scalar
loc TM        transmorphic matrix
loc CC        class
loc VV        void

loc version   "version 15.1:"

mata:
mata set matastrict on

`CC' CLSPL2                                                        /* class   */
{
	public:
		`RS' rcond, tolerance, iteration_limit, r, final, alpha, kappaC,     ///
		     kappaB, kappaA, rmsa, r2_partial, nrmse, nrmse_partial, seed
		`RV' b, zhat, z, y, rmsa_i, rmsa_dkappaC, rmsa_dkappaB,              ///
		     rmsa_dkappaA, rmsa_dnrmse, rmsa_dzhat, rmsa_dz, rmsa_dx,        ///
		     nrmse_ttest, z_lower, z_upper, y_lower, y_upper
		`RM' A, C_idx, Z, x, x_lower, x_upper
		`SS' distribution
		`TM' corr(), ttest()
		`VV' new(), canonize(), solve()
	private:
		`RS' cond()
		`RV' rnormal(), runiform(), rlaplace()
		`RM' uniqrows_tagindex()
}

`VV' CLSPL2::new()                                                 /* public  */
{
	/*
		                                                                        
		(void) new()                                                            
		                                                                        
		class constructor: default order and values of member variables         
		                                                                        
		This class initiates and stores the results of a two-stage estimation   
		process for constrained least-squares problems under a canonical matrix 
		structure. The approach is designed for underdetermined, ill-posed, or  
		structurally constrained systems.                                       
		                                                                        
		This class initiates and stores the results of a two-stage estimation   
		process for constrained least-squares problems under a canonical matrix 
		structure. The approach is designed for underdetermined, ill-posed, or  
		structurally constrained systems.                                       
		                                                                        
		Stage 1 estimates the solution z from the system A·z = b using a        
		pseudoinverse method. Depending on the projection matrix Z, this may    
		be the Moore–Penrose inverse (Z = I) or the Bott–Duffin inverse (Z ≠ I).
		If r > 1, the pseudoinverse solution is iteratively refined using an    
		updated slack matrix Q to improve numerical stability and feasibility.  
		                                                                        
		Stage 2 optionally refines the estimate using a convex optimization     
		step centered around the pseudoinverse estimate ẑ. Regularization is    
		governed by the pre-specified parameter α:                              
		    • α = 1      → Ridge (ℓ₂ norm).
		                                                                        
		The class also supports diagnostic routines such as Monte Carlo         
		hypothesis testing and row-wise deletion sensitivity analysis.          
		                                                                        
		Attributes:                                                             
		----                                                                    
		this.A : real matrix                                                    
		    Canonical design matrix with block structure A = [C | S; M | Q].    
		                                                                        
		this.C_idx : real matrix                                                
		    Pair of integers defining the row and column ranges of the C block  
		    inside A. Used for matrix slicing and partitioning.                 
		                                                                        
		this.b : real colvector                                                 
		    Right-hand side vector for the linear system A·z = b.               
		                                                                        
		this.Z : real matrix                                                    
		    Projection matrix for Bott–Duffin inversion. Must be symmetric and  
		    idempotent. Defaults to identity for Moore–Penrose.                 
		                                                                        
		rcond : real scalar                                                     
		    Regularization parameter for the Moore-Penrose and Bott-Duffin      
		    inverses.                                                           
		                                                                        
		this.tolerance : real scalar                                            
		    Convergence tolerance for NRMSE change between iterations.          
		                                                                        
		this.iteration_limit : real scalar                                      
		    Maximum number of iterations allowed in the refinement loop.        
		                                                                        
		this.zhat : real colvector                                              
		    Unregularized pseudoinverse estimate of z from Step 1.              
		                                                                        
		this.final : real scalar                                                
		    Whether to run the second-stage convex refinement. If 1, the        
		    estimate is corrected via convex programming.                       
		                                                                        
		this.alpha : real scalar                                                
		    Regularization parameter:                                           
		    - α = 1: Ridge (L2 norm).                                           
		                                                                        
		this.z : real colvector                                                 
		    Final estimate of z after correction (if final is 1). If skipped,   
		    z = zhat.                                                           
		                                                                        
		this.x : real matrix                                                    
		    Variable component extracted from z, reshaped into m × p.           
		                                                                        
		this.y : real colvector                                                 
		    Slack component of z representing inequality residuals.             
		                                                                        
		this.r : real scalar                                                    
		    Number of refinement iterations performed during Step 1. Iteration  
		    stops when NRMSE stabilizes or the limit is reached.                
		                                                                        
		this.kappaC : real scalar                                               
		    Condition number of the constraint block C (upper-left block in A). 
		                                                                        
		this.kappaB : real scalar                                               
		    Condition number of the projected estimator B^(r) = A·pinv(C),      
		    calculated during refinement.                                       
		                                                                        
		this.kappaA : real scalar                                               
		    Condition number of the full canonical matrix A^(r) after refinement
		    step r.                                                             
		                                                                        
		this.rmsa : real scalar                                                 
		    Total RMSA (Root Mean Square Adjustment) over all rows.             
		                                                                        
		this.rmsa_i : real colvector                                            
		    Change in RMSA caused by removing one row from [C | S] and          
		    re-estimating the CLSP solution.                                    
		                                                                        
		this.rmsa_dkappaC : real colvector                                      
		    Change in κ(C) caused by removing one row from [C | S] and          
		    re-estimating the CLSP solution.                                    
		                                                                        
		this.rmsa_dkappaB : real colvector                                      
		    Change in κ(B) caused by removing one row from [C | S] and          
		    re-estimating the CLSP solution.                                    
		                                                                        
		this.rmsa_dkappaA : real colvector                                      
		    Change in κ(A) caused by removing one row from [C | S] and          
		    re-estimating the CLSP solution.                                    
		                                                                        
		this.rmsa_dnrmse : real colvector                                       
		    Change in NRMSE caused by removing one row from [C | S] and         
		    re-estimating the CLSP solution.                                    
		                                                                        
		this.rmsa_dzhat : real colvector                                        
		    Change in zhat caused by removing one row from [C | S] and          
		    re-estimating the CLSP solution.                                    
		                                                                        
		this.rmsa_dz : real colvector                                           
		    Change in z caused by removing one row from [C | S] and             
		    re-estimating the CLSP solution.                                    
		                                                                        
		this.rmsa_dx : real colvector                                           
		    Change in x caused by removing one row from [C | S] and             
		    re-estimating the CLSP solution.                                    
		                                                                        
		this.r2_partial : real scalar                                           
		    R² statistic computed over the M block. Reflects partial            
		    goodness-of-fit in structured systems.                              
		                                                                        
		this.nrmse : real scalar                                                
		    Normalized RMSE between A·z and b, computed over full system.       
		                                                                        
		this.nrmse_partial : real scalar                                        
		    NRMSE computed only over M block rows.                              
		                                                                        
		this.nrmse_ttest : real colvector                                       
		    NRMSE samples generated via simulation under the null, used for     
		    empirical t-testing.                                                
		                                                                        
		this.z_lower : real colvector                                           
		    Lower bound of confidence band on z. Computed from κ(A) and residual
		    norm.                                                               
		                                                                        
		this.z_upper : real colvector                                           
		    Upper bound of confidence band on z. Symmetric to z_lower.          
		                                                                        
		this.x_lower : real matrix                                              
		    Lower bound of confidence band on x. Computed from κ(A) and residual
		    norm.                                                               
		                                                                        
		this.x_upper : real matrix                                              
		    Upper bound of confidence band on x. Symmetric to x_lower.          
		                                                                        
		this.y_lower : real colvector                                           
		    Lower bound of confidence band on y. Computed from κ(A) and residual
		    norm.                                                               
		                                                                        
		this.y_upper : real colvector                                           
		    Upper bound of confidence band on y. Symmetric to y_lower.          
		                                                                        
		this.seed : real scalar                                                 
		    Random seed used for reproducible Monte Carlo diagnostics.          
		                                                                        
		this.distribution : string scalar                                       
		    Function that generates random samples for simulation. Must accept a
		    single integer argument `n`.                                        
		                                                                        
	*/
	`TM' tmp, rc

	// default values of member variables                                       
	this.A               = /* design matrix, [C|S;M|Q] */ J(0,0,.)
	this.C_idx           = /* indices of the C block   */ J(1,2,.)
	this.b               = /* right-hand side          */ J(0,1,.)
	this.Z               = /* B-D subspace matrix      */ J(0,0,.)
	this.rcond           = /* regularization parameter */ .
	this.tolerance       =                                sqrt(epsilon(1))
	this.iteration_limit = /* limit  of iterations     */ 50
	this.r               = /* number of iterations     */ 0
	this.zhat            = /* first-step estimate      */ J(0,1,.)
	this.final           = /* inclusion of second step */ 1
	this.alpha           = /* regularization parameter */ 1.0
	this.z               = /* final solution           */ J(0,1,.)
	this.x               = /* variable component of z  */ J(0,0,.)
	this.y               = /* slack component of z     */ J(0,1,.)
	this.kappaC          = /* spectral κ() for C_canon */ .
	this.kappaB          = /* spectral κ() for B^(r)   */ .
	this.kappaA          = /* spectral κ() for A^(r)   */ .
	this.rmsa            = /* total RMSA               */ .
	this.rmsa_i          = /* list of row RMSA         */ J(0,1,.)
	this.rmsa_dkappaC    = /* list of Δκ(C)            */ J(0,1,.)
	this.rmsa_dkappaB    = /* list of Δκ(B)            */ J(0,1,.)
	this.rmsa_dkappaA    = /* list of Δκ(A)            */ J(0,1,.)
	this.rmsa_dnrmse     = /* list of ΔNRMSE           */ J(0,1,.)
	this.rmsa_dzhat      = /* list of Δzhat            */ J(0,1,.)
	this.rmsa_dz         = /* list of Δz               */ J(0,1,.)
	this.rmsa_dx         = /* list of Δx               */ J(0,1,.)
	this.r2_partial      = /* R^2   for the M block    */ .
	this.nrmse           = /* NRMSE for A              */ .
	this.nrmse_partial   = /* NRMSE for the M block    */ .
	this.nrmse_ttest     = /* list of NRMSE            */ J(0,1,.)
	this.z_lower         = /* lower confidence band    */ J(0,1,.)
	this.z_upper         = /* upper confidence band    */ J(0,1,.)
	this.x_lower         = /* lower confidence band    */ J(0,0,.)
	this.x_upper         = /* upper confidence band    */ J(0,0,.)
	this.y_lower         = /* lower confidence band    */ J(0,1,.)
	this.y_upper         = /* upper confidence band    */ J(0,1,.)
	this.seed            = /* Monte Carlo              */ 123456789
	this.distribution    = "normal"
}

`VV' CLSPL2::canonize(|`SS' problem, `RM' C, `RM' S, `RM' M, `RM' Q, `RV' b,    
									 `RS' m, `RS' p, `RS' i, `RS' j,            
									 `RS' zero_diagonal)           /* public  */
{
	/*
		                                                                        
		(void) canonize()                                                       
		                                                                        
		This method assembles the constraint matrix A from user-supplied or     
		internally generated components — C, S, M, and Q — and assigns the      
		corresponding right-hand side vector b. It is a required pre-step before
		solving a Convex Least Squares Programming (CLSP) problem.              
		                                                                        
		Depending on the specified problem type, it can generate allocation,    
		tabular matrix, or modular constraints and enforce optional diagonal    
		exclusions. All missing blocks are padded to ensure conformability.     
		                                                                        
		Parameters:                                                             
		----                                                                    
		problem : string scalar                                                 
		    Structural template for matrix construction. One of:                
		    - "ap"   or "tm" : allocation or tabular matrix problem.            
		    - "cmls" or "rp" : constrained modular least squares or RP-type.    
		    - ""     or other: General CLSP problems (user-defined C and/or M). 
		                                                                        
		C, S, M : real matrix                                                   
		    Blocks of the constraint matrix A = [C | S; M | Q].                 
		    If `C` and/or `M` are provided, the matrix A is constructed         
		    accordingly. If both are J(0,0,.) and A is not yet defined, an error
		    is raised.                                                          
		                                                                        
		Q : real matrix                                                         
		    Externally supplied residual slack matrix used to adjust inequality 
		    constraints in M. Required only when r > 1. Encodes the sign pattern
		    of residuals from the previous iteration and is used to construct   
		    the [C | S; M | Q] canonical form. Defaults to a conformable zero   
		    matrix on the first iteration.                                      
		                                                                        
		b : real colvector                                                      
		    Right-hand side vector. Must have as many rows as A. Required.      
		                                                                        
		m, p : real scalar                                                      
		    Dimensions of X ∈ ℝ^{m×p}, relevant for allocation problems ("ap"). 
		                                                                        
		i, j : real scalar, default = 1                                         
		    Grouping sizes for row and column sum constraints in AP problems.   
		                                                                        
		zero_diagonal : real scalar                                             
		    If 1, enforces structural zero diagonals via identity truncation.   
		                                                                        
		Returns                                                                 
		----                                                                    
		N/A, updates this.A, this.C_idx, this.b                                 
	*/
	`RS' M_flag, k, n_col
	`RM' row_groups, col_groups, M_diag
	`TM' tmp, rc
	i = i < . ? i : 1
	j = j < . ? j : 1

	// (b) Ensure the right-hand side is defined and set `this.b`               
	if (b == J(0,1,.)) {
		errprintf("Right-hand side vector b must be provided.\n")
		exit(198)
	} else this.b = colshape(b, 1)

	// (A) Option 1. AP (TM) problems with an optional zero diagonal            
	if (regexm(strlower(problem),   "ap|tm")) {
		if (m >= . | p >= .) {
			errprintf("Both m and p must be specified.\n")
			exit(198)
		}
		if (mod(m, i) != 0) {
			errprintf("m = " + strofreal(m) + " must be divisible by "         +
			          "i = " + strofreal(i) + ".\n")
			exit(503)
		}
		if (mod(p, j) != 0) {
			errprintf("p = " + strofreal(p) + " must be divisible by "         +
			          "j = " + strofreal(j) + ".\n")
			exit(503)
		}
		/* construct the C block using Kronecker product                      */
		row_groups = (I(trunc(m / i)) # J(1, i, 1)) # J(1, p, 1)
		col_groups =  J(1, m, 1) # (I(trunc(p / j)) # J(1, j, 1))
		if (C != J(0,0,.)) {
			if (cols(C) != cols(row_groups)) {
				errprintf("C must have " + strofreal(cols(row_groups))         +
				          " columns.\n")
				exit(503)
			}
		}
		C = (C != J(0,0,.) ? row_groups \ col_groups \ C                        
			               : row_groups \ col_groups                            
		)
		/* append an optional identity matrix to M, remove duplicates         */
		if (zero_diagonal < . & zero_diagonal) {
			M_flag = M == J(0,0,.)
			M_diag = J(min((m, p)),m * p,0)
			for (k = 1; k <= min((m, p)); k++) {
				M_diag[k, (k - 1) * p + k] = 1
			}
			M = this.uniqrows_tagindex(M_flag ? M_diag : M \ M_diag, tmp)
			this.b = b[1::rows(C)] \ (M_flag ? J(min((m, p)),1,0)               
				                             : (b[|rows(C)+1\.|]               \
				                                J(min((m, p)),1,0))[tmp])
		}
	}

	// (A) Option 2. CMLS and RP problems                                       
	if (regexm(strlower(problem), "cmls|rp")) {
		if (C == J(0,0,.) | M == J(0,0,.)) {
			errprintf("Both C and M must be provided.\n")
			exit(198)
		}
	}

	// (A) Option 3. General problems                                           
	if (C == J(0,0,.) & M == J(0,0,.)) {
		errprintf("At least one of C or M must be provided.\n")
		exit(198)
	}

	// (A) Convert missing blocks to conformable zero matrices                  
	n_col = C != J(0,0,.) ? cols(C) : cols(M)
	C     = C != J(0,0,.) ? C       : J(0,      n_col,0)
	M     = M != J(0,0,.) ? M       : J(0,      n_col,0)
	S     = S != J(0,0,.) ? S       : J(rows(C),0,    0)
	Q     = Q != J(0,0,.) ? Q       : J(rows(M),0,    0)
	if (rows(C) != rows(S)) {
		errprintf("C and S must have the same number of rows: "                +
		          strofreal(rows(C)) + " vs " + strofreal(rows(S)) + ".\n")
		exit(503)
	}
	if (cols(C) != cols(M)) {
		errprintf("C and M must have the same number of columns: "             +
		          strofreal(cols(C)) + " vs " + strofreal(cols(M)) + ".\n")
		exit(503)
	}

	// (A) Pad C and Q with zeros and set `this.A` and `this.C_idx`             
	this.A     = C,S,J(rows(C),cols(Q),0) \ M,J(rows(M),cols(S),0),Q
	this.C_idx =       rows(C),cols(C)
}

`VV' CLSPL2::solve(|`SS' problem, `RM' C, `RM' S, `RM' M, `RV' b,               
								  `RS' m, `RS' p, `RS' i, `RS' j,               
								  `RS' zero_diagonal,     `RS' r, `RM' Z,       
								  `RS' rcond,             `RS' tolerance,       
								  `RS' iteration_limit,                         
								  `RS' final,                                   
								  `RS' cond_tolerance)             /* public  */
{
	/*
		                                                                        
		(void) solve()                                                          
		                                                                        
		Solve the Convex Least Squares Programming (CLSP) problem.              
		                                                                        
		This method performs a two-step estimation:                             
		(1) a pseudoinverse-based solution using either the Moore–Penrose or    
		    Bott–Duffin inverse, optionally iterated for convergence;           
		(2) a convex-programming correction using Lasso, Ridge, or Elastic Net  
		regularization (if enabled).                                            
		                                                                        
		Parameters:                                                             
		----                                                                    
		problem : string scalar                                                 
		    Structural template for matrix construction. One of:                
		    - "ap"   or "tm" : allocation or tabular matrix problem.            
		    - "cmls" or "rp" : constrained modular least squares or RP-type.    
		    - ""     or other: General CLSP problems (user-defined C and/or M). 
		                                                                        
		C, S, M : real matrix                                                   
		    Blocks of the constraint matrix A = [C | S; M | Q].                 
		    If `C` and/or `M` are provided, the matrix A is constructed         
		    accordingly. If both are J(0,0,.) and A is not yet defined, an error
		    is raised.                                                          
		                                                                        
		b : real colvector                                                      
		    Right-hand side vector. Must have as many rows as A. Required.      
		                                                                        
		m, p : real scalar                                                      
		     Dimensions of X ∈ ℝ^{m×p}, relevant for allocation problems ("ap").
		                                                                        
		i, j : real scalar, default = 1                                         
		    Grouping sizes for row and column sum constraints in AP problems.   
		                                                                        
		zero_diagonal : real scalar                                             
		    If 1, enforces structural zero diagonals via identity truncation.   
		                                                                        
		r : real scalar, default = 1                                            
		    Number of refinement iterations for the pseudoinverse-based         
		    estimator. When `r > 1`, the slack block Q is updated iteratively to
		    improve feasibility in underdetermined or ill-posed systems.        
		                                                                        
		Z : real matrix                                                         
		    A symmetric idempotent matrix (projector) defining the subspace for 
		    Bott–Duffin pseudoinversion. If J(0,0,.), the identity matrix is    
		    used, reducing to the Moore–Penrose case.                           
		                                                                        
		rcond : real scalar                                                     
		    Regularization parameter for the Moore-Penrose and Bott-Duffin      
		    inverses, providing numerically stable inversion and ensuring       
		    convergence of singular values. If 1, an automatic tolerance equal  
		    to `tolerance` is applied. If |rcond| < 1, it specifies the relative
		    cutoff below which small singular values are treated as zero.       
		                                                                        
		tolerance : real scalar                                                 
		    Convergence tolerance for NRMSE change between iterations.          
		                                                                        
		iteration_limit : real scalar                                           
		    Maximum number of iterations allowed in the refinement loop.        
		                                                                        
		final : real scalar                                                     
		    If 1, an inverse problem is solved to refine zhat. The resulting    
		    solution zminimizes the L2 distance tozhat subject to Az = b (via a 
		    second pseudoinverse estimation in the absence of an optimizer).    
		                                                                        
		cond_tolerance : real scalar                                            
		    Singular-value cutoff for the custom condition number function. If  
		    cond_tolerance ≥ ., the implementation uses an internal relative    
		    cutoff of 1e-14.                                                    
		                                                                        
		Returns                                                                 
		----                                                                    
		N/A, updates this.A, this.C_idx, this.b, this.Z, this.tolerance,        
		             this.iteration_limit, this.r, this.zhat, this.final,       
		             this.z, this.x, this.y, this.kappaC, this.kappaB,          
		             this.kappaA, this.r2_partial, this.nrmse,                  
		             this.nrmse_partial, this.z_lower, this.z_upper,            
		             this.x_lower, this.x_upper, this.y_lower, this.y_upper     
	*/
	`RS' n_iter, nrmse_prev, Z_delta, dz
	`RV' b_M, residuals_M
	`RM' Q
	`TM' tmp, rc
	i     = i     < . ? i : 1
	j     = j     < . ? j : 1
	r     = r     < . ? r : 1
	final = final < 1 ? 0 : 1

	// (A), (b) Construct a conformable canonical form for the CLSP estimator   
	if ((C != J(0,0,.) | M != J(0,0,.)) | (m < . & p < .)) {
		this.canonize(problem, C, S, M, J(0,0,.), b, m, p, i, j, zero_diagonal)
	} else if (this.A == J(0,0,.)) {
		errprintf("At least one of C, M, m, or p must be provided.\n")
		exit(198)
	}
	if (rows(this.A) != rows(this.b)) {
		errprintf("The matrix A and vector b must have the same number of "    +
		          "rows: A has " + strofreal(rows(this.A)) + ", but b has "    +
		          strofreal(rows(this.b)) + ".\n")
		exit(503)
	}

	// (zhat) (Iterated if r > 1) first-step estimate                           
	if (r < 1) {
		errprintf("Number of refinement iterations r must be ≥ 1.\n")
		exit(198)
	}
	     if (Z      != J(0,0,.))   this.Z =          Z
	else if (this.Z == J(0,0,.))   this.Z =          I(cols(this.A))
	else                           this.Z =          this.Z[1::cols(this.A),    
															1..cols(this.A)]
	if (rcond           < .)       this.rcond           = rcond
	if (tolerance       < .)       this.tolerance       = tolerance
	if (iteration_limit < .)       this.iteration_limit = iteration_limit
																	 rc = .
	if      (rows(this.Z)   !=     cols(this.Z))                     rc = 503
	else if (rows(this.Z)   !=     cols(this.A)                                |
			 max(abs(this.Z -           this.Z')) >  this.tolerance            |
			 max(abs(this.Z * this.Z  - this.Z )) >  this.tolerance) rc = 503
	if (rc < .) {
		errprintf("Matrix Z must be symmetric, idempotent and match the "      +
		          "number of columns in A: expected ("                         +
		          strofreal(cols(this.A)) + ","                                +
		          strofreal(cols(this.A)) + "), got ("                         +
		          strofreal(rows(this.Z)) + ","                                +
		          strofreal(cols(this.Z)) + ").\n")
		exit(rc)
	}
	for (n_iter=1; n_iter <= (rows(this.A) > this.C_idx[1] ? r : 1); n_iter++) {
		// save NRMSE from the previous step, construct Q and Z                 
		if (n_iter > 1) {
				nrmse_prev = norm(this.b - this.A * this.zhat,   2)            /
							 sqrt(rows( this.b))                               /
							 sqrt(mean((this.b :- mean(this.b)):^2))
				Q = diag(-sign(this.b - this.A *this.zhat)[|this.C_idx[1]+1\.|])
				this.canonize(problem,  C, S, M, Q, colshape(this.b, 1),        
									    m, p, i, j, zero_diagonal)
				Z_delta                 = cols(this.A) - rows(this.Z)
				if (Z_delta > 0) this.Z = this.Z, J(rows(this.Z),Z_delta,0)    \
										  J(Z_delta,cols(this.Z),0),I(Z_delta)
		}
		// solve via the Bott–Duffin inverse                                    
		this.zhat  = (pinv(this.Z    * quadcross(this.A,this.A )* this.Z, .,    
				        abs(this.rcond)>=1 & abs(this.rcond)<. ? -this.tolerance
				     : -abs(this.rcond))  * this.Z    * this.A')* this.b
		this.nrmse = norm(this.b - this.A * this.zhat,   2)                    /
					 sqrt(rows( this.b))                                       /
					 sqrt(mean((this.b :- mean(this.b)):^2))
		// break on convergence                                                 
		this.r     = n_iter
		if (n_iter > 1) {
			if (abs(this.nrmse - nrmse_prev) < this.tolerance                  |
				n_iter                       > this.iteration_limit) break
		}
	}
	if (any(this.zhat :>= .)) {
		this.zhat = J(0,1,.)
		errprintf("Pseudoinverse estimate zhat failed.\n")
		exit(504)
	}

	// (z) Final solution (if available), or set this.z = this.zhat             
	this.final = final
	if (final) {
		this.z = this.zhat + this.A' * pinv(this.A * this.A', .,                
										abs(this.rcond)>=1 &  abs(this.rcond)<. 
										?  -this.tolerance : -abs(this.rcond)) *
							(this.b - this.A * this.zhat)
		if (all(this.z :< .))
			this.nrmse = norm(this.b - this.A * this.z,      2)                /
						 sqrt(rows( this.b))                                   /
						 sqrt(mean((this.b :- mean(this.b)):^2))
		else this.z = this.zhat
	}   else this.z = this.zhat

	// (x), (y) Variable and slack components of z                              
	this.x = colshape(this.z[1::this.C_idx[2]], p < . ? p : 1)
	this.y = (cols(this.A) >    this.C_idx[2]
			  ? this.z[|this.C_idx[2]+1\.|] : J(0,1,0))

	// (kappaC), (kappaB), (kappaA) Condition numbers                           
	this.kappaC = this.cond(              this.A[1::this.C_idx[1],.],  "e")
	this.kappaB = this.cond(this.A * pinv(this.A[1::this.C_idx[1],.]), "e")
	this.kappaA = this.cond(              this.A,                      "e")

	// (r2_partial), (nrmse_partial) M-block-based statistics                   
	if (rows(this.A) > this.C_idx[1]) {
		M                  = this.A[this.C_idx[1]+1::rows(this.A),
									1..this.C_idx[2]]
		b_M                = this.b[|rows(this.b)-rows(M)+1\.|]
		residuals_M        = b_M - M * colshape(this.x, 1)
		this.r2_partial    = 1 - norm(residuals_M,      2)^2                   /
								 norm(b_M :- mean(b_M), 2)^2
		this.nrmse_partial = norm(residuals_M,          2)                     /
							 sqrt(rows( b_M))                                  /
							 sqrt(mean((b_M :- mean(b_M)):^2))
	}

	// (z_lower), (z_upper) Condition-weighted confidence band                  
	dz           = this.kappaA * norm(this.b - this.A * this.z, 2)             /
								 norm(this.b,                   2)
	this.z_lower = this.z :* (1 - dz)
	this.z_upper = this.z :* (1 + dz)

	// (x_lower), (x_upper), (y_lower), (y_upper)                               
	this.x_lower = colshape(this.z_lower[1::this.C_idx[2]], p < . ? p : 1)
	this.x_upper = colshape(this.z_upper[1::this.C_idx[2]], p < . ? p : 1)
	this.y_lower = (cols(this.A) >  this.C_idx[2]
					? this.z_lower[|this.C_idx[2]+1\.|] : J(0,1,0))
	this.y_upper = (cols(this.A) >  this.C_idx[2]
					? this.z_upper[|this.C_idx[2]+1\.|] : J(0,1,0))
}

`TM' CLSPL2::corr(|`RS' reset, `RS' threshold)                     /* public  */
{
	/*
		                                                                        
		(transmorphic matrix) corr()                                            
		                                                                        
		Compute the structural correlogram of the CLSP constraint system.       
		                                                                        
		This method performs a row-deletion sensitivity analysis on the         
		canonical constraint matrix [C | S], denoted as C_canon, and evaluates  
		the marginal effect of each constraint row on numerical stability,      
		angular alignment, and estimator sensitivity.                           
		                                                                        
		For each row i in C_canon, it computes:                                 
		- The Root Mean Square Alignment (RMSA_i) with all other rows j ≠ i.    
		- The change in condition numbers κ(C), κ(B), and κ(A) when row i is    
		  deleted.                                                              
		- The effect on estimation quality: changes in NRMSE, zhat, z, and x.   
		                                                                        
		Additionally, it computes the total RMSA statistic across all rows,     
		summarizing the overall angular alignment of the constraint block.      
		                                                                        
		Parameters:                                                             
		----                                                                    
		reset : real scalar, default = 0                                        
		    If 1, forces recomputation of all diagnostic values.                
		                                                                        
		threshold : real scalar, default = 0                                    
		    If positive, limits the output to constraints with                  
		    RMSA_i ≥ threshold.                                                 
		                                                                        
		Returns                                                                 
		----                                                                    
		associative array                                                       
		An associative array containing per-row diagnostic values:              
		{                                                                       
		    "constraint"   : [1, 2, ..., k],  # 1-based indices,                
		    "rmsa_i"       : real colvector of RMSA_i values,                   
		    "rmsa_dkappaC" : real colvector of Δκ(C) after deleting row i,      
		    "rmsa_dkappaB" : real colvector of Δκ(B) after deleting row i,      
		    "rmsa_dkappaA" : real colvector of Δκ(A) after deleting row i,      
		    "rmsa_dnrmse"  : real colvector of ΔNRMSE after deleting row i,     
		    "rmsa_dzhat"   : real colvector of Δ‖zhat‖_2 after deleting row i,  
		    "rmsa_dz"      : real colvector of Δ‖z‖_2 after deleting row i,     
		    "rmsa_dx"      : real colvector of Δ‖x‖_F|2 after deleting row i    
		}                                                                       
	*/
	`RS' k, p, i, j
	`RV' norms
	`RM' C_canon
	`TM' tmp, rc
	`CC' CLSPL2           scalar obj
	`CC' AssociativeArray scalar dict
	reset       = reset       < 1 | reset       >= . ? 0  : 1
	threshold   = threshold   < 0 | threshold   >  1 ? 0  : threshold

	// (RMSA) Total RMSA                                                        
	if (this.rmsa >= . | reset) {
		k       = this.C_idx[1]
		p       = this.C_idx[2]
		C_canon = this.A[1::k,.]
		norms   = sqrt(rowsum(C_canon:^2))
		rc      = 0
		for (i  = 1; i <= k - 1; i++) {
			for (j = i +  1;  j <= k; j++) {
				rc = rc + ((C_canon[i,.] * C_canon[j,.]')                      /
						   (norms[i] * norms[j]))^2
			}
		}
		this.rmsa = sqrt(2 / k / (k - 1) * rc)
	}

	// (RMSA) Constraint-wise RMSA, changes in condition numbers, and GoF       
	if (rows(this.rmsa_i) != this.C_idx[1] | reset) {
		obj               = CLSPL2()
		k                 = this.C_idx[1]
		p                 = this.C_idx[2]
		C_canon           = this.A[1::k,.]
		norms             = sqrt(rowsum(C_canon:^2))
		this.rmsa_i       = J(k,1,.)
		this.rmsa_dkappaC = J(k,1,.)
		this.rmsa_dkappaB = J(k,1,.)
		this.rmsa_dkappaA = J(k,1,.)
		this.rmsa_dnrmse  = J(k,1,.)
		this.rmsa_dzhat   = J(k,1,.)
		this.rmsa_dz      = J(k,1,.)
		this.rmsa_dx      = J(k,1,.)
		for (i = 1; i <= k; i++) {
			obj.A                = select(this.A, (1::rows(this.A)) :!= i)
			obj.b                = select(this.b, (1::rows(this.b)) :!= i)
			obj.C_idx            = k - 1,p
			obj.solve("", J(0,0,.), J(0,0,.), J(0,0,.), J(0,1,.),               
					  ., ., ., ., .,                                            
					  this.r, this.Z, this.rcond, this.tolerance,               
					  this.iteration_limit, this.final)
			rc                   = 0
			for (j = 1; j <= k; j++) {
				if (j != i) rc   = rc + ((C_canon[i,.] * C_canon[j,.]')        /
										 (norms[i] * norms[j]))^2
			}
			this.rmsa_i[i]       = sqrt(1 / (k - 1) * rc)
			this.rmsa_dkappaC[i] = obj.kappaC - this.kappaC
			this.rmsa_dkappaB[i] = obj.kappaB - this.kappaB
			this.rmsa_dkappaA[i] = obj.kappaA - this.kappaA
			this.rmsa_dnrmse[i]  = obj.nrmse  - this.nrmse
			this.rmsa_dzhat[i]   = norm(obj.zhat, 2) - norm(this.zhat, 2)
			this.rmsa_dz[i]      = norm(obj.z,    2) - norm(this.z,    2)
			this.rmsa_dx[i]      = norm(obj.x,  cols(obj.x ) > 1 ? 0 : 2)      -
								   norm(this.x, cols(this.x) > 1 ? 0 : 2)
		}
	}

	// Return the correlogram                                                   
	rc   = selectindex(this.rmsa_i :>= threshold)
	dict = AssociativeArray()
	dict.put("constraint",   (1::rows(this.rmsa_i))[rc])    // 1-based indexing  
	dict.put("rmsa_i",       this.rmsa_i[rc])
	dict.put("rmsa_dkappaC", this.rmsa_dkappaC[rc])
	dict.put("rmsa_dkappaB", this.rmsa_dkappaB[rc])
	dict.put("rmsa_dkappaA", this.rmsa_dkappaA[rc])
	dict.put("rmsa_dnrmse",  this.rmsa_dnrmse[rc])
	dict.put("rmsa_dzhat",   this.rmsa_dzhat[rc])
	dict.put("rmsa_dz",      this.rmsa_dz[rc])
	dict.put("rmsa_dx",      this.rmsa_dx[rc])
	return(dict)
}

`TM' CLSPL2::ttest(|`RS' reset,        `RS' sample_size, `RS' seed,             
					`SS' distribution, `RS' partial,                            
					`RS' simulate)                                 /* public  */
{
	/*
		                                                                        
		(transmorphic matrix) ttest()                                           
		                                                                        
		Perform bootstrap or Monte Carlo t-tests on the NRMSE statistic from the
		CLSP estimator.                                                         
		                                                                        
		This function either (a) resamples residuals via a nonparametric        
		bootstrap to generate an empirical NRMSE sample, or (b) produces        
		synthetic right-hand side vectors `b` from a user-defined or default    
		distribution and re-estimates the model. It tests whether the observed  
		NRMSE significantly deviates from the null distribution of resampled or 
		simulated NRMSE values.                                                 
		                                                                        
		Parameters:                                                             
		----                                                                    
		reset : real scalar, default = 0                                        
		    If 1, forces recomputation of the NRMSE null distribution.          
		                                                                        
		sample_size : real scalar, default = 50                                 
		    Size of the Monte Carlo simulated sample under H₀.                  
		                                                                        
		seed : real scalar                                                      
		    Optional random seed to override the default.                       
		                                                                        
		distribution : string scalar                                            
		    Distribution for generating synthetic b vectors. One of: "normal",  
		    "uniform", "laplace". Defaults to standard normal.                  
		                                                                        
		partial : real scalar, default = 0                                      
		    If 1, runs the t-test on the partial NRMSE: during simulation, the  
		    C-block entries are preserved and the M-block entries are simulated.
		                                                                        
		simulate : bool, default = 0                                            
		    If 1, performs a parametric Monte Carlo simulation by generating    
		    synthetic right-hand side vectors `b`. If 0 (default), executes a   
		    nonparametric bootstrap procedure on residuals without              
		    re-estimation.                                                      
		                                                                        
		Returns                                                                 
		----                                                                    
		associative array                                                       
		A dictionary with test results and null distribution statistics:        
		{                                                                       
		    "p_one_left"  : P(nrmse ≤ null mean),                               
		    "p_one_right" : P(nrmse ≥ null mean),                               
		    "p_two_sided" : 2-sided t-test p-value,                             
		    "nrmse"       : observed value,                                     
		    "mean_null"   : mean of null distribution,                          
		    "std_null"    : std of null distribution                            
		}                                                                       
	*/
	`RS' i, mean_null, std_null, t_stat, p_left, p_right, p_two
	`RV' residuals, b
	`TM' tmp, rc
	`CC' CLSPL2           scalar obj
	`CC' AssociativeArray scalar dict
	reset       = reset       < 1 | reset       >= . ? 0  : 1
	sample_size = sample_size < 1 | sample_size >= . ? 50 : sample_size
	partial     = partial     < 1 | partial     >= . ? 0  : 1
	simulate    = simulate    < 1 | simulate    >= . ? 0  : 1

	// Set the seed, RNG configuration, and distribution                        
	if (seed                  <  . )   this.seed          = seed
	                             rseed(this.seed)
	if (strtrim(distribution) != "")   this.distribution  = distribution
	if (!  regexm(strlower(this.distribution), "normal|uniform|laplace")) {
		errprintf("Unsupported distribution: " + this.distribution + ".\n")
		exit(198)
	}

	// (t-test) Bootstrap-resampled or simulated NRMSE distribution under H0    
	if (rows(this.nrmse_ttest) != sample_size | reset) {
		if (partial & rows(this.A) == this.C_idx[1]) {
			errprintf("No M-block present in A; falling back to full NRMSE "   +
			          "t-test.\n")
			partial = 0
		}
		this.nrmse_ttest = J(sample_size,1,.)
		// (re)generate a nonparametric bootstrap sample                        
		if (! simulate) {
			residuals = this.A * this.z - this.b
			for (i  = 1; i <= sample_size; i++) {
				rc = (tmp=! partial ? residuals :  residuals[|                  
													       this.C_idx[1]+1\.|])[
													   runiformint(rows(tmp),1, 
													              1,rows(tmp))]
				b  = (tmp=! partial ? this.b :     this.b[|this.C_idx[1]+1\.|])
				this.nrmse_ttest[i] = norm(rc,    2)                           /
									  sqrt(rows( b))                           /
									  sqrt(mean((b :- mean(b)):^2))
			}
		// (re)generate a parametric Monte Carlo sample                         
		} else {
			obj       = CLSPL2()
			obj.A     = this.A
			obj.C_idx = this.C_idx
			for (i  = 1; i <= sample_size; i++) {          // simulate b_M only 
				       if (regexm(strlower(this.distribution), "normal" )) {
					obj.b = (! partial ? this.rnormal(rows(this.b     ))
									   : this.b[1::this.C_idx[1]]              \
									     this.rnormal(rows(this.b[|             
									               this.C_idx[1]+1\.|])))
				} else if (regexm(strlower(this.distribution), "uniform")) {
					obj.b = (! partial ? this.runiform(rows(this.b    ))
									   : this.b[1::this.C_idx[1]]              \
									     this.runiform(rows(this.b[|            
									               this.C_idx[1]+1\.|])))
				} else if (regexm(strlower(this.distribution), "laplace")) {
					obj.b = (! partial ? this.rlaplace(rows(this.b    ))
									   : this.b[|1\this.C_idx[1]|]             \
									     this.rlaplace(rows(this.b[|            
									               this.C_idx[1]+1\.|])))
				}
				obj.solve("", J(0,0,.), J(0,0,.), J(0,0,.), J(0,1,.),           
						  ., ., ., ., .,                                        
						  this.r, this.Z, this.rcond, this.tolerance,           
						  this.iteration_limit, this.final)
				this.nrmse_ttest[i] =  ! partial ? obj.nrmse : obj.nrmse_partial
			}
		}
	}

	// Return the t-test                                                        
	mean_null = mean(this.nrmse_ttest)
	std_null  = sqrt(variance(this.nrmse_ttest))
	t_stat    = ((! partial ? this.nrmse : this.nrmse_partial) - mean_null)    /
				 (std_null  / sqrt(sample_size))
	p_left    = ttail(sample_size - 1,    -t_stat)
	p_right   = 1 - p_left
	p_two     = ttail(sample_size - 1, abs(t_stat)) * 2
	dict      = AssociativeArray()
	dict.put("p_one_left",  p_left)
	dict.put("p_one_right", p_right)
	dict.put("p_two_sided", p_two)
	dict.put("nrmse",       ! partial ? this.nrmse : this.nrmse_partial)
	dict.put("mean_null",   mean_null)
	dict.put("std_null",    std_null)
	return(dict)
}

`RM' CLSPL2::uniqrows_tagindex(`RM' X, `TM' idx)                   /* private */
{
	/*
		--                                                                      
		M = `RM' uniqrows_tagindex(M, tmp)                                      
		                      -  -                                              
		unique rows with first-occurrence indices                               
		                                                                        
		returns: (uniqrows(X))                                           (k x p)
	*/
	`RS' r
	`RM' U
	`RV' hit

	U   = uniqrows(X)
	idx = J(0, 1, .)
	for (r  = 1; r <= rows(U); r++) {
		hit = selectindex(rowsum(X :!= U[r,.]) :== 0)
		idx = idx \ hit[1]
	}
	return(U)
}

`RS' CLSPL2::cond(`RM' X, |`SS' n, `RS' tol)                       /* private */
{
	/*
		--                                                                      
		kappaA = `RS' cond(X)                                                   
		                      -  -                                              
		Compute a custom condition number using the singular values of X        
		                                                                        
		returns: (c)                                                     (1 x 1)
	*/
	`RS' smax
	`RV' s
	n    = n != "" ? strlower(substr(strtrim(n), 1, 1)) : "e"
	tol  = tol < . ? tol : 1e-14

	if (! rows(X) | ! cols(X))  return(.)
	s    = select((s= svdsv(X)), s :> tol * (smax=max(s)))
	if (! rows(s))              return(.)
	if (n ==  "e")              return(smax            / min(s)                )
	else                        return(sqrt(sum(s:^2)) * sqrt(sum((1 :/ s):^2)))
}

`RV' CLSPL2::rnormal(`RS'  n)                                      /* private */
{
	/*
		--                                                                      
		this.b = `RM' this.rnormal(n)                                           
		                      -  -                                              
		normal (Gaussian) pseudorandom variates                                 
		                                                                        
		returns: (X)                                                     (n x 1)
	*/
	return(rnormal(n, 1, 0, 1))
}

`RV' CLSPL2::runiform(`RS'  n)                                     /* private */
{
	/*
		--                                                                      
		this.b = `RM' this.runiform(n)                                          
		                      -  -                                              
		matrix of uniformly distributed pseudorandom numbers                    
		                                                                        
		returns: (X)                                                     (n x 1)
	*/
	return(runiform(n, 1, -1, 1))
}

`RV' CLSPL2::rlaplace(`RS'  n)                                     /* private */
{
	/*
		--                                                                      
		this.b = `RM' this.rlaplace(n)                                          
		                      -  -                                              
		Laplace pseudorandom variates                                           
		                                                                        
		returns: (X)                                                     (n x 1)
	*/
	return(rlaplace(n, 1, 0, 1))
}
end

version 15.1: lmbuild lclspl2.mlib, replace size(12)
