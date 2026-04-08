version 15.1

// CERTIFICATION FILE for CLSPL2()
discard
clear
clear mata
mata: CLSP = CLSPL2()

// AP problem with zero diagonal
mata: C = J(0,0,.)
mata: S = J(0,0,.)
mata: M = J(0,0,.)
mata: b = 10 \ 20 \ 15 \ 15
mata: m = 2
mata: p = 2
mata: i = 1
mata: j = 1
mata: zerodiagonal = 1
mata: r = 1
mata: CLSP.solve("ap", C, S, M, b, m, p, i, j, zerodiagonal, r)
/* estimation results */
mata: CLSP.A
mata: CLSP.C_idx
mata: CLSP.b
mata: CLSP.zhat
mata: CLSP.z
mata: CLSP.x
mata: CLSP.y
mata: CLSP.kappaC
mata: CLSP.kappaB
mata: CLSP.kappaA
mata: CLSP.r2_partial
mata: CLSP.nrmse
mata: CLSP.nrmse_partial
/* structural correlogram */
mata: R = CLSP.corr()
mata: R.get("constraint")
mata: R.get("rmsa_i")
/* t-tests on the NRMSE */
mata: T = CLSP.ttest()
mata: T.get("p_two_sided")
mata: T.get("nrmse")
