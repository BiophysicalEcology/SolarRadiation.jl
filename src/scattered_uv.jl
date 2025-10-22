function allocate_scattered_uv()
    (;
        CHX  = zeros(101),
        CHY  = zeros(101),
        AMU  = zeros(101),
        X1   = zeros(101),
        Y1   = zeros(101),
        X2   = zeros(101),
        Y2   = zeros(101),
        AIL  = zeros(101),
        GAMR = zeros(101),
        GAML = zeros(101),
        chandrasekhar_xy_buffers = init_chandrasekhar_xy_buffers(),
    )
end

scattered_uv(TAU1::Float64) = scattered_uv!(allocate_scattered_uv(), TAU1)

function scattered_uv!(buffers, TAU1::Float64)
    # Large arrays (mutable, normal)
    (; CHX, CHY, AMU, X1, Y1, X2, Y2, AIL, GAMR, GAML, chandrasekhar_xy_buffers) = buffers

    # Small fixed-size arrays (use StaticArrays)
    AI  = @MVector zeros(30)

    # Set up AMU array
    AMU[1] = 0.0
    for I in 2:101
        AMU[I] = 0.01 * (I - 1)
    end


    # Compute X1, Y1 using chandrasekhar_xy
    CFA = SVector((0.75, -0.75, 0.0))
    # CHX_, CHY_, _ = chandrasekhar_xy!(chandrasekhar_xy_buffers, TAU1, collect(CFA), 111)
    CHX_, CHY_, _ = chandrasekhar_xy!(chandrasekhar_xy_buffers, TAU1, CFA, 111)
    X1 .= CHX_
    Y1 .= CHY_

    # Compute X2, Y2 using chandrasekhar_xy
    CFA = SVector((0.375, -0.375, 0.0))
    # CHX_, CHY_, _ = chandrasekhar_xy!(chandrasekhar_xy_buffers, TAU1, collect(CFA), 0)
    CHX_, CHY_, _ = chandrasekhar_xy!(chandrasekhar_xy_buffers, TAU1, CFA, 0)
    X2 .= CHX_
    Y2 .= CHY_

    # Compute AIL (quadrature weights)
    AIL[1] = 0.01 / 3.0
    CNU1 = 4.0 * AIL[1]
    CNU2 = 2.0 * AIL[1]
    for I in 2:2:100
        AIL[I] = CNU1
        AIL[I+1] = CNU2
    end
    AIL[101] = AIL[1]

    # Scalar accumulators
    xa1 = xa2 = xa3 = xa4 = 0.0
    xb1 = xb2 = xb3 = xb4 = xb5 = xb6 = xb7 = xb8 = 0.0

    for I in 1:101
        a  = AMU[I]
        a2 = a * a
        a3 = a2 * a

        c1 = AIL[I] * X1[I] * a
        xa1 += c1
        xa2 += c1 * a

        c2 = AIL[I] * Y1[I] * a
        xa3 += c2
        xa4 += c2 * a

        c3 = AIL[I] * X2[I]
        xb1 += c3
        xb2 += c3 * a
        xb3 += c3 * a2
        xb4 += c3 * a3

        c4 = AIL[I] * Y2[I]
        xb5 += c4
        xb6 += c4 * a
        xb7 += c4 * a2
        xb8 += c4 * a3
    end

    # Fill AI vector
    AI[1]  = xb1 + xb5 - 8.0 / 3.0
    AI[2]  = xb2 + xb6
    AI[3]  = xb3 + xb7
    AI[4]  = xb1 - xb5 - 8.0 / 3.0
    AI[5]  = xb2 - xb6
    AI[6]  = xb3 - xb7
    AI[7]  = xb4 - xb8
    AI[8]  = xa1 + xa3
    AI[9]  = xa2 + xa4
    AI[10] = xa1 - xa3
    AI[11] = xa2 - xa4

    AI[12] = (AI[1] - AI[3]) / ((AI[4] - AI[6]) * TAU1 + 2.0 * (AI[5] - AI[7]))
    AI[13] = 1.0 / (AI[4] * AI[10] - AI[5] * AI[11])
    AI[14] = 1.0 / (AI[1] * AI[8] - AI[2] * AI[9] -
                    2.0 * AI[12] * (AI[5] * AI[8] - AI[4] * AI[9]))
    AI[15] = 2.0 * (AI[8] * AI[10] - AI[9] * AI[11])
    AI[16] = AI[13] * AI[15]
    AI[17] = AI[14] * AI[15]

    CNU1 = 0.5 * (AI[16] - AI[17])
    CNU2 = 0.5 * (AI[16] + AI[17])

    AI[15] = AI[13] * (AI[5] * AI[8] - AI[4] * AI[9])
    AI[16] = AI[14] * (AI[2] * AI[10] - AI[1] * AI[11] -
                       2.0 * AI[12] * (AI[4] * AI[10] - AI[5] * AI[11]))

    AI[15] = AI[13] * (AI[2] * AI[10] - AI[1] * AI[11])
    AI[16] = AI[14] * (AI[5] * AI[8] - AI[4] * AI[9])
    CU3   = 0.5 * (AI[15] - AI[16])
    CU4   = 0.5 * (AI[15] + AI[16])

    AI[15] = AI[14] * (AI[1] * AI[8] - AI[2] * AI[9])
    SBAR  = 1.0 - 0.375 * AI[12] * (AI[4] - AI[6]) *
                   ((CNU2 - CNU1) * AI[8] + (CU4 - CU3) * AI[2] - AI[15] * AI[6])

    AI[20] = 0.375 * AI[12] * (CNU2 - CNU1) * (AI[4] - AI[6])
    AI[21] = 0.375 * AI[12] * (AI[4] - AI[6])
    AI[22] = AI[21] * (CU4 - CU3)
    AI[23] = AI[21] * AI[15]

    for I in 1:101
        GAML[I] = AI[20] * (X1[I] + Y1[I])
        GAMR[I] = AI[22] * (X2[I] + Y2[I]) - AMU[I] * AI[23] * (X2[I] - Y2[I])
    end

    return GAMR, GAML, SBAR
end

function init_chandrasekhar_xy_buffers()
    arrays = (;
        PSI = zeros(101),
        AMU = zeros(101),
        XA = zeros(101),
        XB = zeros(101),
        FNPP = zeros(101),
        FNPN = zeros(101),
        FNC0 = zeros(101),
        FNC1 = zeros(101),
        FNX = zeros(101),
        FNY = zeros(101),
        FNW = zeros(101),
        FMC0 = zeros(101),
        FMC1 = zeros(101),
        XD = zeros(101),
        XE = zeros(101),
        CHXA = zeros(101),
        CHYA = zeros(101),
        CHX = zeros(101),
        CHY = zeros(101),
    )
    return arrays
end

"""
    chandrasekhar_xy(TAU1::Float64, CFA::Vector{Float64}, NCASE::Int) 
        -> (CHX::Vector{Float64}, CHY::Vector{Float64}, nomitr::Int)

Compute Chandrasekhar's X and Y functions for radiative transfer.

# Description
This routine evaluates the X- and Y-functions of Chandrasekhar using
double precision arithmetic. The method starts with the fourth
approximation given in Sec. 59 of Chandrasekhar’s *Radiative Transfer*
(Dover Publications, 1960), and iteratively refines the values
according to the procedure in Sec. 60. Iteration terminates when
successive corrected values of the Y-function agree to four significant
figures.

# Inputs
- `TAU1::Float64`:  
  Normal optical thickness of the atmosphere.  
  Must be ≤ 2.0.

- `CFA::NTuple{3,Float64}`:  
  Coefficients of the characteristic function in polynomial form:  
  ```math
  C(μ) = Σⱼ Aⱼ * μ^(2(j-1)),   j = 1,2,3

Outputs

CHX::Vector{Float64}
Values of the X-function at 101 evenly spaced μ values from 0.00 to 1.00 in steps of 0.01.

CHY::Vector{Float64}
Values of the Y-function at the same μ grid.

nomitr::Int
Number of iterations performed before convergence.

Notes

If ncase != 0, a conservative case is assumed and a standard solution is returned.
The program terminates with an error if:
- tau1 > 2.0
- the characteristic function is negative for any μ
- the integral of the characteristic function exceeds 0.5

References

https://en.wikipedia.org/wiki/Chandrasekhar%27s_X-_and_Y-function

McCullough, E. C., & Porter, W. P. (1971). Computing clear day solar radiation 
spectra for the terrestrial ecological environment. Ecology, 52(6), 1008–1015.
     https://doi.org/10.2307/1933806

"""
chandrasekhar_xy(TAU1::Float64, CFA::Vector{Float64}, NCASE::Int) =
    chandrasekhar_xy!(init_chandrasekhar_xy_buffers(), TAU1, CFA, NCASE)

function chandrasekhar_xy!(buffers, TAU1::Float64, CFA::AbstractVector{Float64}, NCASE::Int)
    PSI   = buffers.PSI
    AMU   = buffers.AMU
    XA    = buffers.XA
    XB    = buffers.XB
    UMA  = @MVector zeros(5)
    ACAP = @MVector zeros(5)
    TEMX = @MVector zeros(8)
    TEMY = @MVector zeros(8)
    RTK  = @MVector zeros(5)
    ALAM = @MVector zeros(5)
    FNPP  = buffers.FNPP
    FNPN  = buffers.FNPN
    FNC0  = buffers.FNC0
    FNC1  = buffers.FNC1
    FNX   = buffers.FNX
    FNY   = buffers.FNY
    FNW   = buffers.FNW
    FMC0  = buffers.FMC0
    FMC1  = buffers.FMC1
    XD    = buffers.XD    # equivalence
    XE    = buffers.XE    # equivalence
    CHXA  = buffers.CHXA  # equivalence
    CHYA  = buffers.CHYA  # equivalence
    CHX   = buffers.CHX
    CHY   = buffers.CHY
    XA    = buffers.XA
    XB    = buffers.XB

    # Variables
    PERA = 0.0

    # Terminate if TAU1 is too large or negative
    if TAU1 <= 2.0
        # proceed
    else
        println(" THE PROGRAM IS TERMINATED BECAUSE TAU1 = ", TAU1)
        error("Program terminated due to TAU1 > 2.0")
    end

    if TAU1 < 0.0
        println(" THE PROGRAM IS TERMINATED BECAUSE TAU1 = ", TAU1)
        error("Program terminated due to TAU1 < 0.0")
    end

    PERA = CFA[1] + CFA[2] / 3.0 + 0.2 * CFA[3]
    if NCASE != 0
        PERA = 0.5
    end

    if PERA < 0.0
        println("NO COMPUTATIONS CAN BE DONE AS THE COEFFICIENTS ARE = ", CFA[1], " ", CFA[2], " ", CFA[3])
        error("Program terminated due to PERA < 0.0")
    end

    if PERA > 0.5
        println("NO COMPUTATIONS CAN BE DONE AS THE COEFFICIENTS ARE = ", CFA[1], " ", CFA[2], " ", CFA[3])
        error("Program terminated due to PERA > 0.5")
    end

    # Compute MU, PSI(MU), and weights
    for i in 1:101
        AMU[i] = (i - 1) * 0.01
        TEMA = AMU[i]^2
        PSI[i] = CFA[1] + CFA[2] * TEMA + CFA[3] * (TEMA^2)
        if PSI[i] > -1e-15
            continue
        else
            println("THE PROGRAM IS TERMINATED AS PSI($i) = ", PSI[i])
            println("NO COMPUTATIONS CAN BE DONE AS THE COEFFICIENTS ARE = ", CFA[1], " ", CFA[2], " ", CFA[3])
            error("Program terminated due to PSI(i) < threshold")
        end
    end

    XA[1] = 0.01 / 3.0
    TEMA = 4.0 * XA[1]
    TEMB = 2.0 * XA[1]
    for i in 2:2:100
        XA[i] = TEMA
        XA[i+1] = TEMB
    end
    XA[101] = XA[1]

    # Suppress all intermediate output
    NPRT = 0

    # Compute roots of the characteristic equation
    if NCASE != 0
        KMX = 5
        UMA[1] = 0.97390652851717172
        UMA[2] = 0.86506336668898451
        UMA[3] = 0.67940956829902441
        UMA[4] = 0.43339539412924719
        UMA[5] = 0.14887433898163121

        ACAP[1] = 0.066671344308688138
        ACAP[2] = 0.14945134915058059
        ACAP[3] = 0.21908636251598204
        ACAP[4] = 0.26926671930999636
        ACAP[5] = 0.29552422471475287
    else
        KMX = 4
        N1 = 0
        UMA[1] = 0.96028985649753623
        UMA[2] = 0.79666647741362674
        UMA[3] = 0.52553240991632899
        UMA[4] = 0.18343464249564980

        ACAP[1] = 0.10122853629037626
        ACAP[2] = 0.22238103445337447
        ACAP[3] = 0.31370664587788729
        ACAP[4] = 0.36268378337836198
    end

    for i in 1:KMX
        TEMX[i] = UMA[i]^2
        TEMY[i] = CFA[1] + CFA[2] * TEMX[i] + CFA[3] * TEMX[i]^2
        TEMY[i] = 2.0 * ACAP[i] * TEMY[i]
    end

    if NCASE != 0
        IST = 2
        RTK[1] = 0.0
    else
        IST = 1
    end

    for i in IST:KMX # Fortran line 152
        RTK[i] = (1.0 - TEMY[i]) / TEMX[i]
        if i == 1
            TEMA = 1.0 / UMA[1]^2
            if RTK[1] >= TEMA
                RTK[1] = 0.5 * TEMA
            end
        else
            TEMA = 1.0 / UMA[i-1]^2
            TEMB = 1.0 / UMA[i]^2
            if !(RTK[i] > TEMA && RTK[i] < TEMB)
                RTK[i] = 0.5 * (TEMA + TEMB)
            end
        end
    end

    J = IST # Fortran line 164
    while J <= KMX
        if J == 1
            TEMA = 0.0
            TEMB = 1.0 / UMA[1]^2
            N1 = 0
        else
            TEMA = 1.0 / UMA[J-1]^2
            TEMB = 1.0 / UMA[J]^2
            N1 = 0
        end

        TEMC = 1.0
        for i in 1:KMX
            TEMC -= TEMY[i] / (1.0 - RTK[J] * TEMX[i])
        end
        TEMD = abs(TEMC)
        if TEMD < 1e-14
            J += 1
        else
            N1 += 1
            if N1 > 50
                println("THE PROGRAM IS TERMINATED BECAUSE ROOTS CANNOT BE FOUND TO SATISFY THE CRITERION..")
                println("CFA(1) = $(CFA[1])    CFA(2) = $(CFA[2])    CFA(3) = $(CFA[3])")
                println("THE TROUBLE IS WITH ROOT NUMBER ", J)
                error("Root finding failed.")
            end

            if TEMC > 0.0
                TEMA = RTK[J]
            elseif TEMC < 0.0
                TEMB = RTK[J]
            end

            TEMD = 0.0
            for i in 1:KMX
                TEMD -= (TEMY[i] * TEMX[i]) / (1.0 - RTK[J] * TEMX[i])^2
            end

            TEMC = RTK[J] - TEMC / TEMD

            if TEMC <= TEMA || TEMC >= TEMB
                RTK[J] = 0.5 * (TEMA + TEMB)
            else
                RTK[J] = TEMC
            end
        end
    end

    for i in 1:KMX
        RTK[i] = sqrt(RTK[i])
    end

    if NCASE != 0
        N1 = 11
        KMX = 4
        for j in 1:KMX
            RTK[j] = RTK[j+1]
        end
    end
    UMA[1] = 0.96028985649753623
    UMA[2] = 0.79666647741362674
    UMA[3] = 0.52553240991632899
    UMA[4] = 0.18343464249564980

    ACAP[1] = 0.10122853629037626
    ACAP[2] = 0.22238103445337447
    ACAP[3] = 0.31370664587788729
    ACAP[4] = 0.36268378337836198

    # --- COMPUTE FUNCTIONS LAMDA, P AND W ---
    for j in 1:KMX
        ALAM[j] = 1.0
        for i in 1:KMX
            ALAM[j] *= (RTK[j] * UMA[i] + 1.0) / (RTK[j] * UMA[i] - 1.0)
        end
        ALAM[j] = exp(-RTK[j] * TAU1) / ALAM[j]
    end

    if NPRT != 0
        #Printf.printf("%12.5E %12.5E %12.5E\n", CFA[1], CFA[2], CFA[3])
        #Printf.printf("%12.5E\n", TAU1)
        #Printf.printf("\n")
        for j in 1:KMX
            TEMA = 1.0 / RTK[j]
            # (In the FORTRAN code, TEMA is calculated but not used or printed here)
        end
    end

    for i in 1:101 # Fortran line 225
        FNPP[i] = 1.0
        FNPN[i] = 1.0
        FNW[i] = 1.0
        for j in 1:KMX
            FNPP[i] *= (AMU[i] / UMA[j] - 1.0)
            FNPN[i] *= (-AMU[i] / UMA[j] - 1.0)
            FNW[i] *= (1.0 - RTK[j]^2 * AMU[i]^2)
        end
    end
    # --- COMPUTE C₀ AND C₁ ---

    TEMX[1] = 1.0
    TEMX[8] = 1.0
    for k in 2:7
        TEMX[k] = 1.0
        for i in 1:2
            N1 = NC0[i, k-1]
            for j in 1:2
                N2 = NC0[j+2, k-1]
                TEMX[k] *= (RTK[N1] + RTK[N2]) / (RTK[N1] - RTK[N2])
            end
        end
        TEMX[k] = -TEMX[k]
    end

    for k in 1:4
        TEMY[k] = 1.0
        N2 = NC1[4, k]
        for i in 1:3
            N1 = NC1[i, k]
            TEMY[k] *= (RTK[N1] + RTK[N2]) / (RTK[N1] - RTK[N2])
        end
    end

    for k in 5:8
        TEMY[k] = 1.0
        N1 = NC1[1, k]
        for j in 1:3
            N2 = NC1[j+1, k]
            TEMY[k] *= (RTK[N1] + RTK[N2]) / (RTK[N1] - RTK[N2])
        end
        TEMY[k] = -TEMY[k]
    end

    for i in 1:101 # Fortran line 266
        TEMA = 1.0
        TEMB = 1.0
        for j in 1:4
            TEMA *= (1.0 + RTK[j] * AMU[i])
            TEMB *= (1.0 - RTK[j] * AMU[i])
        end
        FNC0[i] = TEMA
        FMC0[i] = TEMB

        TEMA = 1.0
        TEMB = 1.0
        for j in 1:4
            TEMA *= (1.0 - RTK[j] * AMU[i]) * ALAM[j]
            TEMB *= (1.0 + RTK[j] * AMU[i]) * ALAM[j]
        end
        FNC0[i] += TEMA
        FMC0[i] += TEMB

        IST = 2
        while IST <= 7
            TEMA = 1.0
            TEMB = 1.0
            for k in 1:2
                N2 = NC0[k+2, IST-1]
                TEMA *= (1.0 - RTK[N2] * AMU[i]) * ALAM[N2]
                TEMB *= (1.0 + RTK[N2] * AMU[i]) * ALAM[N2]
            end
            for j in 1:2
                N1 = NC0[j, IST-1]
                TEMA *= (1.0 + RTK[N1] * AMU[i])
                TEMB *= (1.0 - RTK[N1] * AMU[i])
            end
            FNC0[i] += TEMA * TEMX[IST]
            FMC0[i] += TEMB * TEMX[IST]

            IST += 1
        end
    end
    for i in 1:101
        FNC1[i] = 0.0
        FMC1[i] = 0.0
        IST = 1
        while IST <= 4
            N2 = NC1[4, IST]
            TEMA = (1.0 - RTK[N2] * AMU[i]) * ALAM[N2]
            TEMB = (1.0 + RTK[N2] * AMU[i]) * ALAM[N2]
            for j in 1:3
                N1 = NC1[j, IST]
                TEMA *= (1.0 + RTK[N1] * AMU[i])
                TEMB *= (1.0 - RTK[N1] * AMU[i])
            end
            FNC1[i] += TEMY[IST] * TEMA
            FMC1[i] += TEMY[IST] * TEMB
            IST += 1
        end
        while IST <= 8
            N1 = NC1[1, IST]
            TEMA = 1.0 + RTK[N1] * AMU[i]
            TEMB = 1.0 - RTK[N1] * AMU[i]
            for j in 1:3
                N2 = NC1[j+1, IST]
                TEMA *= (1.0 - RTK[N2] * AMU[i]) * ALAM[N2]
                TEMB *= (1.0 + RTK[N2] * AMU[i]) * ALAM[N2]
            end
            FNC1[i] += TEMY[IST] * TEMA
            FMC1[i] += TEMY[IST] * TEMB
            IST += 1
        end
        FNC1[i] = -FNC1[i]
        FMC1[i] = -FMC1[i]
    end

    if NPRT != 0
        # These WRITE statements are formatted outputs; replacing with println for now.
        println()
        println("CFA values: ", CFA[1:3])
        println("TAU1 = ", TAU1)
        println("Table header for output:")

        for i in 1:101
            TEMD = FNC0[i] * FMC0[i] - FNC1[i] * FMC1[i] - (FNC0[1]^2 - FNC1[1]^2) * FNW[i]
            println(AMU[i], " ", FNPP[i], " ", FNPN[i], " ", FNW[i], " ", FNC0[i], " ", FMC0[i], " ", FNC1[i], " ", FMC1[i], " ", TEMD)
        end

        println()
    end
    # COMPUTE THE FOURTH APPROXIMATION OF X AND Y FUNCTIONS
    XB[1] = TAU1 == 0.0 ? 1.0 : 0.0 # Fortran line 345

    for i in 2:101
        XB[i] = exp(-TAU1 / AMU[i])
    end

    TEMA = 1.0 / sqrt(FNC0[1]^2 - FNC1[1]^2)
    for i in 1:101
        TEMC = TEMA / FNW[i]
        FNX[i] = (FNPN[i] * FMC0[i] - XB[i] * FNPP[i] * FNC1[i]) * TEMC
        FNY[i] = (XB[i] * FNPP[i] * FNC0[i] - FNPN[i] * FMC1[i]) * TEMC
    end

    CHXA[1] = 1.0
    CHYA[1] = XB[1]

    converged, nomitr = _chandrasekhar_xy_converge!(FNX, FNY, AMU, PSI, XA, XB, XD, XE, CHX, CHY, CHXA, CHYA)

    # if NCASE ≠ 0, generate standard solution (Fortran 975…990)
    if NCASE != 0
        tsumx = 0.0
        tsumb = 0.0
        tsumc = 0.0
        for i in 1:101
            δ = PSI[i] * AMU[i] * XA[i]
            tsumx += δ * CHX[i]
            tsumb += δ * CHY[i]
            tsumc += PSI[i] * CHY[i] * XA[i]
        end
        ratio = tsumc / (tsumx + tsumb)
        for i in 1:101
            Δ = ratio * AMU[i] * (CHX[i] + CHY[i])
            CHX[i] += Δ
            CHY[i] -= Δ
        end
    end

    return CHX, CHY, nomitr

end

# Separated out from dchxy for easier optimisation
# This algorithm is very expensive
# TODO these argument names are nightmare fuel
@noinline function _chandrasekhar_xy_converge!(FNX, FNY, AMU, PSI, XA, XB, XD, XE, CHX, CHY, CHXA, CHYA)
    nomitr = 1 # Fortran line 362
    TEMC = 0.0 # Initialize before convergence loop
    converged = false

    while !converged
        for I in 2:101
            fnx_i = FNX[I] 
            fny_i = FNY[I] 
            amu_i = AMU[I]

            #######################################################################################################
            # Compute XD and XE for this I
            # The most performance-intensive code of the package: loop inside loop inside while, called from another loop
            # Possibly there is a faster algorithm?
            # works marginally better when each line is separate
            for IC in 1:101
                XD[IC] = PSI[IC] * (fnx_i * FNX[IC] - fny_i * FNY[IC]) / (amu_i + AMU[IC])
            end
            for IC in 1:101
                XE[IC] = PSI[IC] * (fny_i * FNX[IC] - fnx_i * FNY[IC]) / (amu_i - AMU[IC])
            end
            #######################################################################################################

            # Everett's formula / interpolation for XE[I]
            XE[I] = if I <= 3
                0.5 * (XE[I+1] + XE[I-1])
            elseif I <= 5
                0.0625 * (9.0*(XE[I+1] + XE[I-1]) - XE[I+3] - XE[I-3])
            elseif I <= 96
                (3.0*(XE[I+5] + XE[I-5]) + 150.0*(XE[I+1] + XE[I-1]) - 25.0*(XE[I+3] + XE[I-3])) / 256.0
            else
                5.0*XE[I-1] + 10.0*XE[I-3] + XE[I-5] - 10.0*XE[I-2] - 5.0*XE[I-4]
            end


            #########################################################
            # Second most expensive code in the package
            # is a huge performance gain
            sxd = 0.0
            sxe = 0.0
            for ic in 1:101
                sxd += XA[ic] * XD[ic]
                sxe += XA[ic] * XE[ic]
            end
            #########################################################

            CHXA[I] = 1.0 + amu_i * sxd
            CHYA[I] = XB[I] + amu_i * sxe
        end

        # Correction to CHX and CHY
        for i in 1:101
            TEMD = TEMC * AMU[i] * (1.0 - XB[i])
            CHX[i] = CHXA[i] + TEMD
            CHY[i] = CHYA[i] + TEMD
        end

        # Check convergence (same as before)
        if nomitr > 1
            for I in 2:101
                rel_error = abs((CHY[I] - FNY[I]) / CHY[I])
                # TODO this seems wrong? shouldnt it only break if errors are <= 2.0e-4 for all I ?
                if rel_error <= 2.0e-4
                    converged = true
                    break
                end
            end
        end

        # Prepare for next iteration
        for I in 1:101
            FNX[I] = CHX[I]
            FNY[I] = CHY[I]
        end

        nomitr += 1
        nomitr > 15 && break
    end 

    return converged, nomitr
end
