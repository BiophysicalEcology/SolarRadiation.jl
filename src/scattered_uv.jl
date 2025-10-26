function allocate_scattered_radiation()
    (;
        chandrasekhar_X  = zeros(101),
        chandrasekhar_Y  = zeros(101),
        μ  = zeros(101),
        X1   = zeros(101),
        Y1   = zeros(101),
        X2   = zeros(101),
        Y2   = zeros(101),
        quad_weights  = zeros(101),
        γᵣ = zeros(101),
        γₗ = zeros(101),
        chandrasekhar_XY_buffers = init_chandrasekhar_XY_buffers(),
    )
end

scattered_radiation(τ::Float64) = scattered_radiation!(allocate_scattered_radiation(), τ)

function scattered_radiation!(buffers, τ::Float64)
    # Large arrays (mutable, normal)
    (; μ, X1, Y1, X2, Y2, quad_weights, γᵣ, γₗ, chandrasekhar_XY_buffers) = buffers

    # Small fixed-size arrays (use StaticArrays)
    AI  = @MVector zeros(30)

    # Set up μ array
    μ[1] = 0.0
    for i in 2:101
        μ[i] = 0.01 * (i - 1)
    end


    # Compute X1, Y1 using chandrasekhar_XY
    characteristic_function_coeffs = SVector((0.75, -0.75, 0.0))
    # X_, Y_, _ = chandrasekhar_xy!(chandrasekhar_XY_buffers, τ, collect(characteristic_function_coeffs), 111)
    X_, Y_, _ = chandrasekhar_xy!(chandrasekhar_XY_buffers, τ, characteristic_function_coeffs, 111)
    X1 .= X_
    Y1 .= Y_

    # Compute X2, Y2 using chandrasekhar_XY
    characteristic_function_coeffs = SVector((0.375, -0.375, 0.0))
    # X_, Y_, _ = chandrasekhar_xy!(chandrasekhar_XY_buffers, τ, collect(characteristic_function_coeffs), 0)
    X_, Y_, _ = chandrasekhar_xy!(chandrasekhar_XY_buffers, τ, characteristic_function_coeffs, 0)
    X2 .= X_
    Y2 .= Y_

    # Compute quad_weights (quadrature weights)
    quad_weights[1] = 0.01 / 3.0
    weight_multiple_1 = 4.0 * quad_weights[1]
    weight_multiple_2 = 2.0 * quad_weights[1]
    for i in 2:2:100
        quad_weights[i] = weight_multiple_1
        quad_weights[i+1] = weight_multiple_2
    end
    quad_weights[101] = quad_weights[1]

    # Scalar accumulators
    xa1 = xa2 = xa3 = xa4 = 0.0
    xb1 = xb2 = xb3 = xb4 = xb5 = xb6 = xb7 = xb8 = 0.0

    for i in 1:101
        μ_i  = μ[i]
        μ2 = μ_i * μ_i
        μ3 = μ2 * μ_i

        term1 = quad_weights[i] * X1[i] * μ_i
        xa1 += term1
        xa2 += term1 * μ_i

        term2 = quad_weights[i] * Y1[i] * μ_i
        xa3 += term2
        xa4 += term2 * μ_i

        term3 = quad_weights[i] * X2[i]
        xb1 += term3
        xb2 += term3 * μ_i
        xb3 += term3 * μ2
        xb4 += term3 * μ3

        term4 = quad_weights[i] * Y2[i]
        xb5 += term4
        xb6 += term4 * μ_i
        xb7 += term4 * μ2
        xb8 += term4 * μ3
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

    AI[12] = (AI[1] - AI[3]) / ((AI[4] - AI[6]) * τ + 2.0 * (AI[5] - AI[7]))
    AI[13] = 1.0 / (AI[4] * AI[10] - AI[5] * AI[11])
    AI[14] = 1.0 / (AI[1] * AI[8] - AI[2] * AI[9] -
                    2.0 * AI[12] * (AI[5] * AI[8] - AI[4] * AI[9]))
    AI[15] = 2.0 * (AI[8] * AI[10] - AI[9] * AI[11])
    AI[16] = AI[13] * AI[15]
    AI[17] = AI[14] * AI[15]

    term_cnu1 = 0.5 * (AI[16] - AI[17])
    term_cnu2 = 0.5 * (AI[16] + AI[17])

    AI[15] = AI[13] * (AI[5] * AI[8] - AI[4] * AI[9])
    AI[16] = AI[14] * (AI[2] * AI[10] - AI[1] * AI[11] -
                       2.0 * AI[12] * (AI[4] * AI[10] - AI[5] * AI[11]))

    AI[15] = AI[13] * (AI[2] * AI[10] - AI[1] * AI[11])
    AI[16] = AI[14] * (AI[5] * AI[8] - AI[4] * AI[9])
    term_cu3   = 0.5 * (AI[15] - AI[16])
    term_cu4   = 0.5 * (AI[15] + AI[16])

    AI[15] = AI[14] * (AI[1] * AI[8] - AI[2] * AI[9])
    s̄  = 1.0 - 0.375 * AI[12] * (AI[4] - AI[6]) *
                   ((term_cnu2 - term_cnu1) * AI[8] + (term_cu4 - term_cu3) * AI[2] - AI[15] * AI[6])

    AI[20] = 0.375 * AI[12] * (term_cnu2 - term_cnu1) * (AI[4] - AI[6])
    AI[21] = 0.375 * AI[12] * (AI[4] - AI[6])
    AI[22] = AI[21] * (term_cu4 - term_cu3)
    AI[23] = AI[21] * AI[15]

    for i in 1:101
        γₗ[i] = AI[20] * (X1[i] + Y1[i])
        γᵣ[i] = AI[22] * (X2[i] + Y2[i]) - μ[i] * AI[23] * (X2[i] - Y2[i])
    end

    return γᵣ, γₗ, s̄
end

function init_chandrasekhar_XY_buffers()
    arrays = (;
        Ψ = zeros(101),
        μ = zeros(101),
        quad_weights_Xa = zeros(101),
        quad_weights_Xb = zeros(101),
        fn_plus_p = zeros(101),
        fn_plus_n = zeros(101),
        fn_c0 = zeros(101),
        fn_c1 = zeros(101),
        fn_X = zeros(101),
        fn_Y = zeros(101),
        fn_w = zeros(101),
        fm_c0 = zeros(101),
        fm_c1 = zeros(101),
        X_d_values = zeros(101),
        X_e_values = zeros(101),
        X_approx = zeros(101),
        Y_approx = zeros(101),
        X = zeros(101),
        Y = zeros(101),
    )
    return arrays
end

"""
    chandrasekhar_XY(τ::Float64, characteristic_function_coeffs::Vector{Float64}, case_number::Int) 
        -> (X::Vector{Float64}, Y::Vector{Float64}, num_iterations::Int)

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
- `τ::Float64`:  
  Normal optical thickness of the atmosphere.  
  Must be ≤ 2.0.

- `characteristic_function_coeffs::NTuple{3,Float64}`:  
  Coefficients of the characteristic function in polynomial form:  
  ```math
  C(μ) = Σⱼ Aⱼ * μ^(2(j-1)),   j = 1,2,3

Outputs

X::Vector{Float64}
Values of the X-function at 101 evenly spaced μ values from 0.00 to 1.00 in steps of 0.01.

Y::Vector{Float64}
Values of the Y-function at the same μ grid.

num_iterations::Int
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
chandrasekhar_XY(τ::Float64, characteristic_function_coeffs::Vector{Float64}, case_number::Int) =
    chandrasekhar_xy!(init_chandrasekhar_XY_buffers(), τ, characteristic_function_coeffs, case_number)

function chandrasekhar_xy!(buffers, τ::Float64, characteristic_function_coeffs::AbstractVector{Float64}, case_number::Int)
    Ψ   = buffers.Ψ
    μ   = buffers.μ
    quad_weights_Xa    = buffers.quad_weights_Xa
    quad_weights_Xb    = buffers.quad_weights_Xb
    μ_roots  = @MVector zeros(5)
    cap_a_coeffs = @MVector zeros(5)
    temp_X = @MVector zeros(8)
    temp_Y = @MVector zeros(8)
    k_roots  = @MVector zeros(5)
    λ = @MVector zeros(5)
    fn_plus_p  = buffers.fn_plus_p
    fn_plus_n  = buffers.fn_plus_n
    fn_c0  = buffers.fn_c0
    fn_c1  = buffers.fn_c1
    fn_X   = buffers.fn_X
    fn_Y   = buffers.fn_Y
    fn_w   = buffers.fn_w
    fm_c0  = buffers.fm_c0
    fm_c1  = buffers.fm_c1
    X_d_values    = buffers.X_d_values    # equivalence
    X_e_values    = buffers.X_e_values    # equivalence
    X_approx  = buffers.X_approx  # equivalence
    Y_approx  = buffers.Y_approx  # equivalence
    X   = buffers.X
    Y   = buffers.Y
    quad_weights_Xa    = buffers.quad_weights_Xa
    quad_weights_Xb    = buffers.quad_weights_Xb

    # Variables
    integral_of_char_func = 0.0

    # Terminate if τ is too large or negative
    if τ <= 2.0
        # proceed
    else
        println(" THE PROGRAM IS TERMINATED BECAUSE τ = ", τ)
        error("Program terminated due to τ > 2.0")
    end

    if τ < 0.0
        println(" THE PROGRAM IS TERMINATED BECAUSE τ = ", τ)
        error("Program terminated due to τ < 0.0")
    end

    integral_of_char_func = characteristic_function_coeffs[1] + characteristic_function_coeffs[2] / 3.0 + 0.2 * characteristic_function_coeffs[3]
    if case_number != 0
        integral_of_char_func = 0.5
    end

    if integral_of_char_func < 0.0
        println("No computations can be done as the coefficients are = ", characteristic_function_coeffs[1], " ", characteristic_function_coeffs[2], " ", characteristic_function_coeffs[3])
        error("Program terminated due to integral_of_char_func < 0.0")
    end

    if integral_of_char_func > 0.5
        println("No computations can be done as the coefficients are = ", characteristic_function_coeffs[1], " ", characteristic_function_coeffs[2], " ", characteristic_function_coeffs[3])
        error("Program terminated due to integral_of_char_func > 0.5")
    end

    # Compute MU, Ψ(MU), and weights
    for i in 1:101
        μ[i] = (i - 1) * 0.01
        μ² = μ[i]^2
        Ψ[i] = characteristic_function_coeffs[1] + characteristic_function_coeffs[2] * μ² + characteristic_function_coeffs[3] * (μ²^2)
        if Ψ[i] > -1e-15
            continue
        else
            println("The program is terminated as Ψ($i) = ", Ψ[i])
            println("No computations can be done as the coefficients are = ", characteristic_function_coeffs[1], " ", characteristic_function_coeffs[2], " ", characteristic_function_coeffs[3])
            error("Program terminated due to Ψ(i) < threshold")
        end
    end

    quad_weights_Xa[1] = 0.01 / 3.0
    temp_a = 4.0 * quad_weights_Xa[1]
    temp_b = 2.0 * quad_weights_Xa[1]
    for i in 2:2:100
        quad_weights_Xa[i] = temp_a
        quad_weights_Xa[i+1] = temp_b
    end
    quad_weights_Xa[101] = quad_weights_Xa[1]

    # Suppress all intermediate output
    suppress_output = true

    # Compute roots of the characteristic equation
    if case_number != 0
        max_k = 5
        μ_roots[1] = 0.97390652851717172
        μ_roots[2] = 0.86506336668898451
        μ_roots[3] = 0.67940956829902441
        μ_roots[4] = 0.43339539412924719
        μ_roots[5] = 0.14887433898163121

        cap_a_coeffs[1] = 0.066671344308688138
        cap_a_coeffs[2] = 0.14945134915058059
        cap_a_coeffs[3] = 0.21908636251598204
        cap_a_coeffs[4] = 0.26926671930999636
        cap_a_coeffs[5] = 0.29552422471475287
    else
        max_k = 4
        N1 = 0
        μ_roots[1] = 0.96028985649753623
        μ_roots[2] = 0.79666647741362674
        μ_roots[3] = 0.52553240991632899
        μ_roots[4] = 0.18343464249564980

        cap_a_coeffs[1] = 0.10122853629037626
        cap_a_coeffs[2] = 0.22238103445337447
        cap_a_coeffs[3] = 0.31370664587788729
        cap_a_coeffs[4] = 0.36268378337836198
    end

    for i in 1:max_k
        temp_X[i] = μ_roots[i]^2
        temp_Y[i] = characteristic_function_coeffs[1] + characteristic_function_coeffs[2] * temp_X[i] + characteristic_function_coeffs[3] * temp_X[i]^2
        temp_Y[i] = 2.0 * cap_a_coeffs[i] * temp_Y[i]
    end

    if case_number != 0
        start_index = 2
        k_roots[1] = 0.0
    else
        start_index = 1
    end

    for i in start_index:max_k # Fortran line 152
        k_roots[i] = (1.0 - temp_Y[i]) / temp_X[i]
        if i == 1
            temp_a = 1.0 / μ_roots[1]^2
            if k_roots[1] >= temp_a
                k_roots[1] = 0.5 * temp_a
            end
        else
            temp_a = 1.0 / μ_roots[i-1]^2
            temp_b = 1.0 / μ_roots[i]^2
            if !(k_roots[i] > temp_a && k_roots[i] < temp_b)
                k_roots[i] = 0.5 * (temp_a + temp_b)
            end
        end
    end

    j_index = start_index # Fortran line 164
    while j_index <= max_k
        if j_index == 1
            temp_a = 0.0
            temp_b = 1.0 / μ_roots[1]^2
            N1 = 0
        else
            temp_a = 1.0 / μ_roots[j_index-1]^2
            temp_b = 1.0 / μ_roots[j_index]^2
            N1 = 0
        end

        temp_c = 1.0
        for i in 1:max_k
            temp_c -= temp_Y[i] / (1.0 - k_roots[j_index] * temp_X[i])
        end
        temp_d = abs(temp_c)
        if temp_d < 1e-14
            j_index += 1
        else
            N1 += 1
            if N1 > 50
                println("The program is terminated because roots cannot be found to satisfy the criterion..")
                println("cfa(1) = $(characteristic_function_coeffs[1])    cfa(2) = $(characteristic_function_coeffs[2])    cfa(3) = $(characteristic_function_coeffs[3])")
                println("The trouble is with root number ", j_index)
                error("Root finding failed.")
            end

            if temp_c > 0.0
                temp_a = k_roots[j_index]
            elseif temp_c < 0.0
                temp_b = k_roots[j_index]
            end

            temp_d = 0.0
            for i in 1:max_k
                temp_d -= (temp_Y[i] * temp_X[i]) / (1.0 - k_roots[j_index] * temp_X[i])^2
            end

            temp_c = k_roots[j_index] - temp_c / temp_d

            if temp_c <= temp_a || temp_c >= temp_b
                k_roots[j_index] = 0.5 * (temp_a + temp_b)
            else
                k_roots[j_index] = temp_c
            end
        end
    end

    for i in 1:max_k
        k_roots[i] = sqrt(k_roots[i])
    end

    if case_number != 0
        N1 = 11
        max_k = 4
        for j in 1:max_k
            k_roots[j] = k_roots[j+1]
        end
    end
    μ_roots[1] = 0.96028985649753623
    μ_roots[2] = 0.79666647741362674
    μ_roots[3] = 0.52553240991632899
    μ_roots[4] = 0.18343464249564980

    cap_a_coeffs[1] = 0.10122853629037626
    cap_a_coeffs[2] = 0.22238103445337447
    cap_a_coeffs[3] = 0.31370664587788729
    cap_a_coeffs[4] = 0.36268378337836198

    # --- COMPUTE FUNCTIONS LAMDA, P AND W ---
    for j in 1:max_k
        λ[j] = 1.0
        for i in 1:max_k
            λ[j] *= (k_roots[j] * μ_roots[i] + 1.0) / (k_roots[j] * μ_roots[i] - 1.0)
        end
        λ[j] = exp(-k_roots[j] * τ) / λ[j]
    end

    if !suppress_output 
        #Printf.printf("%12.5E %12.5E %12.5E\n", characteristic_function_coeffs[1], characteristic_function_coeffs[2], characteristic_function_coeffs[3])
        #Printf.printf("%12.5E\n", τ)
        #Printf.printf("\n")
        for j in 1:max_k
            temp_a = 1.0 / k_roots[j]
            # (In the FORTRAN code, temp_a is calculated but not used or printed here)
        end
    end

    for i in 1:101 # Fortran line 225
        fn_plus_p[i] = 1.0
        fn_plus_n[i] = 1.0
        fn_w[i] = 1.0
        for j in 1:max_k
            fn_plus_p[i] *= (μ[i] / μ_roots[j] - 1.0)
            fn_plus_n[i] *= (-μ[i] / μ_roots[j] - 1.0)
            fn_w[i] *= (1.0 - k_roots[j]^2 * μ[i]^2)
        end
    end
    # --- COMPUTE C₀ AND C₁ ---

    temp_X[1] = 1.0
    temp_X[8] = 1.0
    for k in 2:7
        temp_X[k] = 1.0
        for i in 1:2
            N1 = NC0[i, k-1]
            for j in 1:2
                N2 = NC0[j+2, k-1]
                temp_X[k] *= (k_roots[N1] + k_roots[N2]) / (k_roots[N1] - k_roots[N2])
            end
        end
        temp_X[k] = -temp_X[k]
    end

    for k in 1:4
        temp_Y[k] = 1.0
        N2 = NC1[4, k]
        for i in 1:3
            N1 = NC1[i, k]
            temp_Y[k] *= (k_roots[N1] + k_roots[N2]) / (k_roots[N1] - k_roots[N2])
        end
    end

    for k in 5:8
        temp_Y[k] = 1.0
        N1 = NC1[1, k]
        for j in 1:3
            N2 = NC1[j+1, k]
            temp_Y[k] *= (k_roots[N1] + k_roots[N2]) / (k_roots[N1] - k_roots[N2])
        end
        temp_Y[k] = -temp_Y[k]
    end

    for i in 1:101 # Fortran line 266
        temp_a = 1.0
        temp_b = 1.0
        for j in 1:4
            temp_a *= (1.0 + k_roots[j] * μ[i])
            temp_b *= (1.0 - k_roots[j] * μ[i])
        end
        fn_c0[i] = temp_a
        fm_c0[i] = temp_b

        temp_a = 1.0
        temp_b = 1.0
        for j in 1:4
            temp_a *= (1.0 - k_roots[j] * μ[i]) * λ[j]
            temp_b *= (1.0 + k_roots[j] * μ[i]) * λ[j]
        end
        fn_c0[i] += temp_a
        fm_c0[i] += temp_b

        start_index = 2
        while start_index <= 7
            temp_a = 1.0
            temp_b = 1.0
            for k in 1:2
                N2 = NC0[k+2, start_index-1]
                temp_a *= (1.0 - k_roots[N2] * μ[i]) * λ[N2]
                temp_b *= (1.0 + k_roots[N2] * μ[i]) * λ[N2]
            end
            for j in 1:2
                N1 = NC0[j, start_index-1]
                temp_a *= (1.0 + k_roots[N1] * μ[i])
                temp_b *= (1.0 - k_roots[N1] * μ[i])
            end
            fn_c0[i] += temp_a * temp_X[start_index]
            fm_c0[i] += temp_b * temp_X[start_index]

            start_index += 1
        end
    end
    for i in 1:101
        fn_c1[i] = 0.0
        fm_c1[i] = 0.0
        start_index = 1
        while start_index <= 4
            N2 = NC1[4, start_index]
            temp_a = (1.0 - k_roots[N2] * μ[i]) * λ[N2]
            temp_b = (1.0 + k_roots[N2] * μ[i]) * λ[N2]
            for j in 1:3
                N1 = NC1[j, start_index]
                temp_a *= (1.0 + k_roots[N1] * μ[i])
                temp_b *= (1.0 - k_roots[N1] * μ[i])
            end
            fn_c1[i] += temp_Y[start_index] * temp_a
            fm_c1[i] += temp_Y[start_index] * temp_b
            start_index += 1
        end
        while start_index <= 8
            N1 = NC1[1, start_index]
            temp_a = 1.0 + k_roots[N1] * μ[i]
            temp_b = 1.0 - k_roots[N1] * μ[i]
            for j in 1:3
                N2 = NC1[j+1, start_index]
                temp_a *= (1.0 - k_roots[N2] * μ[i]) * λ[N2]
                temp_b *= (1.0 + k_roots[N2] * μ[i]) * λ[N2]
            end
            fn_c1[i] += temp_Y[start_index] * temp_a
            fm_c1[i] += temp_Y[start_index] * temp_b
            start_index += 1
        end
        fn_c1[i] = -fn_c1[i]
        fm_c1[i] = -fm_c1[i]
    end

    if !suppress_output
        # These WRITE statements are formatted outputs; replacing with println for now.
        println()
        println("characteristic_function_coeffs values: ", characteristic_function_coeffs[1:3])
        println("τ = ", τ)
        println("Table header for output:")

        for i in 1:101
            temp_d = fn_c0[i] * fm_c0[i] - fn_c1[i] * fm_c1[i] - (fn_c0[1]^2 - fn_c1[1]^2) * fn_w[i]
            println(μ[i], " ", fn_plus_p[i], " ", fn_plus_n[i], " ", fn_w[i], " ", fn_c0[i], " ", fm_c0[i], " ", fn_c1[i], " ", fm_c1[i], " ", temp_d)
        end

        println()
    end
    # COMPUTE THE FOURTH APPROXIMATION OF X AND Y FUNCTIONS
    quad_weights_Xb[1] = τ == 0.0 ? 1.0 : 0.0 # Fortran line 345

    for i in 2:101
        quad_weights_Xb[i] = exp(-τ / μ[i])
    end

    temp_a = 1.0 / sqrt(fn_c0[1]^2 - fn_c1[1]^2)
    for i in 1:101
        temp_c = temp_a / fn_w[i]
        fn_X[i] = (fn_plus_n[i] * fm_c0[i] - quad_weights_Xb[i] * fn_plus_p[i] * fn_c1[i]) * temp_c
        fn_Y[i] = (quad_weights_Xb[i] * fn_plus_p[i] * fn_c0[i] - fn_plus_n[i] * fm_c1[i]) * temp_c
    end

    X_approx[1] = 1.0
    Y_approx[1] = quad_weights_Xb[1]

    converged, num_iterations = _chandrasekhar_xy_converge!(fn_X, fn_Y, μ, Ψ, quad_weights_Xa, quad_weights_Xb, X_d_values, X_e_values, X, Y, X_approx, Y_approx)

    # if case_number ≠ 0, generate standard solution (Fortran 975…990)
    if case_number != 0
        tsumx = 0.0
        tsumb = 0.0
        tsumc = 0.0
        for i in 1:101
            δ = Ψ[i] * μ[i] * quad_weights_Xa[i]
            tsumx += δ * X[i]
            tsumb += δ * Y[i]
            tsumc += Ψ[i] * Y[i] * quad_weights_Xa[i]
        end
        ratio = tsumc / (tsumx + tsumb)
        for i in 1:101
            Δ = ratio * μ[i] * (X[i] + Y[i])
            X[i] += Δ
            Y[i] -= Δ
        end
    end

    return X, Y, num_iterations

end

# Separated out from dchxy for easier optimisation
# This algorithm is very expensive
# TODO these argument names are nightmare fuel
@noinline function _chandrasekhar_xy_converge!(fn_X, fn_Y, μ, Ψ, quad_weights_Xa, quad_weights_Xb, X_d_values, X_e_values, X, Y, X_approx, Y_approx)
    num_iterations = 1 # Fortran line 362
    temp_c = 0.0 # Initialize before convergence loop
    converged = false

    while !converged
        for i in 2:101
            fnx_i = fn_X[i] 
            fny_i = fn_Y[i] 
            amu_i = μ[i]

            #######################################################################################################
            # Compute X_d_values and X_e_values for this i
            # The most performance-intensive code of the package: loop inside loop inside while, called from another loop
            # Possibly there is a faster algorithm?
            # works marginally better when each line is separate
            for j in 1:101
                X_d_values[j] = Ψ[j] * (fnx_i * fn_X[j] - fny_i * fn_Y[j]) / (amu_i + μ[j])
            end
            for j in 1:101
                X_e_values[j] = Ψ[j] * (fny_i * fn_X[j] - fnx_i * fn_Y[j]) / (amu_i - μ[j])
            end
            #######################################################################################################

            # Everett's formula / interpolation for X_e_values[i]
            X_e_values[i] = if i <= 3
                0.5 * (X_e_values[i+1] + X_e_values[i-1])
            elseif i <= 5
                0.0625 * (9.0*(X_e_values[i+1] + X_e_values[i-1]) - X_e_values[i+3] - X_e_values[i-3])
            elseif i <= 96
                (3.0*(X_e_values[i+5] + X_e_values[i-5]) + 150.0*(X_e_values[i+1] + X_e_values[i-1]) - 25.0*(X_e_values[i+3] + X_e_values[i-3])) / 256.0
            else
                5.0*X_e_values[i-1] + 10.0*X_e_values[i-3] + X_e_values[i-5] - 10.0*X_e_values[i-2] - 5.0*X_e_values[i-4]
            end


            #########################################################
            # Second most expensive code in the package
            # is a huge performance gain
            sxd = 0.0
            sxe = 0.0
            for ic in 1:101
                sxd += quad_weights_Xa[ic] * X_d_values[ic]
                sxe += quad_weights_Xa[ic] * X_e_values[ic]
            end
            #########################################################

            X_approx[i] = 1.0 + amu_i * sxd
            Y_approx[i] = quad_weights_Xb[i] + amu_i * sxe
        end

        # Correction to X and Y
        for i in 1:101
            temp_d = temp_c * μ[i] * (1.0 - quad_weights_Xb[i])
            X[i] = X_approx[i] + temp_d
            Y[i] = Y_approx[i] + temp_d
        end

        # Check convergence (same as before)
        if num_iterations > 1
            for i in 2:101
                rel_error = abs((Y[i] - fn_Y[i]) / Y[i])
                # TODO this seems wrong? shouldnt it only break if errors are <= 2.0e-4 for all i ?
                if rel_error <= 2.0e-4
                    converged = true
                    break
                end
            end
        end

        # Prepare for next iteration
        for i in 1:101
            fn_X[i] = X[i]
            fn_Y[i] = Y[i]
        end

        num_iterations += 1
        num_iterations > 15 && break
    end 

    return converged, num_iterations
end
