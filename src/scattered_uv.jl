function allocate_scattered_uv()
    (;
        chandrasekhar_x  = zeros(101),
        chandrasekhar_y  = zeros(101),
        cosine_of_zenith_angle  = zeros(101),
        x1   = zeros(101),
        y1   = zeros(101),
        x2   = zeros(101),
        y2   = zeros(101),
        quadrature_weights  = zeros(101),
        gamma_right = zeros(101),
        gamma_left = zeros(101),
        chandrasekhar_xy_buffers = init_chandrasekhar_xy_buffers(),
    )
end

scattered_uv(optical_thickness::Float64) = scattered_uv!(allocate_scattered_uv(), optical_thickness)

function scattered_uv!(buffers, optical_thickness::Float64)
    # Large arrays (mutable, normal)
    (; chandrasekhar_x, chandrasekhar_y, cosine_of_zenith_angle, x1, y1, x2, y2, quadrature_weights, gamma_right, gamma_left, chandrasekhar_xy_buffers) = buffers

    # Small fixed-size arrays (use StaticArrays)
    intermediate_variables  = @MVector zeros(30)

    # Set up cosine_of_zenith_angle array
    cosine_of_zenith_angle[1] = 0.0
    for i in 2:101
        cosine_of_zenith_angle[i] = 0.01 * (i - 1)
    end


    # Compute x1, y1 using chandrasekhar_xy
    characteristic_function_coeffs = SVector((0.75, -0.75, 0.0))
    # chandrasekhar_x_, chandrasekhar_y_, _ = chandrasekhar_xy!(chandrasekhar_xy_buffers, optical_thickness, collect(characteristic_function_coeffs), 111)
    chandrasekhar_x_, chandrasekhar_y_, _ = chandrasekhar_xy!(chandrasekhar_xy_buffers, optical_thickness, characteristic_function_coeffs, 111)
    x1 .= chandrasekhar_x_
    y1 .= chandrasekhar_y_

    # Compute x2, y2 using chandrasekhar_xy
    characteristic_function_coeffs = SVector((0.375, -0.375, 0.0))
    # chandrasekhar_x_, chandrasekhar_y_, _ = chandrasekhar_xy!(chandrasekhar_xy_buffers, optical_thickness, collect(characteristic_function_coeffs), 0)
    chandrasekhar_x_, chandrasekhar_y_, _ = chandrasekhar_xy!(chandrasekhar_xy_buffers, optical_thickness, characteristic_function_coeffs, 0)
    x2 .= chandrasekhar_x_
    y2 .= chandrasekhar_y_

    # Compute quadrature_weights (quadrature weights)
    quadrature_weights[1] = 0.01 / 3.0
    weight_multiple_1 = 4.0 * quadrature_weights[1]
    weight_multiple_2 = 2.0 * quadrature_weights[1]
    for i in 2:2:100
        quadrature_weights[i] = weight_multiple_1
        quadrature_weights[i+1] = weight_multiple_2
    end
    quadrature_weights[101] = quadrature_weights[1]

    # Scalar accumulators
    integral_xa1 = integral_xa2 = integral_xa3 = integral_xa4 = 0.0
    integral_xb1 = integral_xb2 = integral_xb3 = integral_xb4 = integral_xb5 = integral_xb6 = integral_xb7 = integral_xb8 = 0.0

    for i in 1:101
        mu  = cosine_of_zenith_angle[i]
        mu2 = mu * mu
        mu3 = mu2 * mu

        term1 = quadrature_weights[i] * x1[i] * mu
        integral_xa1 += term1
        integral_xa2 += term1 * mu

        term2 = quadrature_weights[i] * y1[i] * mu
        integral_xa3 += term2
        integral_xa4 += term2 * mu

        term3 = quadrature_weights[i] * x2[i]
        integral_xb1 += term3
        integral_xb2 += term3 * mu
        integral_xb3 += term3 * mu2
        integral_xb4 += term3 * mu3

        term4 = quadrature_weights[i] * y2[i]
        integral_xb5 += term4
        integral_xb6 += term4 * mu
        integral_xb7 += term4 * mu2
        integral_xb8 += term4 * mu3
    end

    # Fill intermediate_variables vector
    intermediate_variables[1]  = integral_xb1 + integral_xb5 - 8.0 / 3.0
    intermediate_variables[2]  = integral_xb2 + integral_xb6
    intermediate_variables[3]  = integral_xb3 + integral_xb7
    intermediate_variables[4]  = integral_xb1 - integral_xb5 - 8.0 / 3.0
    intermediate_variables[5]  = integral_xb2 - integral_xb6
    intermediate_variables[6]  = integral_xb3 - integral_xb7
    intermediate_variables[7]  = integral_xb4 - integral_xb8
    intermediate_variables[8]  = integral_xa1 + integral_xa3
    intermediate_variables[9]  = integral_xa2 + integral_xa4
    intermediate_variables[10] = integral_xa1 - integral_xa3
    intermediate_variables[11] = integral_xa2 - integral_xa4

    intermediate_variables[12] = (intermediate_variables[1] - intermediate_variables[3]) / ((intermediate_variables[4] - intermediate_variables[6]) * optical_thickness + 2.0 * (intermediate_variables[5] - intermediate_variables[7]))
    intermediate_variables[13] = 1.0 / (intermediate_variables[4] * intermediate_variables[10] - intermediate_variables[5] * intermediate_variables[11])
    intermediate_variables[14] = 1.0 / (intermediate_variables[1] * intermediate_variables[8] - intermediate_variables[2] * intermediate_variables[9] -
                    2.0 * intermediate_variables[12] * (intermediate_variables[5] * intermediate_variables[8] - intermediate_variables[4] * intermediate_variables[9]))
    intermediate_variables[15] = 2.0 * (intermediate_variables[8] * intermediate_variables[10] - intermediate_variables[9] * intermediate_variables[11])
    intermediate_variables[16] = intermediate_variables[13] * intermediate_variables[15]
    intermediate_variables[17] = intermediate_variables[14] * intermediate_variables[15]

    term_cnu1 = 0.5 * (intermediate_variables[16] - intermediate_variables[17])
    term_cnu2 = 0.5 * (intermediate_variables[16] + intermediate_variables[17])

    intermediate_variables[15] = intermediate_variables[13] * (intermediate_variables[5] * intermediate_variables[8] - intermediate_variables[4] * intermediate_variables[9])
    intermediate_variables[16] = intermediate_variables[14] * (intermediate_variables[2] * intermediate_variables[10] - intermediate_variables[1] * intermediate_variables[11] -
                       2.0 * intermediate_variables[12] * (intermediate_variables[4] * intermediate_variables[10] - intermediate_variables[5] * intermediate_variables[11]))

    intermediate_variables[15] = intermediate_variables[13] * (intermediate_variables[2] * intermediate_variables[10] - intermediate_variables[1] * intermediate_variables[11])
    intermediate_variables[16] = intermediate_variables[14] * (intermediate_variables[5] * intermediate_variables[8] - intermediate_variables[4] * intermediate_variables[9])
    term_cu3   = 0.5 * (intermediate_variables[15] - intermediate_variables[16])
    term_cu4   = 0.5 * (intermediate_variables[15] + intermediate_variables[16])

    intermediate_variables[15] = intermediate_variables[14] * (intermediate_variables[1] * intermediate_variables[8] - intermediate_variables[2] * intermediate_variables[9])
    s_bar  = 1.0 - 0.375 * intermediate_variables[12] * (intermediate_variables[4] - intermediate_variables[6]) *
                   ((term_cnu2 - term_cnu1) * intermediate_variables[8] + (term_cu4 - term_cu3) * intermediate_variables[2] - intermediate_variables[15] * intermediate_variables[6])

    intermediate_variables[20] = 0.375 * intermediate_variables[12] * (term_cnu2 - term_cnu1) * (intermediate_variables[4] - intermediate_variables[6])
    intermediate_variables[21] = 0.375 * intermediate_variables[12] * (intermediate_variables[4] - intermediate_variables[6])
    intermediate_variables[22] = intermediate_variables[21] * (term_cu4 - term_cu3)
    intermediate_variables[23] = intermediate_variables[21] * intermediate_variables[15]

    for i in 1:101
        gamma_left[i] = intermediate_variables[20] * (x1[i] + y1[i])
        gamma_right[i] = intermediate_variables[22] * (x2[i] + y2[i]) - cosine_of_zenith_angle[i] * intermediate_variables[23] * (x2[i] - y2[i])
    end

    return gamma_right, gamma_left, s_bar
end

function init_chandrasekhar_xy_buffers()
    arrays = (;
        psi_values = zeros(101),
        cosine_of_zenith_angles = zeros(101),
        quadrature_weights_xa = zeros(101),
        quadrature_weights_xb = zeros(101),
        fn_plus_p = zeros(101),
        fn_plus_n = zeros(101),
        fn_c0 = zeros(101),
        fn_c1 = zeros(101),
        fn_x = zeros(101),
        fn_y = zeros(101),
        fn_w = zeros(101),
        fm_c0 = zeros(101),
        fm_c1 = zeros(101),
        x_d_values = zeros(101),
        x_e_values = zeros(101),
        chandrasekhar_x_approx = zeros(101),
        chandrasekhar_y_approx = zeros(101),
        chandrasekhar_x = zeros(101),
        chandrasekhar_y = zeros(101),
    )
    return arrays
end

"""
    chandrasekhar_xy(optical_thickness::Float64, characteristic_function_coeffs::Vector{Float64}, case_number::Int) 
        -> (chandrasekhar_x::Vector{Float64}, chandrasekhar_y::Vector{Float64}, num_iterations::Int)

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
- `optical_thickness::Float64`:  
  Normal optical thickness of the atmosphere.  
  Must be ≤ 2.0.

- `characteristic_function_coeffs::NTuple{3,Float64}`:  
  Coefficients of the characteristic function in polynomial form:  
  ```math
  C(μ) = Σⱼ Aⱼ * μ^(2(j-1)),   j = 1,2,3

Outputs

chandrasekhar_x::Vector{Float64}
Values of the X-function at 101 evenly spaced μ values from 0.00 to 1.00 in steps of 0.01.

chandrasekhar_y::Vector{Float64}
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
chandrasekhar_xy(optical_thickness::Float64, characteristic_function_coeffs::Vector{Float64}, case_number::Int) =
    chandrasekhar_xy!(init_chandrasekhar_xy_buffers(), optical_thickness, characteristic_function_coeffs, case_number)

function chandrasekhar_xy!(buffers, optical_thickness::Float64, characteristic_function_coeffs::AbstractVector{Float64}, case_number::Int)
    psi_values   = buffers.psi_values
    cosine_of_zenith_angles   = buffers.cosine_of_zenith_angles
    quadrature_weights_xa    = buffers.quadrature_weights_xa
    quadrature_weights_xb    = buffers.quadrature_weights_xb
    mu_roots  = @MVector zeros(5)
    cap_a_coeffs = @MVector zeros(5)
    temp_x = @MVector zeros(8)
    temp_y = @MVector zeros(8)
    k_roots  = @MVector zeros(5)
    lambda_values = @MVector zeros(5)
    fn_plus_p  = buffers.fn_plus_p
    fn_plus_n  = buffers.fn_plus_n
    fn_c0  = buffers.fn_c0
    fn_c1  = buffers.fn_c1
    fn_x   = buffers.fn_x
    fn_y   = buffers.fn_y
    fn_w   = buffers.fn_w
    fm_c0  = buffers.fm_c0
    fm_c1  = buffers.fm_c1
    x_d_values    = buffers.x_d_values    # equivalence
    x_e_values    = buffers.x_e_values    # equivalence
    chandrasekhar_x_approx  = buffers.chandrasekhar_x_approx  # equivalence
    chandrasekhar_y_approx  = buffers.chandrasekhar_y_approx  # equivalence
    chandrasekhar_x   = buffers.chandrasekhar_x
    chandrasekhar_y   = buffers.chandrasekhar_y
    quadrature_weights_xa    = buffers.quadrature_weights_xa
    quadrature_weights_xb    = buffers.quadrature_weights_xb

    # Variables
    integral_of_char_func = 0.0

    # Terminate if optical_thickness is too large or negative
    if optical_thickness <= 2.0
        # proceed
    else
        println(" THE PROGRAM IS TERMINATED BECAUSE optical_thickness = ", optical_thickness)
        error("Program terminated due to optical_thickness > 2.0")
    end

    if optical_thickness < 0.0
        println(" THE PROGRAM IS TERMINATED BECAUSE optical_thickness = ", optical_thickness)
        error("Program terminated due to optical_thickness < 0.0")
    end

    integral_of_char_func = characteristic_function_coeffs[1] + characteristic_function_coeffs[2] / 3.0 + 0.2 * characteristic_function_coeffs[3]
    if case_number != 0
        integral_of_char_func = 0.5
    end

    if integral_of_char_func < 0.0
        println("NO COMPUTATIONS CAN BE DONE AS THE COEFFICIENTS ARE = ", characteristic_function_coeffs[1], " ", characteristic_function_coeffs[2], " ", characteristic_function_coeffs[3])
        error("Program terminated due to integral_of_char_func < 0.0")
    end

    if integral_of_char_func > 0.5
        println("NO COMPUTATIONS CAN BE DONE AS THE COEFFICIENTS ARE = ", characteristic_function_coeffs[1], " ", characteristic_function_coeffs[2], " ", characteristic_function_coeffs[3])
        error("Program terminated due to integral_of_char_func > 0.5")
    end

    # Compute MU, psi_values(MU), and weights
    for i in 1:101
        cosine_of_zenith_angles[i] = (i - 1) * 0.01
        mu_sq = cosine_of_zenith_angles[i]^2
        psi_values[i] = characteristic_function_coeffs[1] + characteristic_function_coeffs[2] * mu_sq + characteristic_function_coeffs[3] * (mu_sq^2)
        if psi_values[i] > -1e-15
            continue
        else
            println("THE PROGRAM IS TERMINATED AS psi_values($i) = ", psi_values[i])
            println("NO COMPUTATIONS CAN BE DONE AS THE COEFFICIENTS ARE = ", characteristic_function_coeffs[1], " ", characteristic_function_coeffs[2], " ", characteristic_function_coeffs[3])
            error("Program terminated due to psi_values(i) < threshold")
        end
    end

    quadrature_weights_xa[1] = 0.01 / 3.0
    temp_a = 4.0 * quadrature_weights_xa[1]
    temp_b = 2.0 * quadrature_weights_xa[1]
    for i in 2:2:100
        quadrature_weights_xa[i] = temp_a
        quadrature_weights_xa[i+1] = temp_b
    end
    quadrature_weights_xa[101] = quadrature_weights_xa[1]

    # Suppress all intermediate output
    suppress_output = true

    # Compute roots of the characteristic equation
    if case_number != 0
        max_k = 5
        mu_roots[1] = 0.97390652851717172
        mu_roots[2] = 0.86506336668898451
        mu_roots[3] = 0.67940956829902441
        mu_roots[4] = 0.43339539412924719
        mu_roots[5] = 0.14887433898163121

        cap_a_coeffs[1] = 0.066671344308688138
        cap_a_coeffs[2] = 0.14945134915058059
        cap_a_coeffs[3] = 0.21908636251598204
        cap_a_coeffs[4] = 0.26926671930999636
        cap_a_coeffs[5] = 0.29552422471475287
    else
        max_k = 4
        N1 = 0
        mu_roots[1] = 0.96028985649753623
        mu_roots[2] = 0.79666647741362674
        mu_roots[3] = 0.52553240991632899
        mu_roots[4] = 0.18343464249564980

        cap_a_coeffs[1] = 0.10122853629037626
        cap_a_coeffs[2] = 0.22238103445337447
        cap_a_coeffs[3] = 0.31370664587788729
        cap_a_coeffs[4] = 0.36268378337836198
    end

    for i in 1:max_k
        temp_x[i] = mu_roots[i]^2
        temp_y[i] = characteristic_function_coeffs[1] + characteristic_function_coeffs[2] * temp_x[i] + characteristic_function_coeffs[3] * temp_x[i]^2
        temp_y[i] = 2.0 * cap_a_coeffs[i] * temp_y[i]
    end

    if case_number != 0
        start_index = 2
        k_roots[1] = 0.0
    else
        start_index = 1
    end

    for i in start_index:max_k # Fortran line 152
        k_roots[i] = (1.0 - temp_y[i]) / temp_x[i]
        if i == 1
            temp_a = 1.0 / mu_roots[1]^2
            if k_roots[1] >= temp_a
                k_roots[1] = 0.5 * temp_a
            end
        else
            temp_a = 1.0 / mu_roots[i-1]^2
            temp_b = 1.0 / mu_roots[i]^2
            if !(k_roots[i] > temp_a && k_roots[i] < temp_b)
                k_roots[i] = 0.5 * (temp_a + temp_b)
            end
        end
    end

    j_index = start_index # Fortran line 164
    while j_index <= max_k
        if j_index == 1
            temp_a = 0.0
            temp_b = 1.0 / mu_roots[1]^2
            N1 = 0
        else
            temp_a = 1.0 / mu_roots[j_index-1]^2
            temp_b = 1.0 / mu_roots[j_index]^2
            N1 = 0
        end

        temp_c = 1.0
        for i in 1:max_k
            temp_c -= temp_y[i] / (1.0 - k_roots[j_index] * temp_x[i])
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
                temp_d -= (temp_y[i] * temp_x[i]) / (1.0 - k_roots[j_index] * temp_x[i])^2
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
    mu_roots[1] = 0.96028985649753623
    mu_roots[2] = 0.79666647741362674
    mu_roots[3] = 0.52553240991632899
    mu_roots[4] = 0.18343464249564980

    cap_a_coeffs[1] = 0.10122853629037626
    cap_a_coeffs[2] = 0.22238103445337447
    cap_a_coeffs[3] = 0.31370664587788729
    cap_a_coeffs[4] = 0.36268378337836198

    # --- COMPUTE FUNCTIONS LAMDA, P AND W ---
    for j in 1:max_k
        lambda_values[j] = 1.0
        for i in 1:max_k
            lambda_values[j] *= (k_roots[j] * mu_roots[i] + 1.0) / (k_roots[j] * mu_roots[i] - 1.0)
        end
        lambda_values[j] = exp(-k_roots[j] * optical_thickness) / lambda_values[j]
    end

    if !suppress_output
        #Printf.printf("%12.5E %12.5E %12.5E\n", characteristic_function_coeffs[1], characteristic_function_coeffs[2], characteristic_function_coeffs[3])
        #Printf.printf("%12.5E\n", optical_thickness)
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
            fn_plus_p[i] *= (cosine_of_zenith_angles[i] / mu_roots[j] - 1.0)
            fn_plus_n[i] *= (-cosine_of_zenith_angles[i] / mu_roots[j] - 1.0)
            fn_w[i] *= (1.0 - k_roots[j]^2 * cosine_of_zenith_angles[i]^2)
        end
    end
    # --- COMPUTE C₀ AND C₁ ---

    temp_x[1] = 1.0
    temp_x[8] = 1.0
    for k in 2:7
        temp_x[k] = 1.0
        for i in 1:2
            N1 = NC0[i, k-1]
            for j in 1:2
                N2 = NC0[j+2, k-1]
                temp_x[k] *= (k_roots[N1] + k_roots[N2]) / (k_roots[N1] - k_roots[N2])
            end
        end
        temp_x[k] = -temp_x[k]
    end

    for k in 1:4
        temp_y[k] = 1.0
        N2 = NC1[4, k]
        for i in 1:3
            N1 = NC1[i, k]
            temp_y[k] *= (k_roots[N1] + k_roots[N2]) / (k_roots[N1] - k_roots[N2])
        end
    end

    for k in 5:8
        temp_y[k] = 1.0
        N1 = NC1[1, k]
        for j in 1:3
            N2 = NC1[j+1, k]
            temp_y[k] *= (k_roots[N1] + k_roots[N2]) / (k_roots[N1] - k_roots[N2])
        end
        temp_y[k] = -temp_y[k]
    end

    for i in 1:101 # Fortran line 266
        temp_a = 1.0
        temp_b = 1.0
        for j in 1:4
            temp_a *= (1.0 + k_roots[j] * cosine_of_zenith_angles[i])
            temp_b *= (1.0 - k_roots[j] * cosine_of_zenith_angles[i])
        end
        fn_c0[i] = temp_a
        fm_c0[i] = temp_b

        temp_a = 1.0
        temp_b = 1.0
        for j in 1:4
            temp_a *= (1.0 - k_roots[j] * cosine_of_zenith_angles[i]) * lambda_values[j]
            temp_b *= (1.0 + k_roots[j] * cosine_of_zenith_angles[i]) * lambda_values[j]
        end
        fn_c0[i] += temp_a
        fm_c0[i] += temp_b

        start_index = 2
        while start_index <= 7
            temp_a = 1.0
            temp_b = 1.0
            for k in 1:2
                N2 = NC0[k+2, start_index-1]
                temp_a *= (1.0 - k_roots[N2] * cosine_of_zenith_angles[i]) * lambda_values[N2]
                temp_b *= (1.0 + k_roots[N2] * cosine_of_zenith_angles[i]) * lambda_values[N2]
            end
            for j in 1:2
                N1 = NC0[j, start_index-1]
                temp_a *= (1.0 + k_roots[N1] * cosine_of_zenith_angles[i])
                temp_b *= (1.0 - k_roots[N1] * cosine_of_zenith_angles[i])
            end
            fn_c0[i] += temp_a * temp_x[start_index]
            fm_c0[i] += temp_b * temp_x[start_index]

            start_index += 1
        end
    end
    for i in 1:101
        fn_c1[i] = 0.0
        fm_c1[i] = 0.0
        start_index = 1
        while start_index <= 4
            N2 = NC1[4, start_index]
            temp_a = (1.0 - k_roots[N2] * cosine_of_zenith_angles[i]) * lambda_values[N2]
            temp_b = (1.0 + k_roots[N2] * cosine_of_zenith_angles[i]) * lambda_values[N2]
            for j in 1:3
                N1 = NC1[j, start_index]
                temp_a *= (1.0 + k_roots[N1] * cosine_of_zenith_angles[i])
                temp_b *= (1.0 - k_roots[N1] * cosine_of_zenith_angles[i])
            end
            fn_c1[i] += temp_y[start_index] * temp_a
            fm_c1[i] += temp_y[start_index] * temp_b
            start_index += 1
        end
        while start_index <= 8
            N1 = NC1[1, start_index]
            temp_a = 1.0 + k_roots[N1] * cosine_of_zenith_angles[i]
            temp_b = 1.0 - k_roots[N1] * cosine_of_zenith_angles[i]
            for j in 1:3
                N2 = NC1[j+1, start_index]
                temp_a *= (1.0 - k_roots[N2] * cosine_of_zenith_angles[i]) * lambda_values[N2]
                temp_b *= (1.0 + k_roots[N2] * cosine_of_zenith_angles[i]) * lambda_values[N2]
            end
            fn_c1[i] += temp_y[start_index] * temp_a
            fm_c1[i] += temp_y[start_index] * temp_b
            start_index += 1
        end
        fn_c1[i] = -fn_c1[i]
        fm_c1[i] = -fm_c1[i]
    end

    if !suppress_output
        # These WRITE statements are formatted outputs; replacing with println for now.
        println()
        println("characteristic_function_coeffs values: ", characteristic_function_coeffs[1:3])
        println("optical_thickness = ", optical_thickness)
        println("Table header for output:")

        for i in 1:101
            temp_d = fn_c0[i] * fm_c0[i] - fn_c1[i] * fm_c1[i] - (fn_c0[1]^2 - fn_c1[1]^2) * fn_w[i]
            println(cosine_of_zenith_angles[i], " ", fn_plus_p[i], " ", fn_plus_n[i], " ", fn_w[i], " ", fn_c0[i], " ", fm_c0[i], " ", fn_c1[i], " ", fm_c1[i], " ", temp_d)
        end

        println()
    end
    # COMPUTE THE FOURTH APPROXIMATION OF X AND Y FUNCTIONS
    quadrature_weights_xb[1] = optical_thickness == 0.0 ? 1.0 : 0.0 # Fortran line 345

    for i in 2:101
        quadrature_weights_xb[i] = exp(-optical_thickness / cosine_of_zenith_angles[i])
    end

    temp_a = 1.0 / sqrt(fn_c0[1]^2 - fn_c1[1]^2)
    for i in 1:101
        temp_c = temp_a / fn_w[i]
        fn_x[i] = (fn_plus_n[i] * fm_c0[i] - quadrature_weights_xb[i] * fn_plus_p[i] * fn_c1[i]) * temp_c
        fn_y[i] = (quadrature_weights_xb[i] * fn_plus_p[i] * fn_c0[i] - fn_plus_n[i] * fm_c1[i]) * temp_c
    end

    chandrasekhar_x_approx[1] = 1.0
    chandrasekhar_y_approx[1] = quadrature_weights_xb[1]

    converged, num_iterations = _chandrasekhar_xy_converge!(fn_x, fn_y, cosine_of_zenith_angles, psi_values, quadrature_weights_xa, quadrature_weights_xb, x_d_values, x_e_values, chandrasekhar_x, chandrasekhar_y, chandrasekhar_x_approx, chandrasekhar_y_approx)

    # if case_number ≠ 0, generate standard solution (Fortran 975…990)
    if case_number != 0
        tsumx = 0.0
        tsumb = 0.0
        tsumc = 0.0
        for i in 1:101
            δ = psi_values[i] * cosine_of_zenith_angles[i] * quadrature_weights_xa[i]
            tsumx += δ * chandrasekhar_x[i]
            tsumb += δ * chandrasekhar_y[i]
            tsumc += psi_values[i] * chandrasekhar_y[i] * quadrature_weights_xa[i]
        end
        ratio = tsumc / (tsumx + tsumb)
        for i in 1:101
            Δ = ratio * cosine_of_zenith_angles[i] * (chandrasekhar_x[i] + chandrasekhar_y[i])
            chandrasekhar_x[i] += Δ
            chandrasekhar_y[i] -= Δ
        end
    end

    return chandrasekhar_x, chandrasekhar_y, num_iterations

end

# Separated out from dchxy for easier optimisation
# This algorithm is very expensive
# TODO these argument names are nightmare fuel
@noinline function _chandrasekhar_xy_converge!(fn_x, fn_y, cosine_of_zenith_angles, psi_values, quadrature_weights_xa, quadrature_weights_xb, x_d_values, x_e_values, chandrasekhar_x, chandrasekhar_y, chandrasekhar_x_approx, chandrasekhar_y_approx)
    num_iterations = 1 # Fortran line 362
    temp_c = 0.0 # Initialize before convergence loop
    converged = false

    while !converged
        for i in 2:101
            fnx_i = fn_x[i] 
            fny_i = fn_y[i] 
            amu_i = cosine_of_zenith_angles[i]

            #######################################################################################################
            # Compute x_d_values and x_e_values for this i
            # The most performance-intensive code of the package: loop inside loop inside while, called from another loop
            # Possibly there is a faster algorithm?
            # works marginally better when each line is separate
            for j in 1:101
                x_d_values[j] = psi_values[j] * (fnx_i * fn_x[j] - fny_i * fn_y[j]) / (amu_i + cosine_of_zenith_angles[j])
            end
            for j in 1:101
                x_e_values[j] = psi_values[j] * (fny_i * fn_x[j] - fnx_i * fn_y[j]) / (amu_i - cosine_of_zenith_angles[j])
            end
            #######################################################################################################

            # Everett's formula / interpolation for x_e_values[i]
            x_e_values[i] = if i <= 3
                0.5 * (x_e_values[i+1] + x_e_values[i-1])
            elseif i <= 5
                0.0625 * (9.0*(x_e_values[i+1] + x_e_values[i-1]) - x_e_values[i+3] - x_e_values[i-3])
            elseif i <= 96
                (3.0*(x_e_values[i+5] + x_e_values[i-5]) + 150.0*(x_e_values[i+1] + x_e_values[i-1]) - 25.0*(x_e_values[i+3] + x_e_values[i-3])) / 256.0
            else
                5.0*x_e_values[i-1] + 10.0*x_e_values[i-3] + x_e_values[i-5] - 10.0*x_e_values[i-2] - 5.0*x_e_values[i-4]
            end


            #########################################################
            # Second most expensive code in the package
            # is a huge performance gain
            sxd = 0.0
            sxe = 0.0
            for ic in 1:101
                sxd += quadrature_weights_xa[ic] * x_d_values[ic]
                sxe += quadrature_weights_xa[ic] * x_e_values[ic]
            end
            #########################################################

            chandrasekhar_x_approx[i] = 1.0 + amu_i * sxd
            chandrasekhar_y_approx[i] = quadrature_weights_xb[i] + amu_i * sxe
        end

        # Correction to chandrasekhar_x and chandrasekhar_y
        for i in 1:101
            temp_d = temp_c * cosine_of_zenith_angles[i] * (1.0 - quadrature_weights_xb[i])
            chandrasekhar_x[i] = chandrasekhar_x_approx[i] + temp_d
            chandrasekhar_y[i] = chandrasekhar_y_approx[i] + temp_d
        end

        # Check convergence (same as before)
        if num_iterations > 1
            for i in 2:101
                rel_error = abs((chandrasekhar_y[i] - fn_y[i]) / chandrasekhar_y[i])
                # TODO this seems wrong? shouldnt it only break if errors are <= 2.0e-4 for all i ?
                if rel_error <= 2.0e-4
                    converged = true
                    break
                end
            end
        end

        # Prepare for next iteration
        for i in 1:101
            fn_x[i] = chandrasekhar_x[i]
            fn_y[i] = chandrasekhar_y[i]
        end

        num_iterations += 1
        num_iterations > 15 && break
    end 

    return converged, num_iterations
end
